"""
This pipeline runs TR genotyping tools on the catalog(s) of interest.

Per-sample inputs:
    short read or long read cram/crai file
    filter_vcf output catalogs for the tools of interest
    output directory
"""

import collections
import hailtop.fs as hfs
import os
import pandas as pd
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize
import sys

sys.path.append("../str-truth-set/tool_comparison/hail_batch_pipelines")
from expansion_hunter_pipeline import create_expansion_hunter_steps, DOCKER_IMAGE as EH_DOCKER_IMAGE
from gangstr_pipeline import create_gangstr_steps
from hipstr_pipeline import create_hipstr_steps
from constrain_pipeline import create_constrain_step
from trgt_pipeline import create_trgt_step, DOCKER_IMAGE as TRGT_V5_DOCKER_IMAGE, TRGT_V3_DOCKER_IMAGE
from longtr_pipeline import create_longtr_steps
from inquistr_pipeline import create_inquistr_steps
from vamos_pipeline import create_vamos_step


SHORT_READ_TOOLS = {
    "IlluminaEHv5",
    "EHv5",
    "EHv5-bw2-optimized",
    "GangSTR",
    "HipSTR",
    "constrain"
}

LONG_READ_TOOLS = {
    "TRGTv3",
    "TRGTv5",
    "LongTR",
    "inquiSTR",
    "vamos",
}

# The add-columns and plot steps use the /str-truth-set baked into FILTER_VCFS_DOCKER_IMAGE (rebuild that image and
# update its digest to pick up new str-truth-set scripts).

# Motif size bins (min, max) used to stratify the accuracy plots.
MOTIF_SIZE_BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (2, 6), (7, 24), (25, 1000)]

# EHv5/EHv5-bw2-optimized stream single-threaded; split the catalog into this many parallel 1-cpu jobs
# (bw2-fork --start-with/--n-loci) to cut wall time ~N-fold at ~constant total cost.
EHV5_NUM_SHARDS = 10

SHORT_READ_DATA_TYPES = {
    "illumina",
    "illumina_exome",
    "element",
    "ultima",
    "illumina_rnaseq",
}

LONG_READ_DATA_TYPES = {
    "pacbio",
    "ONT",
    "pacbio_isoseq",
}

# RNA-seq data types are labeled by total bases sequenced (Gbp) rather than genome coverage, since
# genome-wide depth is meaningless for transcript data (reads only cover expressed loci).
RNASEQ_DATA_TYPES = {
    "illumina_rnaseq",
    "pacbio_isoseq",
}

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

FILTER_VCFS_DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:251400db5cc8837029e8ecfcb6bdd525f6a7e84c03516769dd05e8e58fe78cc8"

DEFAULT_OUTPUT_DIR = "gs://str-truth-set-v2/tool_results"

def main():
    sample_table_path = "HPRC_all_aligned_short_read_and_long_read_samples.tsv"
    df = pd.read_table(sample_table_path)

    bp = pipeline("run_genotyping_tools", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("-s", "--sample-id", action="append",
                        help="Process only this sample. Can be specified more than once.")
    parser.add_argument("-t", "--tool", action="append", choices=SHORT_READ_TOOLS|LONG_READ_TOOLS, help="The tool to run.")
    parser.add_argument("--data-type", action="append", choices=SHORT_READ_DATA_TYPES|LONG_READ_DATA_TYPES, help="Which data type(s) to process")
    parser.add_argument("-k", "--filename-keyword", help="If specified, only BAM paths that contain this keyword will be processed", action="append")
    parser.add_argument("--filter-vcf-dir", default="gs://str-truth-set-v2/filter_vcf", help="Base dir for filter_vcf pipeline output files")
    parser.add_argument("--truth-set-genotypes-dir", default="gs://str-truth-set-v2/filter_vcf_v2",
                        help="Base dir for the filter_vcf_to_tandem_repeats genotype step output "
                             "({sample_id}/{sample_id}.tandem_repeat_genotypes.tsv.gz), used as the truth set "
                             "(it carries the per-allele repeat purity used by the purity-stratified plots)")
    parser.add_argument("--custom-catalog-path", help="If specified, use this catalog instead of the filter_vcf catalogs")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    args = bp.parse_known_args()

    if not args.tool:
        args.tool = ["TRGTv5"]
    if not args.data_type:
        args.data_type = ["pacbio"]

    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]

    if args.data_type:
        df = df[df.sequencing_data_type.isin(args.data_type)]

    df = df[df.sample_id == "HG002"]  # only HG002 is used for tool evaluations (CHM1_CHM13 excluded)

    if args.custom_catalog_path and args.output_dir == DEFAULT_OUTPUT_DIR:
        parser.error("--custom-catalog-path is set without also setting --output-dir")

    bp.precache_file_paths(os.path.join(args.output_dir, "**/*.*"))
    # precache the prefiltered IlluminaEHv5 catalog(s) so the prefilter step is skipped once it already exists
    bp.precache_file_paths(os.path.join(args.filter_vcf_dir, "**/*.without_loci_with_flanking_Ns.json"))


    download_to_dir = "../results"

    # IlluminaEHv5 prefilter steps, keyed by source EHv5 catalog path so the catalog is filtered once and the step
    # is reused across all IlluminaEHv5 data types/coverages; values are (step, filtered_catalog_path) tuples.
    illumina_eh_prefilter_steps = {}

    for row_i, (_, row) in enumerate(df.iterrows()):
        if args.filename_keyword:
            if not any(keyword in row.read_data_path for keyword in args.filename_keyword):
                continue

        coverage = int(round(float(row.depth_of_coverage)))
        # RNA-seq rows carry total bases sequenced (Gbp) rather than genome coverage, so label them
        # "{N}G" (e.g. "24G") instead of "{N}x" in the output dir, plots, and the Coverage column.
        coverage_label = f"{coverage}G" if row.sequencing_data_type in RNASEQ_DATA_TYPES else f"{coverage}x"
        for tool in args.tool:
            if tool in SHORT_READ_TOOLS and row.sequencing_data_type not in SHORT_READ_DATA_TYPES:
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                      f"doesn't support {row.sequencing_data_type} data")
                continue
            if tool in LONG_READ_TOOLS and row.sequencing_data_type not in LONG_READ_DATA_TYPES:
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                      f"doesn't support {row.sequencing_data_type} data")
                continue

            # GangSTR and IlluminaEHv5 are only run on illumina and illumina_exome data (GangSTR never completes on
            # ultima, and the original Illumina ExpansionHunter build is only meaningful on Illumina WGS/exome data).
            if tool in ("GangSTR", "IlluminaEHv5") and row.sequencing_data_type not in ("illumina", "illumina_exome"):
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} "
                      f"({tool} is only run on illumina and illumina_exome data)")
                continue

            # vamos is not run on pacbio_isoseq data.
            if tool == "vamos" and row.sequencing_data_type == "pacbio_isoseq":
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} "
                      f"(vamos is not run on pacbio_isoseq data)")
                continue

            if args.custom_catalog_path:
                repeat_catalog_paths = args.custom_catalog_path
            elif tool == "inquiSTR":
                # inquiSTR genotypes from a plain region bed, so use the {sample_id}.bed.gz loci catalog
                # (columns: chrom, start0, end, motif). The inquiSTR step derives its region bed from it.
                repeat_catalog_paths = os.path.join(args.filter_vcf_dir, row.sample_id,
                    f"{row.sample_id}.bed.gz")
            elif tool == "vamos":
                # vamos derives its catalog (inside the step) from the unsharded ExpansionHunter catalog json
                repeat_catalog_paths = os.path.join(args.filter_vcf_dir, row.sample_id,
                    f"{row.sample_id}.EHv5*")
            elif tool in ("EHv5", "EHv5-bw2-optimized", "IlluminaEHv5"):
                # all three ExpansionHunter v5 variants share the EHv5 catalog set; the branch below picks the
                # unsharded catalog for EHv5/EHv5-bw2-optimized and the sharded catalog(s) for IlluminaEHv5
                repeat_catalog_paths = os.path.join(args.filter_vcf_dir, row.sample_id,
                    f"{row.sample_id}.EHv5*")
            elif tool in ("TRGTv3", "TRGTv5"):
                # both TRGT versions read the same TRGT BED catalog ({sample_id}.TRGT*)
                repeat_catalog_paths = os.path.join(args.filter_vcf_dir, row.sample_id,
                    f"{row.sample_id}.TRGT*")
            else:
                catalog_path_suffix = tool

                repeat_catalog_paths = os.path.join(args.filter_vcf_dir, row.sample_id,
                    f"{row.sample_id}.{catalog_path_suffix}*")

            print(f"Listing catalogs {repeat_catalog_paths}")
            repeat_catalog_paths = [x.path for x in hfs.ls(repeat_catalog_paths)]
            output_dir = os.path.join(args.output_dir, row.sample_id, row.sequencing_data_type, tool, f"{coverage_label}_coverage")
            if tool in ("EHv5", "EHv5-bw2-optimized", "IlluminaEHv5"):
                # Three ExpansionHunter v5 variants, all genotyped with create_expansion_hunter_steps:
                #   EHv5               - bw2 fork, --analysis-mode low-mem-streaming
                #   EHv5-bw2-optimized - bw2 fork, --analysis-mode optimized-streaming --improved-genotyping
                #   IlluminaEHv5       - original Illumina build, --analysis-mode streaming
                use_illumina_expansion_hunter = (tool == "IlluminaEHv5")
                if tool == "EHv5":
                    analysis_mode = "low-mem-streaming"
                elif tool == "EHv5-bw2-optimized":
                    # optimized-streaming implies --improved-genotyping inside create_expansion_hunter_steps
                    analysis_mode = "optimized-streaming"
                else:  # IlluminaEHv5
                    analysis_mode = "streaming"

                catalog_prefilter_step = None
                if args.custom_catalog_path:
                    variant_catalog_file_paths = repeat_catalog_paths
                elif tool == "IlluminaEHv5":
                    # The official Illumina ExpansionHunter build aborts on the first catalog locus with >5 Ns in its
                    # +/-1000bp flanks ("Flanks can contain at most 5 characters N but found x Ns"); the bw2 fork
                    # tolerates them. So prefilter the single unsharded EHv5 catalog to drop those loci before
                    # genotyping. The prefilter step is built once per source catalog (cached in
                    # illumina_eh_prefilter_steps) and reused across all IlluminaEHv5 data types/coverages.
                    source_eh_catalog_path = next(p for p in repeat_catalog_paths if "001_of_001" in p)
                    if source_eh_catalog_path not in illumina_eh_prefilter_steps:
                        illumina_eh_prefilter_steps[source_eh_catalog_path] = create_illumina_eh_catalog_prefilter_step(
                            bp,
                            eh_catalog_path=source_eh_catalog_path,
                            reference_fasta=REFERENCE_FASTA_PATH,
                            reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                            output_dir=os.path.join(args.filter_vcf_dir, row.sample_id))
                    catalog_prefilter_step, filtered_catalog_path = illumina_eh_prefilter_steps[source_eh_catalog_path]
                    variant_catalog_file_paths = [filtered_catalog_path]
                else:
                    # the streaming EHv5 / EHv5-bw2-optimized variants use the single unsharded catalog
                    variant_catalog_file_paths = [p for p in repeat_catalog_paths if "001_of_001" in p]

                current_step = create_expansion_hunter_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    variant_catalog_file_paths=variant_catalog_file_paths,
                    output_dir=output_dir,
                    output_prefix= f"{row.sample_id}.{tool}",
                    analysis_mode=analysis_mode,
                    loci_to_exclude=None,
                    min_locus_coverage=None,
                    use_illumina_expansion_hunter=use_illumina_expansion_hunter,
                    # IlluminaEHv5 genotypes the single prefiltered catalog and must wait for the prefilter step
                    catalog_prefilter_step=catalog_prefilter_step,
                    # EHv5/EHv5-bw2-optimized stream single-threaded, so split the catalog into EHV5_NUM_SHARDS
                    # parallel 1-cpu jobs (no-op for IlluminaEHv5, which genotypes one prefiltered catalog)
                    num_shards=EHV5_NUM_SHARDS)
            elif tool == "GangSTR":
                current_step = create_gangstr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    repeat_spec_file_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}")
            elif tool == "HipSTR":
                current_step = create_hipstr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    regions_bed_file_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}")
            elif tool == "constrain":
                current_step = create_constrain_step(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    constrain_catalog_bed_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}",
                    cpu=1,
                )
            elif tool in ("TRGTv3", "TRGTv5"):
                current_step = create_trgt_step(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    trgt_catalog_bed_paths=repeat_catalog_paths,
                    parse_reference_region_from_locus_id=True,
                    output_dir=output_dir,
                    output_prefix= f"{row.sample_id}.{tool}",
                    docker_image=TRGT_V3_DOCKER_IMAGE if tool == "TRGTv3" else TRGT_V5_DOCKER_IMAGE)
            elif tool == "LongTR":
                current_step = create_longtr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    regions_bed_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}")
            elif tool == "inquiSTR":
                current_step = create_inquistr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    inquistr_catalog_bed_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}")
            elif tool == "vamos":
                # use the unsharded ExpansionHunter catalog json; the vamos step converts it to a vamos catalog
                current_step = create_vamos_step(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    # vamos needs the single unsharded EHv5 catalog; pick it out of the filter_vcf catalogs by the
                    # "001_of_001" shard name, but pass a --custom-catalog-path through unfiltered (its filename
                    # won't contain that token, so filtering would leave an empty list and crash the vamos step).
                    expansion_hunter_catalog_paths=(
                        repeat_catalog_paths if args.custom_catalog_path
                        else [p for p in repeat_catalog_paths if "001_of_001" in p]),
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}")
            else:
                raise ValueError(f"Unknown tool: {tool}")


            # EHv5, EHv5-bw2-optimized, and IlluminaEHv5 each keep their own label downstream (all three are
            # registered in add_tool_results_columns.py / add_concordance_columns.py / plot_tool_accuracy_by_allele_size.py).
            add_columns_step = add_tool_comparison_columns_step(
                bp,
                current_step,
                tool=tool,
                coverage_label=coverage_label,
                sample_id=row.sample_id,
                output_dir=output_dir,
                truth_set_genotypes_path=os.path.join(
                    args.truth_set_genotypes_dir, row.sample_id, f"{row.sample_id}.tandem_repeat_genotypes.tsv.gz"),
                tool2="Truth",
                download_to_dir=download_to_dir)

            plot_tool_accuracy_step = create_plot_tool_accuracy_steps(
                bp,
                add_columns_step,
                tool=tool,
                coverage_label=coverage_label,
                sample_id=row.sample_id,
                output_dir=output_dir)
    bp.run()


def create_illumina_eh_catalog_prefilter_step(bp, *, eh_catalog_path, reference_fasta, reference_fasta_fai, output_dir):
    """Build a step that drops loci with >5 Ns in their flanks from an ExpansionHunter catalog, for IlluminaEHv5.

    The official Illumina ExpansionHunter v5 build aborts on the first catalog locus with more than 5 Ns in its
    +/-1000bp flanks ("Flanks can contain at most 5 characters N but found x Ns"), while the bw2 fork tolerates them.
    str_analysis.filter_out_loci_with_Ns_in_flanks scans each locus's +/-1000bp reference flanks (the EH default
    --region-extension-length) and writes a catalog with the offending loci removed. The official build can't read a
    gzipped catalog, so the filtered catalog is written as plain .json.

    Args:
        bp: the step_pipeline pipeline object.
        eh_catalog_path: gs:// path of the source ExpansionHunter (EHv5) variant catalog json.
        reference_fasta: gs:// path of the reference fasta (used to read each locus's flanking sequence).
        reference_fasta_fai: gs:// path of the reference fasta .fai index.
        output_dir: directory to write the filtered catalog into (the filter_vcf dir, so it is reused across runs).

    Returns:
        A (step, filtered_catalog_path) tuple. filtered_catalog_path is the plain-json catalog the IlluminaEHv5
        genotyping step should read; the genotyping step must depend on the returned step.
    """
    filtered_catalog_filename = re.sub(r"\.json(\.gz)?$", "", os.path.basename(eh_catalog_path)) + \
        ".without_loci_with_flanking_Ns.json"
    filtered_loci_filename = filtered_catalog_filename.replace(".json", ".filtered_loci.txt")
    filtered_catalog_path = os.path.join(output_dir, filtered_catalog_filename)

    # the EH image (weisburd/str-analysis-with-expansion-hunter, pinned by digest) bakes in str_analysis, so the
    # filter runs reproducibly without a runtime pip install
    step = bp.new_step(
        name=f"Prefilter EHv5 catalog (drop flanking-N loci) for IlluminaEHv5: {os.path.basename(eh_catalog_path)}",
        arg_suffix="prefilter-illumina-eh-catalog-step",
        image=EH_DOCKER_IMAGE,
        cpu=2,
        memory="standard",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)
    step.command("set -ex")
    local_fasta = step.input(reference_fasta)
    step.input(reference_fasta_fai)
    local_catalog = step.input(eh_catalog_path)
    step.command(
        f"python3 -m str_analysis.filter_out_loci_with_Ns_in_flanks "
        f"-R {local_fasta} --region-extension-length 1000 "
        f"-o {filtered_catalog_filename} "
        f"-f {filtered_loci_filename} "
        f"{local_catalog}")
    step.command("ls -lhrt")
    step.output(filtered_catalog_filename)
    step.output(filtered_loci_filename)
    return step, filtered_catalog_path


def add_tool_comparison_columns_step(bp, tool_results_step, *, tool, coverage_label, sample_id, output_dir, truth_set_genotypes_path, tool2="Truth", download_to_dir=None):
    tool_results_path = None
    for output_spec in tool_results_step.get_outputs():
        if output_spec.output_path.endswith(".variants.tsv.gz"):
            tool_results_path = output_spec.output_path

    if tool_results_path is None:
        raise ValueError(f"Couldn't find a .variants.tsv.gz table among the output files for {tool}: "
                         f"{[s.output_path for s in tool_results_step.get_outputs()]}")

    add_columns_step = bp.new_step(
        name=f"Add {sample_id} {tool} results columns to {tool2} table for {os.path.basename(output_dir)}",
        arg_suffix=f"add-columns-step",
        image=FILTER_VCFS_DOCKER_IMAGE,
        cpu=2,
        memory="highmem",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)

    add_columns_step.depends_on(tool_results_step)

    add_columns_step.command("set -ex")

    local_tool_results_input = add_columns_step.input(tool_results_path)
    # the truth set is the v2 genotype table; compute_truth_set_tsv_for_comparisons.py normalizes its columns and
    # carries the per-allele repeat purity (RepeatPurity: Allele 1/2) used by the purity-stratified plots
    local_truth_set_genotypes = add_columns_step.input(truth_set_genotypes_path)

    add_columns_step.command(f"""python3 <<EOF
import pandas as pd
print("Adding columns to {local_tool_results_input}")
df = pd.read_table("{local_tool_results_input}", dtype=str)
df.loc[:, "Coverage"] = "{coverage_label}"
print(f"Writing {{len(df):,d}} rows to {local_tool_results_input}")
df.to_csv("{local_tool_results_input}", sep="\\t", index=False, header=True)
EOF
""")

    # matches the output name compute_truth_set_tsv_for_comparisons.py derives from the input basename
    for_comparison_filename = os.path.basename(truth_set_genotypes_path).replace(".tsv", ".for_comparison.tsv")
    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py "
               f"--output-dir . "
               f"{local_truth_set_genotypes}")

    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/add_tool_results_columns.py "
               f"--tool {tool} "
               f"{local_tool_results_input} "
               f"{for_comparison_filename} ")

    add_columns_step.command("ls -lhrt")

    local_tsv_file_path = for_comparison_filename.replace(".tsv.gz", "") + f".with_{tool}_results.tsv.gz"
    output_filename = for_comparison_filename.replace(".tsv.gz", "") + f".with_{tool}_vs_{tool2}_columns.tsv.gz"
    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/add_concordance_columns.py "
               f"--tool {tool} "
               f"--compare-to {tool2} "
               f"--output-tsv {output_filename} "
               f"{local_tsv_file_path}")

    add_columns_step.command("ls -lhrt")

    add_columns_step.output(output_filename, download_to_dir=download_to_dir)
    add_columns_step.output(output_filename.replace(".tsv", ".alleles.tsv"))

    return add_columns_step


def create_plot_tool_accuracy_steps(bp, add_columns_step, *, tool, coverage_label, sample_id, output_dir):
    plot_tool_accuracy_step = bp.new_step(
        name=f"Plot {sample_id} {tool} accuracy for {os.path.basename(output_dir)}",
        arg_suffix=f"plot-accuracy-step",
        image=FILTER_VCFS_DOCKER_IMAGE,
        cpu=1,
        memory="highmem",
        storage="20Gi",
        output_dir=output_dir)

    local_variants_tsv, local_alleles_tsv = plot_tool_accuracy_step.use_previous_step_outputs_as_inputs(add_columns_step)

    plot_tool_accuracy_step.command("set -ex")

    # The plot script stratifies internally by purity bin, IsPureRepeat, and no-call loci, so each invocation produces
    # many svg files. They're all captured below with a wildcard output.
    for min_motif_size, max_motif_size in MOTIF_SIZE_BINS:
        plot_tool_accuracy_step.command(
            f"python3 -u /str-truth-set/figures_and_tables/plot_tool_accuracy_by_allele_size.py "
            "--verbose "
            f"--tool {tool} "
            f"--coverage {coverage_label} "
            "--q-threshold 0 "
            f"--min-motif-size {min_motif_size} "
            f"--max-motif-size {max_motif_size} "
            "--image-type svg "
            "--show-title "
            f"{local_alleles_tsv} ")
        plot_tool_accuracy_step.command("ls -lhrt")

    # Also generate the unstratified "all motif sizes" plot (".all_motifs", surfaced as the "all" bin in the viewer)
    # by running the script in its default mode with no --min/--max-motif-size. That mode also emits a few default-bin
    # svgs (.2bp_motifs/.3to6bp_motifs/.7to24bp_motifs/.25to50bp_motifs) whose tokens the viewer doesn't use — harmless
    # extras captured by the same wildcard output below.
    plot_tool_accuracy_step.command(
        f"python3 -u /str-truth-set/figures_and_tables/plot_tool_accuracy_by_allele_size.py "
        "--verbose "
        f"--tool {tool} "
        f"--coverage {coverage_label} "
        "--q-threshold 0 "
        "--image-type svg "
        "--show-title "
        f"{local_alleles_tsv} ")
    plot_tool_accuracy_step.command("ls -lhrt")

    # gzip each svg in place (keeping the .svg name) and serve it with the right content headers. The plots are
    # uploaded with a wildcard via gcloud storage cp (Delocalize.GSUTIL_COPY), since Delocalize.COPY needs an explicit
    # filename per output.
    plot_tool_accuracy_step.command('for f in tool_accuracy_by_true_allele_size.*.svg; do gzip "$f"; mv "$f.gz" "$f"; done')
    plot_tool_accuracy_step.output("tool_accuracy_by_true_allele_size.*.svg", delocalize_by=Delocalize.GSUTIL_COPY)

    image_headers_step = bp.new_step(
        name=f"Set image headers for {sample_id} {tool} accuracy plots",
        arg_suffix="image-headers-step",
        image=FILTER_VCFS_DOCKER_IMAGE,
        cpu=1,
        output_dir=output_dir)

    image_headers_step.depends_on(plot_tool_accuracy_step)

    image_headers_step.command("set -ex")
    # gcloud storage objects update on a freshly-uploaded set of objects can intermittently return
    # "HTTPError 409 ... edited during the operation", so retry a few times.
    svg_glob = os.path.join(output_dir, 'tool_accuracy_by_true_allele_size.*.svg')
    image_headers_step.command(
        f"for attempt in 1 2 3 4 5; do "
        f"if gcloud storage objects update --content-type 'image/svg+xml' --content-encoding 'gzip' {svg_glob}; "
        f"then break; fi; "
        f"echo \"image headers update attempt $attempt failed; retrying\"; sleep 15; "
        f"done")

    #plot_tool_accuracy_step.command(f"python3 -u /str-truth-set/figures_and_tables/plot_tool_accuracy_vs_Q.py "
    #                                "--verbose "
    #                                f"--coverage {coverage_label} "
    #                                f"--min-motif-size {min_motif_size} "
    #                                f"--max-motif-size {max_motif_size} "
    #                                "--genotype all "
    #                                "--show-no-call-loci "
    #                                "--image-type svg "
    #                                "--show-title "
    #                                f"{local_alleles_tsv} ")

    plot_tool_accuracy_step.command("ls -lhrt")

    return plot_tool_accuracy_step

#
#    python3 -u plot_tool_accuracy_percent_exactly_right.py --show-title --output-dir ${output_dir}
#    python3 -u plot_tool_accuracy_percent_exactly_right.py --show-title --exclude-hipstr-no-call-loci --output-dir ${output_dir}
#    python3 -u plot_tool_accuracy_by_motif_size.py --output-dir ${output_dir}
#
#    # generate tool accuracy vs Q and tool accuracy by num repeats plots
#    python3 figures_pipeline.py --force --batch-size 25
#
#    gsutil -m cp -r gs://str-truth-set/hg38/figures/accuracy_vs_Q .
#    gsutil -m cp -r gs://str-truth-set/hg38/figures/accuracy_by_allele_size  .
#
#    #python3 -u plot_tool_accuracy_vs_Q.py --verbose
#    #python3 -u plot_tool_accuracy_by_allele_size.py --verbose
#
#    ./generate_figure_panels_for_paper.sh
#
#    TODO plot fraction of loci that are polymorphic by size in reference

if __name__ == "__main__":
    main()
