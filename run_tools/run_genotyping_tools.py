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
from step_pipeline import pipeline, Backend, Localize, Delocalize
import sys

sys.path.append("../str-truth-set/tool_comparison/hail_batch_pipelines")
from expansion_hunter_pipeline import create_expansion_hunter_steps, create_expansion_hunter_dev_steps
from gangstr_pipeline import create_gangstr_steps
from hipstr_pipeline import create_hipstr_steps
from constrain_pipeline import create_constrain_step
from trgt_pipeline import create_trgt_step
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
    "TRGT",
    "LongTR",
    "inquiSTR",
    "vamos",
}

# The add-columns and plot steps use the /str-truth-set baked into FILTER_VCFS_DOCKER_IMAGE (rebuild that image and
# update its digest to pick up new str-truth-set scripts).

# Motif size bins (min, max) used to stratify the accuracy plots.
MOTIF_SIZE_BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (2, 6), (7, 24), (25, 1000)]

SHORT_READ_DATA_TYPES = {
    "illumina",
    "illumina_exome",
    "element",
    "ultima",
}

LONG_READ_DATA_TYPES = {
    "pacbio",
    "ONT",
}

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

FILTER_VCFS_DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:ceea479fcadac72813986411be2fc50549a4e80e3a47f97d431e3a68330956e4"

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
        args.tool = ["TRGT"]
    if not args.data_type:
        args.data_type = ["pacbio"]

    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]

    if args.data_type:
        df = df[df.sequencing_data_type.isin(args.data_type)]

    df = df[df.sample_id.isin({"HG002", "CHM1_CHM13"})]  # only use these samples for tool evaluations

    if args.custom_catalog_path and args.output_dir == DEFAULT_OUTPUT_DIR:
        parser.error("--custom-catalog-path is set without also setting --output-dir")

    bp.precache_file_paths(os.path.join(args.output_dir, "**/*.*"))


    download_to_dir = "../results"

    for row_i, (_, row) in enumerate(df.iterrows()):
        if args.filename_keyword:
            if not any(keyword in row.read_data_path for keyword in args.filename_keyword):
                continue

        coverage = int(round(float(row.depth_of_coverage)))
        for tool in args.tool:
            if tool in SHORT_READ_TOOLS and row.sequencing_data_type not in SHORT_READ_DATA_TYPES:
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                      f"doesn't support {row.sequencing_data_type} data")
                continue
            if tool in LONG_READ_TOOLS and row.sequencing_data_type not in LONG_READ_DATA_TYPES:
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                      f"doesn't support {row.sequencing_data_type} data")
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
            else:
                catalog_path_suffix = tool

                repeat_catalog_paths = os.path.join(args.filter_vcf_dir, row.sample_id,
                    f"{row.sample_id}.{catalog_path_suffix}*")

            print(f"Listing catalogs {repeat_catalog_paths}")
            repeat_catalog_paths = [x.path for x in hfs.ls(repeat_catalog_paths)]
            output_dir = os.path.join(args.output_dir, row.sample_id, row.sequencing_data_type, tool, f"{coverage}x_coverage")
            if tool in ("EHv5", "EHv5-bw2-optimized", "IlluminaEHv5"):
                # Three ExpansionHunter v5 variants, all genotyped with create_expansion_hunter_steps:
                #   EHv5               - bw2 fork, --analysis-mode low-mem-streaming
                #   EHv5-bw2-optimized - bw2 fork, --analysis-mode optimized-streaming --improved-genotyping
                #   IlluminaEHv5       - original Illumina build, --analysis-mode streaming
                use_illumina_expansion_hunter = (tool == "IlluminaEHv5")
                improved_genotyping = False
                if tool == "EHv5":
                    analysis_mode = "low-mem-streaming"
                elif tool == "EHv5-bw2-optimized":
                    analysis_mode = "optimized-streaming"
                    improved_genotyping = True
                else:  # IlluminaEHv5
                    analysis_mode = "streaming"

                if args.custom_catalog_path:
                    variant_catalog_file_paths = repeat_catalog_paths
                elif tool == "IlluminaEHv5":
                    # the unoptimized official build runs over the sharded catalogs (one parallel step per shard);
                    # fall back to the unsharded catalog if the converter only wrote a single shard
                    variant_catalog_file_paths = [p for p in repeat_catalog_paths if "001_of_001" not in p] \
                        or repeat_catalog_paths
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
                    improved_genotyping=improved_genotyping,
                    loci_to_exclude=None,
                    min_locus_coverage=None,
                    use_illumina_expansion_hunter=use_illumina_expansion_hunter)
            elif tool == "GangSTR":
                if row.sequencing_data_type == "ultima":
                    # for some reason GangSTR never completes on ultima data
                    print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                          f"doesn't support {row.sequencing_data_type} data")
                    continue
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
            elif tool == "TRGT":
                if row.sequencing_data_type != "pacbio":
                    print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                          f"doesn't support {row.sequencing_data_type} data")
                    continue

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
                    output_prefix= f"{row.sample_id}.{tool}")
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
                    male_or_female=row.male_or_female,
                    expansion_hunter_catalog_paths=[p for p in repeat_catalog_paths if "001_of_001" in p],
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
                coverage=coverage,
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
                coverage=coverage,
                sample_id=row.sample_id,
                output_dir=output_dir)
    bp.run()


def add_tool_comparison_columns_step(bp, tool_results_step, *, tool, coverage, sample_id, output_dir, truth_set_genotypes_path, tool2="Truth", download_to_dir=None):
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
df.loc[:, "Coverage"] = "{coverage}x"
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


def create_plot_tool_accuracy_steps(bp, add_columns_step, *, tool, coverage, sample_id, output_dir):
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
            f"--coverage {coverage}x "
            "--q-threshold 0 "
            f"--min-motif-size {min_motif_size} "
            f"--max-motif-size {max_motif_size} "
            "--genotype all "
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
    #                                f"--coverage {coverage}x "
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
