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
from step_pipeline import pipeline, Backend, Localize
import sys

sys.path.append("../str-truth-set/tool_comparison/hail_batch_pipelines")
from trgt_pipeline import create_trgt_step
from longtr_pipeline import create_longtr_steps
from straglr_pipeline import create_straglr_steps
from expansion_hunter_pipeline import create_expansion_hunter_steps
from gangstr_pipeline import create_gangstr_steps
from hipstr_pipeline import create_hipstr_steps


SHORT_READ_TOOLS = {
    "EHv5",
    "GangSTR",
    "HipSTR",
}

LONG_READ_TOOLS = {
    "TRGT",
    "LongTR",
    "straglr",
}

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

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

FILTER_VCFS_DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:ba23c9ddab4f1c9e599fe2b57d83a37e54aa5a5e38d4dfa27c277e4d6d83d225"

DEFAULT_OUTPUT_DIR = "gs://str-truth-set-v2/tool_results"

def main():
    sample_table_path = "HPRC_all_aligned_short_read_and_long_read_samples.tsv"
    df = pd.read_table(sample_table_path)

    bp = pipeline("run_genotyping_tools", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--only-pure-repeats", action="store_true")
    parser.add_argument("--exclude-homopolymers", action="store_true")
    parser.add_argument("--skip-combine-steps", action="store_true")
    parser.add_argument("-s", "--sample-id", action="append",
                        help="Process only this sample. Can be specified more than once.")
    parser.add_argument("-t", "--tool", action="append", choices=SHORT_READ_TOOLS|LONG_READ_TOOLS, help="The tool to run.")
    parser.add_argument("--data-type", action="append", choices=SHORT_READ_DATA_TYPES|LONG_READ_DATA_TYPES, help="Which data type(s) to process")
    parser.add_argument("-k", "--filename-keyword", help="If specified, only BAM paths that contain this keyword will be processed", action="append")
    parser.add_argument("--filter-vcf-dir", default="gs://str-truth-set-v2/filter_vcf", help="Base dir for filter_vcf pipeline output files")
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
        paser.error("--custom-catalog-path is set without also setting --output-dir")

    bp.precache_file_paths(os.path.join(args.output_dir, "**/*.*"))

    suffix = ".STRs"
    if args.only_pure_repeats:
        suffix += ".only_pure"
    if args.exclude_homopolymers:
        suffix += ".excluding_homopolymers"

    pure_repeats_or_all_repeats = "pure_repeats" if args.only_pure_repeats else "all_repeats"
    including_or_excluding_homopolymers = "including_homopolymers" if not args.exclude_homopolymers else "excluding_homopolymers"
    excluding_homopolymers_string = ".excluding_homopolymers" if args.exclude_homopolymers else ""
    output_dir_suffix = f"{pure_repeats_or_all_repeats}_{including_or_excluding_homopolymers}"

    filter_steps = []
    figures_to_download = collections.defaultdict(list)
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
            else:
                repeat_catalog_paths = os.path.join(args.filter_vcf_dir, output_dir_suffix, row.sample_id,
                    f"{row.sample_id}.STRs{excluding_homopolymers_string}.positive_loci.{tool}*")

            print(f"Listing catalogs {repeat_catalog_paths}")
            repeat_catalog_paths = [x.path for x in hfs.ls(repeat_catalog_paths)]
            output_dir = os.path.join(args.output_dir, output_dir_suffix, row.sample_id, row.sequencing_data_type, tool, f"{coverage}x_coverage")
            if tool == "EHv5":
                current_step = create_expansion_hunter_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    variant_catalog_file_paths=[p for p in repeat_catalog_paths if "001_of_001" not in p], # exclude the unsharded catalog
                    output_dir=output_dir,
                    output_prefix= f"{row.sample_id}.STRs.positive_loci.{tool}",
                    use_streaming_mode=False,
                    loci_to_exclude=None,
                    min_locus_coverage=None,
                    use_illumina_expansion_hunter=False)
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
                    output_prefix=f"{row.sample_id}.STRs.positive_loci.{tool}")
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
                    output_prefix=f"{row.sample_id}.STRs.positive_loci.{tool}")
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
                    output_prefix= f"{row.sample_id}.STRs.positive_loci.{tool}")
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
                    output_prefix=f"{row.sample_id}.STRs.positive_loci.{tool}")
            elif tool == "straglr":
                current_step = create_straglr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    straglr_catalog_bed_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.STRs.positive_loci.{tool}")
            else:
                raise ValueError(f"Unknown tool: {tool}")


            add_columns_step = add_tool_comparison_columns_step(
                bp,
                current_step,
                tool=tool,
                coverage=coverage,
                sample_id=row.sample_id,
                output_dir=output_dir,
                filter_vcf_dir=os.path.join(args.filter_vcf_dir, output_dir_suffix, row.sample_id),
                suffix=suffix,
                tool2="Truth")

            plot_tool_accuracy_step = create_plot_tool_accuracy_steps(
                bp,
                add_columns_step,
                tool=tool,
                coverage=coverage,
                sample_id=row.sample_id,
                output_dir=output_dir)
    bp.run()


def add_tool_comparison_columns_step(bp, tool_results_step, *, tool, coverage, sample_id, output_dir, filter_vcf_dir, suffix, tool2="Truth"):
    if tool == "EHv5":
        tool = "ExpansionHunter"

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
        cpu=1,
        memory="highmem",
        storage="20Gi",
        output_dir=output_dir)

    add_columns_step.depends_on(tool_results_step)

    add_columns_step.command("set -ex")

    local_tool_results_input = add_columns_step.input(tool_results_path)
    local_truth_set_input = add_columns_step.input(
        os.path.join(filter_vcf_dir, f"{sample_id}{suffix}.annotated.variants.for_comparison.tsv.gz"))

    add_columns_step.command(f"""python3 <<EOF
import pandas as pd
print("Adding columns to {local_tool_results_input}")
df = pd.read_table("{local_tool_results_input}", dtype=str)
df.loc[:, "Coverage"] = "{coverage}x"
print(f"Writing {{len(df):,d}} rows to {local_tool_results_input}")
df.to_csv("{local_tool_results_input}", sep="\\t", index=False, header=True)
EOF
""")

    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/add_tool_results_columns.py "
               f"--tool {tool} "
               f"{local_tool_results_input} " 
               f"{local_truth_set_input} ")

    add_columns_step.command("ls -lhrt")

    local_tsv_file_path = str(local_truth_set_input).replace(".tsv.gz", "") + f".with_{tool}_results.tsv.gz"
    output_filename = local_truth_set_input.filename.replace(".tsv.gz", "") + f".with_{tool}_vs_{tool2}_columns.tsv.gz"
    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/add_concordance_columns.py "
               f"--tool {tool} "
               f"--compare-to {tool2} "
               f"--output-tsv {output_filename} "
               f"{local_tsv_file_path}")

    add_columns_step.command("ls -lhrt")

    add_columns_step.output(output_filename)
    add_columns_step.output(output_filename.replace(".tsv", ".alleles.tsv"))

    return add_columns_step


def create_plot_tool_accuracy_steps(bp, add_columns_step, *, tool, coverage, sample_id, output_dir):
    if tool == "EHv5":
        tool = "ExpansionHunter"


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

    for min_motif_size, max_motif_size in [(2, 6), (7, 1000)]:
        plot_tool_accuracy_step.command(
            f"python3 -u /str-truth-set/figures_and_tables/plot_tool_accuracy_by_allele_size.py "
            "--verbose "
            f"--tool {tool} "
            f"--coverage {coverage}x "
            "--q-threshold 0 "
            f"--min-motif-size {min_motif_size} " +
            (f"--max-motif-size {max_motif_size} " if max_motif_size is not None else "") +
            "--genotype all "
            "--show-no-call-loci "
            "--image-type svg "
            "--show-title "
            f"{local_alleles_tsv} ")

        plot_tool_accuracy_step.command("ls -lhrt")
        for output_filename in [
            f"tool_accuracy_by_true_allele_size.{min_motif_size}to{max_motif_size}bp_motifs.all_genotypes.{coverage}x.{tool}.svg",
            f"tool_accuracy_by_true_allele_size.{min_motif_size}to{max_motif_size}bp_motifs.all_genotypes.pure_repeats.{coverage}x.{tool}.svg",
            f"tool_accuracy_by_true_allele_size.{min_motif_size}to{max_motif_size}bp_motifs.all_genotypes.with_interruptions.{coverage}x.{tool}.svg"
        ]:
            plot_tool_accuracy_step.command(f"gzip {output_filename}")
            plot_tool_accuracy_step.command(f"mv {output_filename}.gz {output_filename}")
            plot_tool_accuracy_step.output(output_filename)


    # create step to run gcloud storage objects update --content-type 'image/svg+xml' --content-encoding 'gzip' gs://str-truth-set-v2/tool_results/all_repeats_excluding_homopolymers/HG002/**/*.svg.gz
    # on each image
    image_headers_step = bp.new_step(
        name=f"Set image headers for {sample_id} {tool} accuracy plots",
        arg_suffix="image-headers-step",
        image=FILTER_VCFS_DOCKER_IMAGE,
        cpu=1,
        output_dir=output_dir)

    image_headers_step.depends_on(plot_tool_accuracy_step)

    image_headers_step.command("set -ex")
    for previous_step_output in plot_tool_accuracy_step.get_outputs():
        if previous_step_output.filename.endswith(".svg") or previous_step_output.filename.endswith(".svg.gz"):
            image_headers_step.command(
                f"gcloud storage objects update --content-type 'image/svg+xml' --content-encoding 'gzip' {previous_step_output.output_path}")

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


#%%

# python3 trgt_pipeline.py \
#   --input-bam gs://bw2-delete-after-60-days/long-reads/CHM1_CHM13.subreads.bam \
#   --input-bai gs://bw2-delete-after-60-days/long-reads/CHM1_CHM13.subreads.bam.bai \
#   --output-dir gs://bw2-delete-after-60-days/excluding_homopolymers/trgt \
#   --trgt-catalog-bed gs://str-truth-set-v2/filter_vcf/all_repeats_excluding_homopolymers/CHM1_CHM13/positive_loci.TRGT_repeat_catalog.bed \
#   --check-file-last-modified-times

# python3 trgt_pipeline.py \
#   --input-bam gs://str-truth-set-v2/raw_data/HG002/long_reads/HG002.bam \
#   --input-bai gs://str-truth-set-v2/raw_data/HG002/long_reads/HG002.bam.bai \
#   --output-dir gs://bw2-delete-after-60-days/excluding_homopolymers/trgt \
#   --trgt-catalog-bed gs://str-truth-set-v2/filter_vcf/all_repeats_excluding_homopolymers/HG002/positive_loci.TRGT_repeat_catalog.bed \
#   --check-file-last-modified-times

# python3 longtr_pipeline.py \
#   --regions-bed gs://str-truth-set-v2/filter_vcf/all_repeats_excluding_homopolymers/CHM1_CHM13/positive_loci.LongTR.001_of_001.bed \
#   --input-bam gs://bw2-delete-after-60-days/long-reads/CHM1_CHM13.subreads.bam \
#   --input-bai gs://bw2-delete-after-60-days/long-reads/CHM1_CHM13.subreads.bam.bai \
#   --output-dir gs://bw2-delete-after-60-days/excluding_homopolymers/longtr

# python3 longtr_pipeline.py  \
#        --regions-bed gs://str-truth-set-v2/filter_vcf/all_repeats_excluding_homopolymers/HG002/positive_loci.HipSTR.001_of_001.bed   \
#        --input-bam gs://str-truth-set-v2/raw_data/HG002/long_reads/HG002.bam   \
#        --input-bai gs://str-truth-set-v2/raw_data/HG002/long_reads/HG002.bam.bai   \
#        --output-dir gs://bw2-delete-after-60-days/excluding_homopolymers/longtr_original_bam

# python3 longtr_pipeline.py  \
#        --regions-bed gs://str-truth-set-v2/filter_vcf/all_repeats_excluding_homopolymers/HG002/positive_loci.HipSTR.001_of_001.bed   \
#        --input-bam gs://bw2-delete-after-60-days/long-reads/HG002.aligned.bam \
#        --input-bai gs://bw2-delete-after-60-days/long-reads/HG002.aligned.bam.bai   \
#        --output-dir gs://bw2-delete-after-60-days/excluding_homopolymers/longtr_egor_bam

# python3 trgt_pipeline.py \
#       --trgt-catalog-bed gs://str-truth-set-v2/filter_vcf/all_repeats_excluding_homopolymers/HG002/positive_loci.TRGT_repeat_catalog.bed \
#       --input-bam gs://bw2-delete-after-60-days/long-reads/HG002.aligned.bam \
#       --input-bai gs://bw2-delete-after-60-days/long-reads/HG002.aligned.bam.bai \
#       --output-dir gs://bw2-delete-after-60-days/excluding_homopolymers/trgt_egor_bam \

"""
n = 6
expansion_hunter_variant_catalogs = " ".join([
    f"--variant-catalog gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.EHv5.00{i}_of_00{n}.json "
    for i in range(1, n+1)
])

gangstr_repeat_specs = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.GangSTR.001_of_001.bed"
hipstr_regions_bed = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.HipSTR.001_of_001.bed"
trgt_catalog_bed = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.TRGT_repeat_catalog.bed"
longtr_catalog_bed = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.LongTR.001_of_001.bed"

input_cram = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
input_30x_cram = "gs://bw2-delete-after-30-days/CHM1_CHM13_WGS2.downsampled_to_30x.bam"
input_long_read_bam = "gs://bw2-delete-after-30-days/long-reads/CHM1_CHM13.subreads.bam"
output_dir = "gs://bw2-delete-after-30-days"


def run(cmd):
    print(cmd)
    #os.system(cmd)

run(f"
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/downsample_bam_pipeline.py \
    --verbose \
    --target-coverage 30 \
    --output-dir gs://bw2-delete-after-30-days/ \
    --input-bam {input_cram}
")


run(f"
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py \
    --input-bam {input_cram} \
    --input-bai {input_cram}.crai \
    --output-dir {output_dir}/eh \
    --use-streaming-mode \
    {expansion_hunter_variant_catalogs}
")

run(f"
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/gangstr_pipeline.py \
    --input-bam {input_cram} \
    --input-bai {input_cram}.crai \
    --output-dir {output_dir}/gangstr \
    --repeat-specs {gangstr_repeat_specs}
")


run(f"
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/hipstr_pipeline.py \
    --input-bam {input_cram} \
    --input-bai {input_cram}.crai \
    --output-dir {output_dir}/hipstr \
    --regions-bed {hipstr_regions_bed}
")

run(f"
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/trgt_pipeline.py \
    --input-bam {input_long_read_bam} \
    --input-bai {input_long_read_bam}.bai \
    --output-dir {output_dir}/trgt \
    --trgt-catalog-bed {trgt_catalog_bed}
")


truth_set=HG00438.STRs.excluding_homopolymers.annotated.variants.tsv.gz
python3 ../../str-truth-set/tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py \
    --output-dir . \
    $truth_set

truth_set_for_comparison=HG00438.STRs.excluding_homopolymers.annotated.variants.for_comparison.tsv.gz
tool_results=combined.positive_loci.EHv5.001_of_001.1_json_files.variants.tsv.gz 
python3 -u ../../str-truth-set/tool_comparison/scripts/add_tool_results_columns.py \
    --tool ExpansionHunter \
    $tool_results \
    $truth_set_for_comparison 

python3 -u ../../str-truth-set/tool_comparison/scripts/add_concordance_columns.py \
    --tool ExpansionHunter \
    --compare-to Truth \
    HG00438.STRs.excluding_homopolymers.annotated.variants.for_comparison.with_ExpansionHunter_results.tsv.gz
"""



#######
"""
python3 -u ../../str-truth-set/tool_comparison/scripts/add_tool_results_columns.py \
    --tool NewTruthSet \
    --filter-to-regions CHM1_CHM13.dip.bed.gz \
    CHM1_CHM13.STRs.excluding_homopolymers.variants.tsv.gz \
    STR_truth_set.v1.for_comparison.variants.tsv.gz

python3 -u ../../str-truth-set/tool_comparison/scripts/add_concordance_columns.py \
    --tool NewTruthSet \
    --compare-to Truth \
    STR_truth_set.v1.for_comparison.variants.with_NewTruthSet_results.tsv.gz

python3 ../../str-truth-set/figures_and_tables/plot_tool_accuracy_by_allele_size.py \
    STR_truth_set.v1.for_comparison.variants.with_NewTruthSet_results.with_concordance.alleles.tsv.gz \
    --tool NewTruthSet \
    --coverage 30x \
    --q-threshold 0 \
    --min-motif-size 2 \
    --max-motif-size 6 \
    --genotype all \
    --show-title
"""


"""
# Run tool pipelines using Hail Batch


force=""
#force="--force"

debug=""
#debug="echo"
if [ -z "$debug" ]; then
  set -ex
fi

function run_pipelines {
    input_bam=$1
    input_bai=$2
    output_dir=$3
    local_dir=$4
    run_illumina_expansion_hunter=$5
    min_locus_coverage_arg=$6

    # ExpansionHunter
    $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg} ${force} &
    $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg} ${force} &

    if [ "$run_illumina_expansion_hunter" == "yes" ]; then
	      $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --use-illumina-expansion-hunter --loci-to-exclude ./tool_comparison/hail_batch_pipelines/truth_set_loci_that_cause_illumina_expansion_hunter_error.txt \
            --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg}  ${force} &
	      $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --use-illumina-expansion-hunter --loci-to-exclude ./tool_comparison/hail_batch_pipelines/negative_loci_that_cause_illumina_expansion_hunter_error.txt \
	          --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg}  ${force} &
    fi

    # GangSTR
    $debug python3 ./tool_comparison/hail_batch_pipelines/gangstr_pipeline.py  --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/gangstr ${force} &
    $debug python3 ./tool_comparison/hail_batch_pipelines/gangstr_pipeline.py  --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/gangstr ${force} &

    # HipSTR
    $debug python3 ./tool_comparison/hail_batch_pipelines/hipstr_pipeline.py  --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/hipstr ${force} &
    $debug python3 ./tool_comparison/hail_batch_pipelines/hipstr_pipeline.py  --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/hipstr ${force} &

    # wait for all pipelines except ExpansionHunterDenovo because that one takes a long time and uses few jobs, so it
    # can just run in the background
    wait_for_pids="$(jobs -p)"

    # ExpansionHunterDenovo
    $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_denovo_pipeline.py --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter_denovo &

    wait < <(echo ${wait_for_pids})

    # download results
    $debug mkdir -p ${local_dir}/expansion_hunter_denovo
    $debug gsutil -m cp ${output_dir}/expansion_hunter_denovo/CHM*.locus.tsv ${local_dir}/expansion_hunter_denovo/

    $debug mkdir -p ${local_dir}/expansion_hunter/positive_loci/      ${local_dir}/expansion_hunter/negative_loci/
    $debug gsutil -m cp ${output_dir}/expansion_hunter/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/positive_loci/
    $debug gsutil -m cp ${output_dir}/expansion_hunter/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/negative_loci/

    $debug mkdir -p ${local_dir}/gangstr/positive_loci/      ${local_dir}/gangstr/negative_loci/
    $debug gsutil -m cp ${output_dir}/gangstr/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/positive_loci/
    $debug gsutil -m cp ${output_dir}/gangstr/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/negative_loci/

    $debug mkdir -p ${local_dir}/hipstr/positive_loci/      ${local_dir}/hipstr/negative_loci/
    $debug gsutil -m cp ${output_dir}/hipstr/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/hipstr/positive_loci/
    $debug gsutil -m cp ${output_dir}/hipstr/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/hipstr/negative_loci/
}

# Downsample coverage
for coverage in 30 20 10 5; do
    $debug python3 ./tool_comparison/hail_batch_pipelines/downsample_bam_pipeline.py --verbose --target-coverage ${coverage} \
	--output-dir gs://bw2-delete-after-30-days/ --input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram &
done


# Run pipelines on original bam
run_pipelines \
  "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram" \
  "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai" \
  "gs://str-truth-set/hg38/tool_results" \
  "./tool_comparison/results" \
  "yes" \
  ""


# Run pipelines on exome data
run_pipelines \
  "gs://broad-public-datasets/CHM1_CHM13_WES/CHMI_CHMI3_Nex1.cram" \
  "gs://broad-public-datasets/CHM1_CHM13_WES/CHMI_CHMI3_Nex1.cram.crai" \
  "gs://str-truth-set/hg38/tool_results_for_exome" \
  "./tool_comparison/results_for_exome" \
  "no" \
  ""

wait   # wait for downsampling to finish

# Process other coverage levels
for coverage in 30 20 10 5; do
    if [ "$coverage" == "10" ] || [ "$coverage" == "5" ]; then
	      min_locus_coverage="--min-locus-coverage 3"
    else
	      min_locus_coverage=""
    fi

    # Rerun pipelines on downsampled bam
    run_pipelines \
      "gs://bw2-delete-after-30-days/CHM1_CHM13_WGS2.downsampled_to_${coverage}x.bam" \
      "gs://bw2-delete-after-30-days/CHM1_CHM13_WGS2.downsampled_to_${coverage}x.bam.bai" \
      "gs://str-truth-set/hg38/tool_results_for_downsampled_${coverage}x_bam" \
      "./tool_comparison/results_for_downsampled_${coverage}x_bam" \
      "no" \
      "${min_locus_coverage}"
done

# Run REViewer
python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --positive-loci --run-reviewer \
	--input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram --input-bai gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai \
	--output-dir gs://str-truth-set/hg38/tool_results/expansion_hunter &
python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --negative-loci --run-reviewer \
	--input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram --input-bai gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai \
	--output-dir gs://str-truth-set/hg38/tool_results/expansion_hunter &

python3 ./tool_comparison/scripts/add_reviewer_image_url_to_bed.py -i ./STR_truth_set.v1.variants.bed.gz
python3 ./tool_comparison/scripts/add_reviewer_image_url_to_bed.py -i ./tool_comparison/variant_catalogs/negative_loci.bed.gz

gsutil -m cp ./STR_truth_set.v1.variants.with_reviewer_image_urls.bed.gz* gs://str-truth-set/hg38/
gsutil -m cp ./tool_comparison/variant_catalogs/negative_loci.with_reviewer_image_urls.bed.gz* gs://str-truth-set/hg38/tool_comparison/variant_catalogs/

wait  # wait for any remaining jobs to finish

echo Done with step C

"""



# ===============================
"""
set -x
set -e
set -u

# ExpansionHunter, GangSTR, HipSTR
for results_folder in \
  results_for_exome \
  results \
  results_for_downsampled_30x_bam \
  results_for_downsampled_20x_bam \
  results_for_downsampled_10x_bam \
  results_for_downsampled_5x_bam
do
  # compute STR_truth_set.v1.variants.for_comparison.tsv.gz
  python3 -u tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py \
      --output-dir ./tool_comparison/${results_folder}/  \
      STR_truth_set.v1.variants.tsv.gz \
      ./tool_comparison/variant_catalogs/negative_loci.tsv.gz

  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.variants.for_comparison.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz

  for i in 1 2 3;
  do
      if [ $i == 1 ]; then
        tool=ExpansionHunter
        subfolder=expansion_hunter
      elif [ $i == 2 ]; then
        tool=GangSTR
        subfolder=gangstr
      elif [ $i == 3 ]; then
        tool=HipSTR
        subfolder=hipstr
      else
        echo ERROR: unexpected i == $i
        exit 1;
      fi


      set +x
      echo ============================================
      echo Processing $results_folder $tool results
      set -x

      python3 -u tool_comparison/scripts/add_tool_results_columns.py \
	      --tool $tool \
	      ./tool_comparison/${results_folder}/${subfolder}/positive_loci/combined.positive_loci.*_json_files.variants.tsv.gz \
	      ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz

      python3 -u tool_comparison/scripts/add_tool_results_columns.py \
	      --tool $tool \
	      ./tool_comparison/${results_folder}/${subfolder}/negative_loci/combined.negative_loci.*_json_files.variants.tsv.gz \
	      ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

      mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_${tool}_results.tsv.gz \
         ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz
      mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_${tool}_results.tsv.gz \
         ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz
  done

  # add tool vs. truth set concordance columns
  python3 -u tool_comparison/scripts/add_concordance_columns.py \
    ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz

  python3 -u tool_comparison/scripts/add_concordance_columns.py \
    ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

  # rename files
  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_concordance.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.variants.tsv.gz
  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_concordance.alleles.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.alleles.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_concordance.tsv.gz \
     ./tool_comparison/${results_folder}/negative_loci.for_comparison.variants.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_concordance.alleles.tsv.gz \
     ./tool_comparison/${results_folder}/negative_loci.for_comparison.alleles.tsv.gz

  gunzip -f ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.variants.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/negative_loci.for_comparison.variants.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.alleles.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/negative_loci.for_comparison.alleles.tsv.gz

  # clean up intermediate files
  rm ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz
  rm ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz
done

# generate EHdn comparison tables for each depth-of-coverage
python3 -u ./tool_comparison/scripts/intersect_expansion_hunter_denovo_results_with_truth_set.py \
  --window-size 600 --truth-set-variants-tsv ./STR_truth_set.v1.variants.tsv.gz  \
  ./tool_comparison/results*/expansion_hunter_denovo/CHM1_CHM13_WGS2.*locus.tsv

# combine all
python3 -u ./tool_comparison/scripts/combine_all_results_tables.py

gsutil -m cp ./tool_comparison/combined.results.*.tsv.gz  gs://str-truth-set/hg38/

echo Done with step D
"""