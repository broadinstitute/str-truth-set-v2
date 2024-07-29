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
    parser.add_argument("--output-dir", default="gs://str-truth-set-v2/tool_results")
    args = bp.parse_known_args()

    if not args.tool:
        args.tool = ["TRGT"]
    unknown_data_types = set(df.sequencing_data_type) - (SHORT_READ_DATA_TYPES | LONG_READ_DATA_TYPES)
    if unknown_data_types:
        raise ValueError(f"Unknown value(s) in the sequencing_data_type column of {sample_table_path}: {unknown_data_types}")

    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]
    if args.data_type:
        df = df[df.sequencing_data_type.isin(args.data_type)]

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

            repeat_catalog_paths = f"gs://str-truth-set-v2/filter_vcf/{output_dir_suffix}/{row.sample_id}/{row.sample_id}.STRs{excluding_homopolymers_string}.positive_loci.{tool}*"
            print(f"Listing catalogs {repeat_catalog_paths}")
            repeat_catalog_paths = [x.path for x in hfs.ls(repeat_catalog_paths)]
            output_dir = os.path.join(args.output_dir, output_dir_suffix, row.sample_id, row.sequencing_data_type, tool, f"{coverage}x_coverage")
            if tool == "EHv5":
                create_expansion_hunter_steps(
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
                create_gangstr_steps(
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
                create_hipstr_steps(
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


    bp.run()


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