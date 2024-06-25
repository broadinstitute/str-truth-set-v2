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

SHORT_READ_TOOLS = [
    "eh",
    "gangstr",
    "hipstr",
]

PACBIO_TOOLS = [
    "trgt",
    "longtr",
]

def main():
    bp = pipeline("filter_vcf_to_STRs", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--only-pure-repeats", action="store_true")
    parser.add_argument("--exclude-homopolymers", action="store_true")
    parser.add_argument("--skip-combine-steps", action="store_true")
    parser.add_argument("-s", "--sample-id", action="append",
                        help="Process only this sample. Can be specified more than once.")
    parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf")
    args = bp.parse_known_args()

    bp.precache_file_paths("gs://str-truth-set-v2/filter_vcf/**/*.*")

    df = pd.read_table("../run_tools/broad_short_read_cram_paths_for_HPRC_samples.txt")
    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]

    suffix = ".STRs"
    if args.only_pure_repeats:
        suffix += ".only_pure"
    if args.exclude_homopolymers:
        suffix += ".excluding_homopolymers"

    if args.only_pure_repeats and args.exclude_homopolymers:
        output_dir_suffix = "pure_repeats_excluding_homopolymers"
    elif args.only_pure_repeats:
        output_dir_suffix = "pure_repeats_including_homopolymers"
    elif args.exclude_homopolymers:
        output_dir_suffix = "all_repeats_excluding_homopolymers"
    else:
        output_dir_suffix = "all_repeats_including_homopolymers"

    filter_steps = []
    figures_to_download = collections.defaultdict(list)
    for row_i, (_, row) in enumerate(df.iterrows()):
        output_dir = os.path.join(args.output_dir, output_dir_suffix, row.sample_id)

        filter_step = create_filter_step(bp, row, suffix, output_dir,
                                         exclude_homopolymers=args.exclude_homopolymers,
                                         only_pure_repeats=args.only_pure_repeats)
        filter_steps.append(filter_step)


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

run(f"""
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/downsample_bam_pipeline.py \
    --verbose \
    --target-coverage 30 \
    --output-dir gs://bw2-delete-after-30-days/ \
    --input-bam {input_cram}
""")


run(f"""
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py \
    --input-bam {input_cram} \
    --input-bai {input_cram}.crai \
    --output-dir {output_dir}/eh \
    --use-streaming-mode \
    {expansion_hunter_variant_catalogs}
""")

run(f"""
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/gangstr_pipeline.py \
    --input-bam {input_cram} \
    --input-bai {input_cram}.crai \
    --output-dir {output_dir}/gangstr \
    --repeat-specs {gangstr_repeat_specs}
""")


run(f"""
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/hipstr_pipeline.py \
    --input-bam {input_cram} \
    --input-bai {input_cram}.crai \
    --output-dir {output_dir}/hipstr \
    --regions-bed {hipstr_regions_bed}
""")

run(f"""
python3 ../str-truth-set/tool_comparison/hail_batch_pipelines/trgt_pipeline.py \
    --input-bam {input_long_read_bam} \
    --input-bai {input_long_read_bam}.bai \
    --output-dir {output_dir}/trgt \
    --trgt-catalog-bed {trgt_catalog_bed}
""")

"""
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
