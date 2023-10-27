"""
This pipeline runs the TR genotyping tools on the catalog(s) of interest.
"""

import os


#%%

n = 6
expansion_hunter_variant_catalogs = " ".join([
    f"--variant-catalog gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.EHv5.00{i}_of_00{n}.json "
    for i in range(1, n+1)
])

gangstr_repeat_specs = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.GangSTR.001_of_001.bed"
hipstr_regions_bed = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.HipSTR.001_of_001.bed"
trgt_catalog_bed = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.TRGT_repeat_catalog.bed"

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