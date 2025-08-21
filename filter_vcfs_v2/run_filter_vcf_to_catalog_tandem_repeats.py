"""Hail Batch pipeline for running str_analysis.filter_vcf_to_catalog_tandem_repeats on DipCall ouptut VCFs generated from aligning HPRC assemblies to hg38"""

#import hailtop.fs as hfs
import os
import pandas as pd

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/str-analysis@sha256:46fd53c06059b89a874c0eb438c2f761f1e96c98b9fa7bf9d83eb88e2756c3f8"

def create_filter_step(bp, row, input_dir, output_dir, cpu=4, memory="lowmem"):

    filter_step = bp.new_step(
        f"filter_vcf_to_catalog_tandem_repeats (cpu={cpu}): {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="filter-step",
        cpu=cpu,
        memory=memory,
        output_dir=output_dir)

    hg38_fasta_input, _ = filter_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    dipcall_vcf_input, dipcall_high_confidence_regions_bed_input = filter_step.inputs(
        os.path.join(input_dir, f"{row.sample_id}.dip.vcf.gz"),
        os.path.join(input_dir, f"{row.sample_id}.dip.bed.gz"),
    )

    filter_step.command("set -exuo pipefail")

    filter_step.command(f"[ -s {dipcall_high_confidence_regions_bed_input} ] || exit 1")  # check that the bed file isn't emtpy

    filter_step.command(f"bedtools intersect -header -f 1 -wa -u \
            -a {dipcall_vcf_input}  \
            -b {dipcall_high_confidence_regions_bed_input} \
            | bgzip > {row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.command(f"tabix -f {row.sample_id}.high_confidence_regions.vcf.gz")

    filter_step.command(f"python3 -u -m str_analysis.filter_vcf_to_catalog_tandem_repeats \
            -R {hg38_fasta_input} \
            --min-repeat-unit-length 1 \
            --min-repeats 3 \
            --min-tandem-repeat-length 9 \
            --trf-executable-path /usr/bin/trf \
            --trf-threads {int(cpu*1.5)} \
            --write-tsv \
            --write-vcf \
            --output-prefix {row.sample_id} \
            {row.sample_id}.high_confidence_regions.vcf.gz |& tee {row.sample_id}.filter_vcf.log")

    filter_step.command(f"tabix -f {row.sample_id}.tandem_repeats.vcf.gz")
    filter_step.command("ls -lhtr")

    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.bed.gz")
    filter_step.output(f"{row.sample_id}.tandem_repeats.bed.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.vcf.gz")
    filter_step.output(f"{row.sample_id}.tandem_repeats.vcf.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.tsv.gz")
    filter_step.output(f"{row.sample_id}.filter_vcf.log")

    return filter_step




def main():
    bp = pipeline("filter_vcf_to_catalog_tandem_repeats", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("-s", "--sample-id", action="append",
                        help="Process only this sample. Can be specified more than once.")
    parser.add_argument("--metadata-tsv", default="../dipcall_pipeline/hprc_assemblies_release2.tsv")
    parser.add_argument("--input-dir", default="gs://str-truth-set-v2/dipcall_pipeline/HPRC_release2")
    parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf_v2/HPRC_release2")
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--memory", default="lowmem", choices=["lowmem", "standard", "highmem"])
    args = bp.parse_known_args()

    bp.precache_file_paths(f"{args.output_dir}/**/*.*")

    df = pd.read_table(args.metadata_tsv)
    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]
    
    for row_i, (_, row) in enumerate(df.iterrows()):
        input_dir = os.path.join(args.input_dir, row.sample_id)

        output_dir = os.path.join(args.output_dir, row.sample_id)
        filter_step = create_filter_step(bp, row, input_dir, output_dir, cpu=args.cpu, memory=args.memory)


    bp.run()


if __name__ == "__main__":
    main()


#%%
