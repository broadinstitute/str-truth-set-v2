"""Hail Batch pipeline for running str_analysis.filter_vcf_to_catalog_tandem_repeats on DipCall ouptut VCFs generated from aligning HPRC assemblies to hg38"""

#import hailtop.fs as hfs
import os
import pandas as pd

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/str-analysis@sha256:17b9c6b289eb6042c4e9d0053f62e94d539d45cdd63b5f9b84b181453d1aea3b"

def parse_args(bp):
    parser = bp.get_config_arg_parser()
    parser.add_argument("-s", "--sample-id", action="append",
                        help="Process only this sample. Can be specified more than once.")
    
    parser.add_argument("--metadata-tsv", default="../dipcall_pipeline/all_assemblies.tsv")
    parser.add_argument("--input-dir", default="gs://str-truth-set-v2/dipcall_pipeline")
    parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf_v2")

    #parser.add_argument("--metadata-tsv", default="../dipcall_pipeline/hprc_assemblies_release2.tsv")
    #parser.add_argument("--input-dir", default="gs://str-truth-set-v2/dipcall_pipeline/HPRC_release2")
    #parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf_v2/HPRC_release2")

    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--memory", default="standard", choices=["lowmem", "standard", "highmem"])
    args = bp.parse_known_args()

    return args


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

    filter_step.command(f"python3 -u -m str_analysis.filter_vcf_to_tandem_repeats catalog \
            -R {hg38_fasta_input} \
            --min-repeat-unit-length 1 \
            --min-repeats 3 \
            --min-tandem-repeat-length 9 \
            --trf-executable-path /usr/bin/trf \
            --trf-threads {2*cpu} \
            --write-detailed-bed \
            --write-tsv \
            --write-vcf \
            --verbose \
            --output-prefix {row.sample_id} \
            {row.sample_id}.high_confidence_regions.vcf.gz |& tee {row.sample_id}.filter_vcf.log")

    filter_step.command("ls -lhtr")

    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.bed.gz")
    filter_step.output(f"{row.sample_id}.tandem_repeats.bed.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.detailed.bed.gz")
    filter_step.output(f"{row.sample_id}.tandem_repeats.detailed.bed.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.vcf.gz")
    filter_step.output(f"{row.sample_id}.tandem_repeats.vcf.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.tsv.gz")
    filter_step.output(f"{row.sample_id}.filter_vcf.log")

    return filter_step


def create_combine_step(bp, filter_steps, data_dir, cpu=1, memory="highmem"):

    output_prefix = f"combined.{len(filter_steps)}_catalogs"
    combine_step = bp.new_step(
        f"combine (cpu={cpu}): {len(filter_steps):,d} catalogs",
        image=DOCKER_IMAGE,
        arg_suffix="combine-step",
        cpu=cpu,
        memory=memory,
        output_dir=data_dir)

    catalog_bed_files = []
    for filter_step in filter_steps:
        combine_step.depends_on(filter_step)

        filter_step_outputs = filter_step.get_outputs()
        local_bed_path = combine_step.input(filter_step_outputs[4].output_path)
        catalog_bed_files.append(local_bed_path)

    combine_step.command("set -exuo pipefail")

    combine_step.command(f"python3 -u -m str_analysis.filter_vcf_to_tandem_repeats merge \
            --write-detailed-bed \
            --verbose \
            --output-prefix {output_prefix} \
            {' '.join([str(path) for path in catalog_bed_files])} |& tee {output_prefix}.log")

    combine_step.command("ls -lhtr")

    combine_step.output(f"{output_prefix}.tandem_repeats.detailed.bed.gz")
    combine_step.output(f"{output_prefix}.tandem_repeats.detailed.bed.gz.tbi")
    combine_step.output(f"{output_prefix}.log")

    return combine_step



def main():
    bp = pipeline("filter_vcf_to_catalog_tandem_repeats", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    args = parse_args(bp)
    bp.precache_file_paths(f"{args.output_dir}/**/*.*")

    df = pd.read_table(args.metadata_tsv)
    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]
    
    filter_steps = []
    for row_i, (_, row) in enumerate(df.iterrows()):
        input_dir = os.path.join(args.input_dir, row.sample_id)

        output_dir = os.path.join(args.output_dir, row.sample_id)
        filter_step = create_filter_step(bp, row, input_dir, output_dir, cpu=args.cpu, memory=args.memory)
        
        filter_steps.append(filter_step)


    combine_step = create_combine_step(bp, filter_steps, args.output_dir)

    bp.run()


if __name__ == "__main__":
    main()


#%%
