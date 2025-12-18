"""Run the new version of the filtering method:  str_analysis.filter_vcf_to_tandem_repeats"""

import collections
import hailtop.fs as hfs
import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/str-analysis@sha256:4de1c779f80ece7f538db344af74ea590afbf010ef4526caac4a9a8bbb1118d0"


def create_filter_step(bp, row, input_dir, output_dir, exclude_homopolymers=False, use_preemptibles=True, cpu=4):
    
    filter_step = bp.new_step(
        f"filter_vcf_to_tandem_repeats: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="filter-step",
        preemptible=use_preemptibles,
        cpu=cpu,
        memory="highmem",
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

    min_repeat_unit_length = 2 if exclude_homopolymers else 1
    filter_step.command(f"python3 -u -m str_analysis.filter_vcf_to_tandem_repeats catalog \
            -R {hg38_fasta_input} \
            --trf-executable-path trf \
            --trf-threads {cpu} \
            --min-indel-size-to-run-trf 7 \
            --min-tandem-repeat-length 9 \
            --min-repeats 3 \
            --min-repeat-unit-length {min_repeat_unit_length} \
            --output-prefix {row.sample_id} \
            --verbose \
            {row.sample_id}.high_confidence_regions.vcf.gz |& tee {row.sample_id}.filter_vcf.log")

    filter_step.command("ls -lhtr")

    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz.tbi")
    filter_step.output(f"{row.sample_id}.tandem_repeats.bed.gz")
    filter_step.output(f"{row.sample_id}.tandem_repeats.bed.gz.tbi")
    filter_step.output(f"{row.sample_id}.filter_vcf.log")

    return filter_step


def main():
    bp = pipeline("filter_vcf_to_tandem_repeats", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--exclude-homopolymers", action="store_true")
    parser.add_argument("--skip-combine-steps", action="store_true")
    parser.add_argument("--use-nonpreemptibles", action="store_true")
    parser.add_argument("-n", type=int, help="Number of samples to process")
    parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample. Can be specified more than once.")
    parser.add_argument("--metadata-tsv", default="../dipcall_pipeline/all_assemblies.tsv")
    parser.add_argument("--input-dir", default="gs://str-truth-set-v2/dipcall_pipeline")
    parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf_to_tandem_repeats")
    args = bp.parse_known_args()

    bp.precache_file_paths(f"{args.output_dir}/**/*.*")

    df = pd.read_table(args.metadata_tsv)
    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]

    if args.n:
        df = df.iloc[:args.n]

    filter_steps = []
    for row_i, (_, row) in enumerate(df.iterrows()):
        input_dir = os.path.join(args.input_dir, row.sample_id)
        output_dir = os.path.join(args.output_dir, row.sample_id)

        filter_step = create_filter_step(bp, row, input_dir, output_dir,
                                         exclude_homopolymers=args.exclude_homopolymers,
                                         use_preemptibles=not args.use_nonpreemptibles)

        filter_steps.append(filter_step)

    bp.run()


if __name__ == "__main__":
    main()


#%%
