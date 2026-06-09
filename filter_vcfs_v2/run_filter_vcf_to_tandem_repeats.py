"""Hail Batch pipeline for running str_analysis.filter_vcf_to_tandem_repeats on DipCall ouptut VCFs generated from
aligning T2T assemblies to hg38. This is the next iteration of the str_analysis.filter_vcf_to_STRs method.
"""

#import hailtop.fs as hfs
import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/str-analysis@sha256:477380c41637a26f5e19d3c4ff6bae60de0cf1db620ebba8549c44df48d586b8"
#DOCKER_IMAGE = "us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis@sha256:16191eb046706d19f2cc031f06e12c4da65e3e5f2e6d2a606b1aa8331bc2acae"

def parse_args(bp):
    parser = bp.get_config_arg_parser()
    parser.add_argument("--exclude-homopolymers", action="store_true")
    parser.add_argument("--skip-combine-steps", action="store_true")
    parser.add_argument("--use-nonpreemptibles", action="store_true")
    parser.add_argument("--allow-multiple-trf-results-per-locus", action="store_true")
    parser.add_argument("--genotype-catalog", help="If specified, genotype each sample against this catalog BED instead "
                        "of the combined catalog produced by the merge step. Lets you --skip-filter-step --skip-combine-step "
                        "and genotype a subset of samples (via -s) against an existing catalog.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar in the genotype step.")
    parser.add_argument("-n", type=int, help="Number of samples to process")
    parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample. Can be specified more than once.")
    parser.add_argument("--metadata-tsv", default="../dipcall_pipeline/all_assemblies.tsv")
    parser.add_argument("--input-dir", default="gs://str-truth-set-v2/dipcall_pipeline")
    parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf_v2")
    parser.add_argument("--cpu", type=float, default=4)
    parser.add_argument("--memory", default="standard", choices=["lowmem", "standard", "highmem"])
    args = bp.parse_known_args()

    return args


def create_filter_step(bp, row, input_dir, output_dir,
                       allow_multiple_trf_results_per_locus=False,
                       exclude_homopolymers=False,
                       use_preemptibles=True,
                       cpu=4,
                       memory="lowmem"):

    filter_step = bp.new_step(
        f"filter_vcf_to_tandem_repeats (cpu={cpu}): {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="filter-step",
        preemptible=use_preemptibles,
        cpu=cpu,
        storage="10G",
        memory=memory,
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)

    hg38_fasta_input, _ = filter_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai")

    dipcall_input_dir = input_dir
    if row.get("subdirectory") and not pd.isna(row.get("subdirectory")):
        dipcall_input_dir = os.path.join(input_dir, row.subdirectory)
    dipcall_input_dir = os.path.join(dipcall_input_dir, row.sample_id)

    dipcall_vcf_input, dipcall_high_confidence_regions_bed_input = filter_step.inputs(
        os.path.join(dipcall_input_dir, f"{row.sample_id}.dip.vcf.gz"),
        os.path.join(dipcall_input_dir, f"{row.sample_id}.dip.bed.gz"))

    filter_step.command("set -exuo pipefail")

    filter_step.command(f"[ -s {dipcall_high_confidence_regions_bed_input} ] || exit 1")  # check that the bed file isn't emtpy

    filter_step.command(f"bedtools intersect -header -f 1 -wa -u \
            -a {dipcall_vcf_input}  \
            -b {dipcall_high_confidence_regions_bed_input} \
            | bgzip > {row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.command(f"tabix -f {row.sample_id}.high_confidence_regions.vcf.gz")

    filter_step.command(f"python3 -m pip uninstall -y str-analysis")
    filter_step.command(f"python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis")
    #filter_step.command(f"python3 -u -m str_analysis.filter_vcf_to_tandem_repeats catalog -h || true")

    min_repeat_unit_length = 2 if exclude_homopolymers else 1
    allow_multiple_arg = "--allow-multiple-trf-results-per-locus" if allow_multiple_trf_results_per_locus else ""
    filter_step.command(f"python3 -u -m str_analysis.filter_vcf_to_tandem_repeats catalog \
            -R {hg38_fasta_input} \
            --min-repeat-unit-length {min_repeat_unit_length} \
            --min-repeats 3 \
            --min-tandem-repeat-length 9 \
            --min-indel-size-to-run-trf 7 \
            --trf-min-repeats-in-reference 2 \
            --trf-min-purity 0.2 \
            --trf-executable-path /usr/bin/trf \
            --trf-threads {int(2*cpu)} \
            --write-detailed-bed \
            --write-tsv \
            --write-vcf \
            --verbose  {allow_multiple_arg} \
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


def create_combine_step(bp, filter_steps, data_dir, cpu=2, memory="highmem"):

    output_prefix = f"combined.{len(filter_steps)}_catalogs"
    combine_step = bp.new_step(
        f"combine (cpu={cpu}): {len(filter_steps):,d} catalogs",
        image=DOCKER_IMAGE,
        arg_suffix="combine-step",
        localize_by=Localize.COPY,
        cpu=cpu,
        memory=memory,
        storage="20G",
        output_dir=data_dir)

    catalog_bed_files = []
    for filter_step in filter_steps:
        combine_step.depends_on(filter_step)

        filter_step_outputs = filter_step.get_outputs()
        local_bed_path = combine_step.input(filter_step_outputs[4].output_path)
        catalog_bed_files.append(local_bed_path)

    combine_step.command("set -exuo pipefail")

    combine_step.command(f"python3 -m pip uninstall -y str-analysis")
    combine_step.command(f"python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis")
    hg38_fasta_input, _ = combine_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai")

    combine_step.command(f"python3 -u -m str_analysis.filter_vcf_to_tandem_repeats merge \
            -R {hg38_fasta_input} \
            --write-detailed-bed \
            --verbose \
            --output-prefix {output_prefix} \
            {' '.join([str(path) for path in catalog_bed_files])} |& tee {output_prefix}.log")

    combine_step.command("ls -lhtr")

    combine_step.output(f"{output_prefix}.tandem_repeats.bed.gz")
    combine_step.output(f"{output_prefix}.tandem_repeats.bed.gz.tbi")
    combine_step.output(f"{output_prefix}.tandem_repeats.detailed.bed.gz", download_to_dir="results")
    combine_step.output(f"{output_prefix}.tandem_repeats.detailed.bed.gz.tbi", download_to_dir="results")
    combine_step.output(f"{output_prefix}.log")

    return combine_step


def create_genotype_step(bp, row, combined_catalog_bed_path, filter_step, combine_step, output_dir,
                         cpu=4, memory="standard", use_preemptibles=True, show_progress_bar=False):

    genotype_step = bp.new_step(
        f"genotype (cpu={cpu}): {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="genotype-step",
        preemptible=use_preemptibles,
        cpu=cpu,
        storage="10G",
        memory=memory,
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)

    genotype_step.depends_on(filter_step)
    if combine_step is not None:
        genotype_step.depends_on(combine_step)

    hg38_fasta_input, _ = genotype_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai")

    catalog_bed_input = genotype_step.input(combined_catalog_bed_path)
    high_confidence_regions_vcf_input, _ = genotype_step.inputs(
        os.path.join(output_dir, f"{row.sample_id}.high_confidence_regions.vcf.gz"),
        os.path.join(output_dir, f"{row.sample_id}.high_confidence_regions.vcf.gz.tbi"))

    genotype_step.command("set -exuo pipefail")

    genotype_step.command(f"python3 -m pip uninstall -y str-analysis")
    genotype_step.command(f"python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis")

    genotype_step.command(f"python3 -u -m str_analysis.filter_vcf_to_tandem_repeats genotype \
            -R {hg38_fasta_input} \
            --catalog-bed {catalog_bed_input} \
            --write-json \
            --add-motif-composition trf \
            --trf-executable-path /usr/bin/trf \
            --trf-threads {int(2*cpu)} \
            {'--show-progress-bar' if show_progress_bar else ''} \
            --output-prefix {row.sample_id} \
            {high_confidence_regions_vcf_input} |& tee {row.sample_id}.genotype.log")

    genotype_step.command("ls -lhtr")

    genotype_step.output(f"{row.sample_id}.tandem_repeat_genotypes.tsv.gz")
    genotype_step.output(f"{row.sample_id}.tandem_repeat_genotypes.json.gz")
    genotype_step.output(f"{row.sample_id}.genotype.log")

    return genotype_step


def main():
    bp = pipeline("filter_vcf_to_tandem_repeats", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    args = parse_args(bp)
    bp.precache_file_paths(f"{args.output_dir}/**/*.*")

    df = pd.read_table(args.metadata_tsv)
    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]

    if args.n:
        df = df.iloc[:args.n]

    filter_steps = []
    for row_i, (_, row) in enumerate(df.iterrows()):
        output_dir = os.path.join(args.output_dir, row.sample_id)
        filter_step = create_filter_step(bp, row, args.input_dir, output_dir,
                                         allow_multiple_trf_results_per_locus=args.allow_multiple_trf_results_per_locus,
                                         exclude_homopolymers=args.exclude_homopolymers,
                                         use_preemptibles=not args.use_nonpreemptibles,
                                         cpu=args.cpu,
                                         memory=args.memory)
        
        filter_steps.append(filter_step)


    combine_step = create_combine_step(bp, filter_steps, args.output_dir)

    genotype_catalog_bed_path = args.genotype_catalog or combine_step.get_outputs()[0].output_path
    for (_, row), filter_step in zip(df.iterrows(), filter_steps):
        create_genotype_step(bp, row, genotype_catalog_bed_path, filter_step,
                             None if args.genotype_catalog else combine_step,
                             output_dir=os.path.join(args.output_dir, row.sample_id),
                             cpu=args.cpu,
                             memory=args.memory,
                             use_preemptibles=not args.use_nonpreemptibles,
                             show_progress_bar=args.show_progress_bar)

    bp.run()


if __name__ == "__main__":
    main()

#%%
