"""Hail Batch pipeline for running dipcall on HPRC assemblies.

Relevant links:

HPRC assemblies
  https://projects.ensembl.org/hprc/

The design and construction of reference pangenome graphs with minigraph by Li et al. 2020
  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z

Increased mutation and gene conversion within human segmental duplications by Vollger et al.
  https://www.nature.com/articles/s41586-023-05895-y

"""
import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize


DOCKER_IMAGE = "weisburd/hprc-pipeline@sha256:656e50b22fb8a5c28cd8aa66678f467f280a0097e0d26121d0eadaacc97649bb"

bp = pipeline("filter_vcf_to_STRs", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

parser = bp.get_config_arg_parser()
parser.add_argument("--only-pure-repeats", action="store_true")
parser.add_argument("--exclude-homopolymers", action="store_true")
parser.add_argument("-s", "--sample-id", action="append",
                    help="Process only this sample. Can be specified more than once.")
parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf")
args = bp.parse_known_args()

df = pd.read_table("hprc_assemblies.tsv")
if args.sample_id:
    df = df[df.sample_id.isin(args.sample_id)]

suffix = ".STRs"
if args.only_pure_repeats:
    suffix += ".only_pure"
if args.exclude_homopolymers:
    suffix += ".excluding_homopolymers"

filter_steps = []
for i, (_, row) in enumerate(df.iterrows()):

    filter_step = bp.new_step(
        f"filter_vcf: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="step1",
        cpu=1,
        memory="highmem",
        #storage="50Gi",
        output_dir=os.path.join(args.output_dir, row.sample_id))

    hg38_fasta_input, _ = filter_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.COPY)

    dipcall_vcf_input, dipcall_high_confidence_regions_bed_input = filter_step.inputs(
        f"gs://str-truth-set-v2/hprc_dipcall/{row.sample_id}/{row.sample_id}.dip.vcf.gz",
        f"gs://str-truth-set-v2/hprc_dipcall/{row.sample_id}/{row.sample_id}.dip.bed",
    )

    filter_step.command("set -exuo pipefail")

    filter_step.command(f"[ -s {dipcall_high_confidence_regions_bed_input} ] || exit 1")  # check that the bed file isn't emtpy

    filter_step.command(f"bedtools intersect -header -f 1 -wa -u \
        -a {dipcall_vcf_input}  \
        -b {dipcall_high_confidence_regions_bed_input} \
        | bgzip > {row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.command(f"tabix -f {row.sample_id}.high_confidence_regions.vcf.gz")

    min_repeat_unit_length = 2 if args.exclude_homopolymers else 1
    allow_interruptions = "no" if args.only_pure_repeats else "only-if-pure-repeats-not-found"
    filter_step.command(f"python3 -u -m str_analysis.filter_vcf_to_STR_variants \
        -R {hg38_fasta_input} \
        --allow-interruptions {allow_interruptions} \
        --write-bed-file \
        --min-str-length 9 \
        --min-str-repeats 3 \
        --min-repeat-unit-length {min_repeat_unit_length} \
        --max-repeat-unit-length 50 \
        --output-prefix {row.sample_id}{suffix} \
        --verbose \
        {row.sample_id}.high_confidence_regions.vcf.gz |& tee {row.sample_id}{suffix}.filter_vcf.log")

    filter_step.command("ls -lhtr")

    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.output(f"{row.sample_id}.high_confidence_regions.vcf.gz.tbi")
    filter_step.output(f"{row.sample_id}{suffix}.vcf.gz")
    filter_step.output(f"{row.sample_id}{suffix}.vcf.gz.tbi")
    filter_step.output(f"{row.sample_id}{suffix}.variants.tsv.gz")
    filter_step.output(f"{row.sample_id}{suffix}.alleles.tsv.gz")
    filter_step.output(f"{row.sample_id}{suffix}.variants.bed.gz")
    filter_step.output(f"{row.sample_id}{suffix}.variants.bed.gz.tbi")
    filter_step.output(f"{row.sample_id}{suffix}.filter_vcf.log")
    filter_steps.append(filter_step)

    # annotate the output tables
    annotate_step = bp.new_step(
        f"annotate: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="step2",
        cpu=1,
        output_dir=os.path.join(args.output_dir, row.sample_id),
        depends_on=[filter_step],
    )
    annotate_step.command("set -exuo pipefail")
    annotate_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/annotate_variant_table.py ")

    # generate variant catalogs step
    variant_catalogs_step = bp.new_step(
        f"variant catalogs: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="step2",
        cpu=1,
        output_dir=os.path.join(args.output_dir, row.sample_id),
        depends_on=[filter_step],
    )

    high_confidence_regions_bed_input = variant_catalogs_step.input(
        os.path.join("gs://str-truth-set-v2/hprc_dipcall", row.sample_id, f"{row.sample_id}.dip.bed.gz"))
    variants_tsv_input = variant_catalogs_step.input(
        os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.variants.tsv.gz"))
    variants_bed_input = variant_catalogs_step.input(
        os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.variants.bed.gz"))
    all_hg38_repeats_bed_input, _ = variant_catalogs_step.inputs(
        "gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz",
        "gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz.tbi")

    variant_catalogs_step.command(
        "python3 -u /str-truth-set/tool_comparison/scripts/convert_truth_set_to_variant_catalogs.py "
        "--output-negative-loci "
        "--expansion-hunter-loci-per-run 1000000 "
        "--skip-gangstr-catalog "
        "--skip-hipstr-catalog "
        #"--gangstr-loci-per-run 1000000 "
        f"--high-confidence-regions-bed {high_confidence_regions_bed_input} "
        f"--all-hg38-repeats-bed {all_hg38_repeats_bed_input} "
        f"--truth-set-bed {variants_bed_input} "
        f"{variants_tsv_input}")

    variant_catalogs_step.command("ls -l")

    # generate plots
    plot_step = bp.new_step(
        f"plots: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="step3",
        output_dir=os.path.join(args.output_dir, row.sample_id),
        depends_on=[filter_step],
    )

    plot_step.command("set -exuo pipefail")

    alleles_tsv_input = plot_step.input(os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.alleles.tsv.gz"))

    #"/ref/other/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz"
    #"/ref/other/gencode.v42.annotation.gtf.gz"

    # ./scripts/compute_overlap_with_other_catalogs_using_bedtools.sh ${output_prefix}.variants.bed.gz "step8:overlap:${STR_type}"
    # python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${output_prefix}.variants.tsv.gz  ${output_prefix}.variants.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}:  "
    # python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${output_prefix}.alleles.tsv.gz  ${output_prefix}.alleles.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}: "
    # python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/gencode.v42.annotation.sorted.gtf.gz  ${output_prefix}.variants.tsv.gz  ${output_prefix}.variants.${suffix}.tsv.gz "
    # python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/gencode.v42.annotation.sorted.gtf.gz  ${output_prefix}.alleles.tsv.gz   ${output_prefix}.alleles.${suffix}.tsv.gz "

    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 2 --width 8 --height 5.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 3 --width 30 --height 6.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 4 --width 16 --height 5.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_mutation_rates.py --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")


    """
    python3 -u /str-truth-set/figures_and_tables/plot_syndip_indel_size_distribution.py --output-dir ${output_dir}
    
    python3 -u /str-truth-set/figures_and_tables/plot_summary_stats.py --output-dir ${output_dir}
    python3 -u /str-truth-set/figures_and_tables/plot_gene_constraint_info.py --output-dir ${output_dir}
    
    python3 -u /str-truth-set/figures_and_tables/plot_mutation_rates.py --output-dir ${output_dir}
    """


combine_step = bp.new_step(
    f"concatenate {len(df)} tables",
    image=DOCKER_IMAGE,
    arg_suffix="step4",
    output_dir=os.path.join(args.output_dir, "combined"),
    depends_on=filter_steps,
)

combine_step.command("set -exuo pipefail")

input_tsvs = [
    combine_step.input(os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.variants.tsv.gz"))
    for _, row in df.iterrows()
]

concat_tsv_output_filename = f"concat.{len(df)}_samples.variants.tsv.gz"
combine_step.command(
    f"python3 /hprc_pipeline/scripts/concat_per_sample_tables.py -o {concat_tsv_output_filename} " +
    " ".join(i.local_path for i in input_tsvs))

joined_tsv_output_filename = f"joined.{len(df)}_samples.variants.tsv.gz"
joined_tsv_output_stats_filename = f"joined.{len(df)}_samples.variants.stats.tsv.gz"
combine_step.command(
    f"python3 /hprc_pipeline/scripts/join_per_sample_variant_tables.py "
    f"--output-stats-tsv {joined_tsv_output_stats_filename} "
    f"-o {joined_tsv_output_filename} " +
    " ".join(i.local_path for i in input_tsvs))

combine_step.output(concat_tsv_output_filename)
combine_step.output(joined_tsv_output_filename)
combine_step.output(joined_tsv_output_stats_filename)
bp.run()


#%%
