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

STR_ANALYSIS_DOCKER_IMAGE = "weisburd/str-analysis@sha256:e13cf6e945bf04f1fbfbe1da880f543a7bb223026e995b2682324cebc8c18649"
DOCKER_IMAGE = "weisburd/hprc-pipeline@sha256:89d59d62ce0d3c1129207bda0246501735e60a471495f9aa86311441201a0105"


def create_filter_step(bp, row, args, suffix):
    filter_step = bp.new_step(
        f"filter_vcf: {row.sample_id}",
        image=STR_ANALYSIS_DOCKER_IMAGE,
        arg_suffix="filter-step",
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
        f"gs://str-truth-set-v2/hprc_dipcall/{row.sample_id}/{row.sample_id}.dip.bed.gz",
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

    return filter_step


def create_annotate_steps(bp, row, args, suffix):

    annotate_steps = []
    for table_type in "variants", "alleles":

        annotate_step = bp.new_step(
            f"annotate {table_type}: {row.sample_id}",
            image=DOCKER_IMAGE,
            arg_suffix=f"annotate-{table_type}-step",
            cpu=4,
            memory="standard",
            output_dir=os.path.join(args.output_dir, row.sample_id),
        )
        bed_input = annotate_step.input(
            os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.variants.bed.gz"))

        for input_bed_path in [
            "gs://str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/hg38_ver13.adjusted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/hg38_ver17.adjusted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/hg38.hipstr_reference.adjusted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/known_disease_associated_STR_loci.GRCh38.bed.gz",
            "gs://str-truth-set/hg38/ref/other/GRCh38GenomicSuperDup.without_decoys.sorted.bed.gz",
        ]:
            input_catalog, _ = annotate_step.inputs(input_bed_path, f"{input_bed_path}.tbi")

        for input_gtf in [
            "gs://str-truth-set/hg38/ref/other/gencode.v43.annotation.gtf.gz",
            "gs://str-truth-set/hg38/ref/other/MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz",
        ]:
            input_gtf, _ = annotate_step.inputs(input_gtf, f"{input_gtf}.tbi")

        for i in range(6, 31, 3):
            annotate_step.inputs(
                f"gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_{i}bp.bed.gz",
                f"gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_{i}bp.bed.gz.tbi")

        tsv_input_table = annotate_step.input(
            os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.{table_type}.tsv.gz"))

        annotate_step.command("mkdir /ref")
        annotate_step.command(f"ln -s {input_catalog.local_dir} /ref/other")

        annotate_step.command("set -exuo pipefail")
        annotate_step.command(f"/hprc_pipeline/scripts/compute_overlap_with_other_catalogs_using_bedtools.sh {bed_input}")

        annotate_step.command(f"""
            python3 -u /hprc_pipeline/scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs {tsv_input_table}  {row.sample_id}{suffix}.{table_type}.with_overlap_columns.tsv.gz     
            mv {row.sample_id}{suffix}.{table_type}.with_overlap_columns.tsv.gz {row.sample_id}{suffix}.{table_type}.tsv.gz
            
            python3 -u /hprc_pipeline/scripts/compute_overlap_with_gene_models.py /ref/other/gencode.v43.annotation.sorted.gtf.gz {row.sample_id}{suffix}.{table_type}.tsv.gz  {row.sample_id}{suffix}.variants.with_gencode_columns.tsv.gz
            mv {row.sample_id}{suffix}.{table_type}.with_gencode_columns.tsv.gz {row.sample_id}{suffix}.{table_type}.tsv.gz
            
            python3 -u /hprc_pipeline/scripts/compute_overlap_with_gene_models.py /ref/other/MANE.v1.2.ensembl_genomic.sorted.gtf.gz {row.sample_id}{suffix}.{table_type}.tsv.gz {row.sample_id}{suffix}.variants.with_MANE_columns.tsv.gz
            mv {row.sample_id}{suffix}.{table_type}.with_MANE_columns.tsv.gz {row.sample_id}{suffix}.annotated.{table_type}.tsv.gz
            """)

        annotate_step.output(f"{row.sample_id}{suffix}.annotated.{table_type}.tsv.gz")

        annotate_steps.append(annotate_step)

    return annotate_steps


def create_variant_catalogs_step(bp, row, args, suffix):
    variant_catalogs_step = bp.new_step(
        f"variant catalogs: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="variant-catalogs-step",
        cpu=1,
        output_dir=os.path.join(args.output_dir, row.sample_id),
    )

    high_confidence_regions_bed_input = variant_catalogs_step.input(
        os.path.join("gs://str-truth-set-v2/hprc_dipcall", row.sample_id, f"{row.sample_id}.dip.bed.gz"))
    variants_tsv_input = variant_catalogs_step.input(
        os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.variants.tsv.gz"))
    variants_bed_input = variant_catalogs_step.input(
        os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.variants.bed.gz"))
    all_hg38_repeats_bed_input, _ = variant_catalogs_step.inputs(
        "gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_9bp.bed.gz",
        "gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_9bp.bed.gz.tbi")

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

    variant_catalogs_step.command(f"ls -lh ./tool_comparison/variant_catalogs/expansion_hunter/positive_loci.EHv5.001_of_001.json")
    variant_catalogs_step.command(f"ls -lh ./tool_comparison/variant_catalogs/expansion_hunter/negative_loci.EHv5.001_of_001.json")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/expansion_hunter/positive_loci.EHv5.001_of_001.json")
    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/expansion_hunter/negative_loci.EHv5.001_of_001.json")
    return variant_catalogs_step


def create_plot_step(bp, row, args, suffix):
    # generate plots
    plot_step = bp.new_step(
        f"plots: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="plots-step",
        output_dir=os.path.join(args.output_dir, row.sample_id),
    )

    plot_step.command("set -exuo pipefail")

    alleles_tsv_input = plot_step.input(os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.annotated.alleles.tsv.gz"))

    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 2 --width 8 --height 5.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 3 --width 30 --height 6.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 4 --width 16 --height 5.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_mutation_rates.py --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")

    return plot_step


def create_combine_results_step(bp, df, args, suffix, filter_steps):
    for combine_step_type, combine_step_prefix, combine_step_suffix in [
        ("concat tsvs", "concat", "tsv"),
        ("join tsvs", "joined", "tsv"),
        ("combine beds", "combined", "bed"),
    ]:
        combine_step = bp.new_step(
            f"{combine_step_type}: {len(df)} samples",
            arg_suffix="combine-step",
            image=DOCKER_IMAGE,
            cpu=2,
            memory="highmem",
            output_dir=os.path.join(args.output_dir, "combined"),
            depends_on=filter_steps,
        )

        combine_step.command("set -exuo pipefail")

        input_files = [
            combine_step.input(
                os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}{suffix}.variants.{combine_step_suffix}.gz")
            ) for _, row in df.iterrows()
        ]

        if combine_step_type == "concat tsvs":
            concat_tsv_output_filename = f"{combine_step_prefix}.{len(df)}_samples.variants.{combine_step_suffix}.gz"
            combine_step.command(
                f"python3 /hprc_pipeline/scripts/concat_per_sample_tables.py -o {concat_tsv_output_filename} " +
                " ".join(i.local_path for i in input_files))
            combine_step.output(concat_tsv_output_filename)
        elif combine_step_type == "join tsvs":
            joined_tsv_output_filename = f"joined.{len(df)}_samples.variants.tsv.gz"
            joined_tsv_output_stats_filename = f"joined.{len(df)}_samples.variants.stats.tsv.gz"
            combine_step.command(
                f"python3 /hprc_pipeline/scripts/join_per_sample_variant_tables.py "
                f"--output-stats-tsv {joined_tsv_output_stats_filename} "
                f"-o {joined_tsv_output_filename} " +
                " ".join(i.local_path for i in input_files))
            combine_step.output(joined_tsv_output_filename)
            combine_step.output(joined_tsv_output_stats_filename)
        elif combine_step_type == "combine beds":
            combined_bed_output_filename = f"combined.{len(df)}_samples.variants.bed.gz"
            combined_bed_output_stats_filename = f"combined.{len(df)}_bed_files.variants.stats.tsv.gz"
            combine_step.command(
                f"python3 /hprc_pipeline/scripts/combine_per_sample_variant_bed_files.py "
                f"--output-stats-tsv {combined_bed_output_stats_filename} "
                f"-o {combined_bed_output_filename} " +
                " ".join(i.local_path for i in input_files))
            combine_step.output(combined_bed_output_filename)
            combine_step.output(combined_bed_output_stats_filename)

    return combine_step


def main():
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
    for row_i, (_, row) in enumerate(df.iterrows()):

        filter_step = create_filter_step(bp, row, args, suffix)
        filter_steps.append(filter_step)

        annotate_variants_step, annotate_alleles_step = create_annotate_steps(bp, row, args, suffix)
        annotate_variants_step.depends_on(filter_step)
        annotate_alleles_step.depends_on(filter_step)

        variant_catalogs_step = create_variant_catalogs_step(bp, row, args, suffix)
        variant_catalogs_step.depends_on(filter_step)

        plot_step = create_plot_step(bp, row, args, suffix)
        plot_step.depends_on(annotate_alleles_step)

    create_combine_results_step(bp, df, args, suffix, filter_steps)

    bp.run()


if __name__ == "__main__":
    main()


#%%
