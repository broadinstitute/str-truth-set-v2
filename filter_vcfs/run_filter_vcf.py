"""Hail Batch pipeline for running dipcall on HPRC assemblies.

Relevant links:

HPRC assemblies
  https://projects.ensembl.org/hprc/

The design and construction of reference pangenome graphs with minigraph by Li et al. 2020
  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z

Increased mutation and gene conversion within human segmental duplications by Vollger et al.
  https://www.nature.com/articles/s41586-023-05895-y

"""
import collections
import hailtop.fs as hfs
import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize, Delocalize

STR_ANALYSIS_DOCKER_IMAGE = "weisburd/str-analysis@sha256:6e3128a78fe2125e5c12a999a64479084bdec2ed69ee95072c87a846d2a550ff"
FILTER_VCFS_DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:d327d8164fb2eb155967a55eb28e95ed5cb37375ec016d9767ac589b605bf133"


EXPANSION_HUNTER_LOCI_PER_RUN = 10_000 # if exclude_homopolymers else 100000
GANGSTR_LOCI_PER_RUN = 1_000_000
STRAGLR_LOCI_PER_RUN = 200_000

def create_filter_step(bp, row, suffix, output_dir, exclude_homopolymers=False, only_pure_repeats=False,
                       keep_loci_that_have_overlapping_variants=False, use_preemptibles=True):
    filter_step = bp.new_step(
        f"filter_vcf: {row.sample_id}" + (" (keep loci that have overlapping variants)" if keep_loci_that_have_overlapping_variants else ""),
        image=STR_ANALYSIS_DOCKER_IMAGE,
        arg_suffix="filter-step",
        preemptible=use_preemptibles,
        cpu=2,
        memory="highmem",
        output_dir=output_dir)

    hg38_fasta_input, _ = filter_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    dipcall_vcf_input, dipcall_high_confidence_regions_bed_input = filter_step.inputs(
        f"gs://str-truth-set-v2/dipcall_pipeline/{row.sample_id}/{row.sample_id}.dip.vcf.gz",
        f"gs://str-truth-set-v2/dipcall_pipeline/{row.sample_id}/{row.sample_id}.dip.bed.gz",
    )

    filter_step.command("set -exuo pipefail")

    filter_step.command(f"[ -s {dipcall_high_confidence_regions_bed_input} ] || exit 1")  # check that the bed file isn't emtpy

    filter_step.command(f"bedtools intersect -header -f 1 -wa -u \
            -a {dipcall_vcf_input}  \
            -b {dipcall_high_confidence_regions_bed_input} \
            | bgzip > {row.sample_id}.high_confidence_regions.vcf.gz")
    filter_step.command(f"tabix -f {row.sample_id}.high_confidence_regions.vcf.gz")

    min_repeat_unit_length = 2 if exclude_homopolymers else 1
    allow_interruptions = "no" if only_pure_repeats else "only-if-pure-repeats-not-found"
    keep_loci_that_have_overlapping_variants_arg = "--keep-loci-that-have-overlapping-variants " if keep_loci_that_have_overlapping_variants else ""
    filter_step.command(f"python3 -u -m str_analysis.filter_vcf_to_STR_variants \
            -R {hg38_fasta_input} \
            --allow-interruptions {allow_interruptions} \
            --write-bed-file \
            --min-str-length 9 \
            --min-str-repeats 3 \
            --min-repeat-unit-length {min_repeat_unit_length} \
            {keep_loci_that_have_overlapping_variants_arg} \
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


def create_annotate_steps(bp, row, suffix, output_dir, exclude_homopolymers=False, use_preemptibles=True, negative_loci=False):

    # TODO implement negative_loci arg
    annotate_steps = []
    for table_type in "variants", "alleles":

        annotate_step = bp.new_step(
            f"annotate {table_type}: {row.sample_id}",
            image=FILTER_VCFS_DOCKER_IMAGE,
            arg_suffix=f"annotate-{table_type}-step",
            cpu=4,
            preemptible=use_preemptibles,
            memory="highmem",
            output_dir=output_dir)
        bed_input = annotate_step.input(
            os.path.join(output_dir, f"{row.sample_id}{suffix}.variants.bed.gz"))

        for input_bed_path in [
            "gs://str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/hg38_ver13.adjusted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/hg38_ver17.adjusted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/hg38.hipstr_reference.adjusted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/known_disease_associated_STR_loci.GRCh38.bed.gz",
            "gs://str-truth-set/hg38/ref/other/GRCh38GenomicSuperDup.without_decoys.sorted.bed.gz",
            "gs://str-truth-set/hg38/ref/other/popstr_catalog_v2.bed.gz",
            "gs://str-truth-set/hg38/ref/other/trgt_repeat_catalog.hg38.reformatted_to_motif_only.bed.gz",
            "gs://str-truth-set/hg38/ref/other/adotto_tr_catalog_v1.2.bed.gz",
        ]:
            input_catalog, _ = annotate_step.inputs(input_bed_path, f"{input_bed_path}.tbi")

        for input_gtf in [
            "gs://str-truth-set/hg38/ref/other/gencode.v42.annotation.sorted.gtf.gz",
            "gs://str-truth-set/hg38/ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz",
        ]:
            input_gtf, _ = annotate_step.inputs(input_gtf, f"{input_gtf}.tbi")

        for i in range(6, 31, 3):
            annotate_step.inputs(
                f"gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_{i}bp.bed.gz",
                f"gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_{i}bp.bed.gz.tbi")
            annotate_step.inputs(
                f"gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_{i}bp.bed.gz",
                f"gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_{i}bp.bed.gz.tbi")

        tsv_input_table = annotate_step.input(
            os.path.join(output_dir, f"{row.sample_id}{suffix}.{table_type}.tsv.gz"))

        annotate_step.command("mkdir /ref")
        annotate_step.command(f"ln -s {input_catalog.local_dir} /ref/other")

        annotate_step.command("set -exuo pipefail")

        which_motifs = "--exclude-homopolymers" if exclude_homopolymers else "--all-motifs"
        annotate_step.command(f"""
            python3 -u /str-truth-set/scripts/compute_overlap_with_gene_models.py /ref/other/gencode.v42.annotation.sorted.gtf.gz {tsv_input_table} {row.sample_id}{suffix}.{table_type}.with_gencode_columns.tsv.gz |& tee {row.sample_id}{suffix}.overlap_with_gencode.{table_type}.log 
            mv {row.sample_id}{suffix}.{table_type}.with_gencode_columns.tsv.gz {row.sample_id}{suffix}.{table_type}.tsv.gz
            
            python3 -u /str-truth-set/scripts/compute_overlap_with_gene_models.py /ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz {row.sample_id}{suffix}.{table_type}.tsv.gz {row.sample_id}{suffix}.{table_type}.with_MANE_columns.tsv.gz |& tee {row.sample_id}{suffix}.overlap_with_MANE.{table_type}.log
            mv {row.sample_id}{suffix}.{table_type}.with_MANE_columns.tsv.gz {row.sample_id}{suffix}.{table_type}.tsv.gz

            python3 -u /str-truth-set/scripts/compute_overlap_with_other_catalogs.py --all-repeats {which_motifs} {row.sample_id}{suffix}.{table_type}.tsv.gz {row.sample_id}{suffix}.{table_type}.with_overlap_columns.tsv.gz |& tee {row.sample_id}{suffix}.overlap_with_other_catalogs.{table_type}.log
            mv {row.sample_id}{suffix}.{table_type}.with_overlap_columns.tsv.gz {row.sample_id}{suffix}.{table_type}.tsv.gz
            
            mv {row.sample_id}{suffix}.{table_type}.tsv.gz {row.sample_id}{suffix}.annotated.{table_type}.tsv.gz            
        """)
        annotate_step.command(f"/str-truth-set/scripts/compute_overlap_with_other_catalogs_using_bedtools.sh {bed_input}")

        annotate_step.output(f"{row.sample_id}{suffix}.annotated.{table_type}.tsv.gz")
        annotate_step.output(f"{row.sample_id}{suffix}.overlap_with_gencode.{table_type}.log")
        annotate_step.output(f"{row.sample_id}{suffix}.overlap_with_MANE.{table_type}.log")
        annotate_step.output(f"{row.sample_id}{suffix}.overlap_with_other_catalogs.{table_type}.log")
        annotate_steps.append(annotate_step)

    return annotate_steps


def create_variant_catalogs_step(bp, row, suffix, output_dir, exclude_homopolymers=False, output_negative_loci=False, use_preemptibles=True):
    variant_catalogs_step = bp.new_step(
        f"variant catalogs: {row.sample_id}",
        image=FILTER_VCFS_DOCKER_IMAGE,
        arg_suffix="variant-catalog-step",
        cpu=1,
        preemptible=use_preemptibles,
        output_dir=output_dir,
        delocalize_by=Delocalize.GSUTIL_COPY,
    )

    hg38_fasta_input, _ = variant_catalogs_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    high_confidence_regions_bed_input = variant_catalogs_step.input(
        os.path.join("gs://str-truth-set-v2/dipcall_pipeline", row.sample_id, f"{row.sample_id}.dip.bed.gz"))
    if output_negative_loci:
        dipcall_vcf_input = variant_catalogs_step.input(
            os.path.join("gs://str-truth-set-v2/dipcall_pipeline", row.sample_id, f"{row.sample_id}.dip.vcf.gz"))
    variants_tsv_input = variant_catalogs_step.input(
        os.path.join(output_dir, f"{row.sample_id}{suffix}.variants.tsv.gz"))
    variants_bed_input = variant_catalogs_step.input(
        os.path.join(output_dir, f"{row.sample_id}{suffix}.variants.bed.gz"))
    all_hg38_repeats_bed_input, _ = variant_catalogs_step.inputs(
        "gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_9bp.bed.gz",
        "gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_9bp.bed.gz.tbi")

    variant_catalogs_step.command("set -exuo pipefail")
    variant_catalogs_step.command(
        f"python3 -u /str-truth-set/tool_comparison/scripts/convert_truth_set_to_variant_catalogs.py "
        f"--ref-fasta {hg38_fasta_input} " +
        (f"--all-indels-vcf {dipcall_vcf_input} " if output_negative_loci else "--skip-negative-loci ") +
        f"--skip-popstr "
        f"--expansion-hunter-loci-per-run {EXPANSION_HUNTER_LOCI_PER_RUN} "
        f"--gangstr-loci-per-run {GANGSTR_LOCI_PER_RUN} "
        f"--straglr-loci-per-run {STRAGLR_LOCI_PER_RUN} "
        f"--high-confidence-regions-bed {high_confidence_regions_bed_input} "
        f"--all-hg38-repeats-bed {all_hg38_repeats_bed_input} "
        f"--output-filename-prefix {row.sample_id}{suffix} "
        f"--truth-set-bed {variants_bed_input} "
        f"{variants_tsv_input}")

    variant_catalogs_step.command(
        f"find tool_comparison -name '*.json'")
    variant_catalogs_step.command(
        f"find tool_comparison -name '*.bed*'")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/{row.sample_id}{suffix}.positive_loci.bed.gz")
    if output_negative_loci:
        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/{row.sample_id}{suffix}.negative_loci.bed.gz")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/expansion_hunter/{row.sample_id}{suffix}.positive_loci.EHv5.*.json")
    if output_negative_loci:
        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/expansion_hunter/{row.sample_id}{suffix}.negative_loci.EHv5.*.json")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/gangstr/{row.sample_id}{suffix}.positive_loci.GangSTR.*.bed")
    if output_negative_loci:
        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/gangstr/{row.sample_id}{suffix}.negative_loci.GangSTR.*.bed")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/hipstr/{row.sample_id}{suffix}.positive_loci.HipSTR.*.bed")
    if output_negative_loci:
        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/hipstr/{row.sample_id}{suffix}.negative_loci.HipSTR.*.bed")

    #for chrom in [*range(1, 23), "X"]:
    #    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/popstr/{row.sample_id}{suffix}.positive_loci.popSTR.chr{chrom}.markerInfo.gz")
    #    if output_negative_loci:
    #        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/popstr/{row.sample_id}{suffix}.negative_loci.popSTR.chr{chrom}.markerInfo.gz")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/trgt/{row.sample_id}{suffix}.positive_loci.TRGT_repeat_catalog.bed")
    if output_negative_loci:
        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/trgt/{row.sample_id}{suffix}.negative_loci.TRGT_repeat_catalog.bed")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/longtr/{row.sample_id}{suffix}.positive_loci.LongTR.*.bed")
    if output_negative_loci:
        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/longtr/{row.sample_id}{suffix}.negative_loci.LongTR.*.bed")

    variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/straglr/{row.sample_id}{suffix}.positive_loci.straglr_catalog.*.bed")
    if output_negative_loci:
        variant_catalogs_step.output(f"./tool_comparison/variant_catalogs/straglr/{row.sample_id}{suffix}.negative_loci.straglr_catalog.*.bed")

    return variant_catalogs_step


def create_table_for_tool_comparisons_step(bp, row, suffix, output_dir, exclude_homopolymers=False,
                                           output_negative_loci=False, use_preemptibles=True):
    table_for_tool_comparisons_step = bp.new_step(
        f"table for tool comparisons: {row.sample_id}",
        image=FILTER_VCFS_DOCKER_IMAGE,
        arg_suffix="table-for-tool-comparisons-step",
        cpu=1,
        memory="highmem",
        preemptible=use_preemptibles,
        output_dir=output_dir,
        delocalize_by=Delocalize.GSUTIL_COPY,
    )

    annotated_variants_table_input = table_for_tool_comparisons_step.input(
        os.path.join(output_dir, f"{row.sample_id}{suffix}.annotated.variants.tsv.gz"))

    if output_negative_loci:
        negative_loci_bed_input = table_for_tool_comparisons_step.input(
            os.path.join(output_dir, f"{row.sample_id}{suffix}.negative_loci.bed"))

    if not output_negative_loci:
        table_for_tool_comparisons_step.command(
            f"python3 -u /str-truth-set/tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py --output-dir . "
                f"{annotated_variants_table_input}")

        table_for_tool_comparisons_step.command("ls -lhtr")

        table_for_tool_comparisons_step.output(
            annotated_variants_table_input.filename.replace(".tsv", ".for_comparison.tsv"))
    else:
        # compute overlap with other catalogs
        #variant_catalog_step.command(
        #    f"python3 -u /str-truth-set/scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs {row.sample_id}{suffix}.negative_loci.bed "
        #    f"{row.sample_id}{suffix}.negative_loci.with_overlap_columns.tsv.gz")
        #table_for_tool_comparisons_step.command(
        #    f"python3 -u /str-truth-set/tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py --output-dir . "
        #    f"{annotated_variants_table_input} "
        #    f"{row.sample_id}{suffix}.negative_loci.with_overlap_columns.tsv.gz")

        table_for_tool_comparisons_step.command("ls -lhtr")

    return table_for_tool_comparisons_step

def create_plot_step(bp, suffix, output_dir, row=None, alleles_tsv_step=None, exclude_homopolymers=False):
    """Must specify either row or alleles_tsv_step"""

    plot_step = bp.new_step(
        f"plots: {row.sample_id}" if row is not None else f"plots: combined",
        image=FILTER_VCFS_DOCKER_IMAGE,
        cpu=1 if row is not None else 16,
        memory="highmem",
        arg_suffix="plot-step",
        output_dir=output_dir,
        #localize_by=Localize.GSUTIL_COPY,
    )

    plot_step.command("set -exuo pipefail")
    #plot_step.switch_gcloud_auth_to_user_account()
    if row is not None:
        alleles_tsv_input = plot_step.input(os.path.join(output_dir, f"{row.sample_id}{suffix}.annotated.alleles.tsv.gz"))
    elif alleles_tsv_step is not None:
        alleles_tsv_input = plot_step.use_previous_step_outputs_as_inputs(alleles_tsv_step)
    else:
        raise ValueError("Must specify row or alleles_tsv_step")

    constraint_table_input = plot_step.input(f"gs://str-truth-set/hg38/ref/other/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
    disease_associated_loci_tsv_input = plot_step.input(f"gs://str-truth-set/hg38/ref/other/known_disease_associated_STR_loci.tsv")
    repeat_spec_for_motif_distribution_plot_input = plot_step.input("gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_12bp.bed.gz")
    ref_dir_local_path = os.path.dirname(repeat_spec_for_motif_distribution_plot_input.local_dir)

    # figure 1 panels
    dipcall_vcf_input = None
    if row is not None:
        dipcall_vcf_input = plot_step.input(f"gs://str-truth-set-v2/dipcall_pipeline/{row.sample_id}/{row.sample_id}.dip.vcf.gz")
        plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_syndip_indel_size_distribution.py --width 12 --height 4 --syndip-vcf {dipcall_vcf_input} --image-type png")

    # figure 2 panels
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 2 --width 8 --height 5.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 3 --width 30 --height 6.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png --only-pure-repeats")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 4 --width 16 --height 5.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png --only-pure-repeats")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py  --only-plot 8 --width 50 --height 6.5  --truth-set-alleles-table {alleles_tsv_input} --image-type png --only-pure-repeats")

    # figure 3 panels
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py --only-plot 6 --width 11 --height 12 --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_motif_distribution.py --width 9 --height 8  --truth-set-alleles-table {alleles_tsv_input} --ref-dir {ref_dir_local_path} --image-type png")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_motif_distribution.py --width 9 --height 8  --only-pure-repeats --truth-set-alleles-table {alleles_tsv_input} --ref-dir {ref_dir_local_path} --image-type png")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_motif_distribution.py --width 9 --height 8  --only-pure-repeats --truth-set-alleles-table {alleles_tsv_input} --ref-dir {ref_dir_local_path} --image-type png --include-homopolymers")

    # figure 7 panels
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_summary_stats.py --image-type png --only-plot 5 --only-pure-repeats --width 8 --height 7 --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_gene_constraint_info.py --constraint-table {constraint_table_input} --disease-associated-loci-table {disease_associated_loci_tsv_input} --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command(f"python3 /str-truth-set/figures_and_tables/plot_mutation_rates.py --truth-set-alleles-table {alleles_tsv_input} --image-type png")

    # supp. figure 2
    #plot_step.command(f"python3 /str-truth-set/figures_and_tables/paper_generate_table3_and_supp_fig2_repeat_catalogs.py  --truth-set-alleles-table {alleles_tsv_input} --image-type png")
    plot_step.command("ls -l")

    files_to_download = {}
    for filename in [
        "allele_size_distribution_by_number_of_repeats.color_by_interruptions.png",
        "allele_size_distribution_by_number_of_repeats.2-6bp_motifs.color_by_interruptions.png",
        "allele_size_distribution_by_number_of_repeats.7-24bp_motifs.color_by_interruptions.png",
        "allele_size_distribution_by_number_of_repeats.25-50bp_motifs.color_by_interruptions.png",
        "allele_size_distribution_by_number_of_repeats_and_motif_size.only_pure_repeats.color_by_multiallelic.png",
        "allele_size_distribution_by_number_of_repeats_and_motif_size.only_pure_repeats.2-6bp_motifs.color_by_multiallelic.png",
        "allele_size_distribution_by_number_of_repeats_and_motif_size.only_pure_repeats.7-24bp_motifs.color_by_multiallelic.png",
        "allele_size_distribution_by_number_of_repeats_and_motif_size.only_pure_repeats.25-50bp_motifs.color_by_multiallelic.png",
        "allele_size_distribution_by_number_of_repeats_and_motif_size.only_pure_repeats.color_by_overlapssegdupintervals.png",
        "gene_constraint_metrics_and_STRs.png",
        "gene_constraint_metrics_vs_STR_allele_size.png",
        "motif_distribution.png",
        "motif_distribution.pure_repeats.png",
        "motif_distribution.including_homopolymers.pure_repeats.png",
        "motif_distribution_in_hg38_with_atleast_12bp.png",
        "mutation_rates_by_allele_size.2bp_motifs.png",
        "mutation_rates_by_allele_size.3-6bp_motifs.png",
        "mutation_rates_by_fraction_interrupted_repeats.png",
        "reference_locus_size_distribution.2_to_6bp_motifs.png",
        "reference_locus_size_distribution.7_to_24bp_motifs.png",
        "reference_locus_size_distribution.25_to_50bp_motifs.png",
        "reference_locus_size_distribution.2bp_motifs.png",
        "reference_locus_size_distribution.3bp_motifs.png",
        "reference_locus_size_distribution.4bp_motifs.png",
        "reference_locus_size_distribution.5bp_motifs.png",
        "reference_locus_size_distribution.6bp_motifs.png",
        "reference_locus_size_distribution.png",
        "truth_set_gene_overlap.only_pure_repeats.all_regions.MANE_v1.png",
        "truth_set_gene_overlap.only_pure_repeats.all_regions.gencode_v42.png",
        "truth_set_gene_overlap.only_pure_repeats.excluding_introns_and_intergenic.MANE_v1.png",
        "truth_set_gene_overlap.only_pure_repeats.excluding_introns_and_intergenic.gencode_v42.png",
    ] + ([
        "allele_size_distribution_by_number_of_repeats.x3.only_pure_repeats.including_homopolymers.png",
        "allele_size_distribution_by_repeat_size_in_base_pairs.x3.only_pure_repeats.including_homopolymers.png",
        "STR_allele_size_distribution_by_number_of_repeats.x6.only_pure_repeats.including_homopolymers.png",
        "STR_allele_size_distributions_by_repeat_size_in_base_pairs.x6.only_pure_repeats.including_homopolymers.png",
        "allele_size_distribution_by_number_of_repeats.1bp_motifs.color_by_interruptions.png",
        "allele_size_distribution_by_number_of_repeats_and_motif_size.only_pure_repeats.1bp_motifs.color_by_multiallelic.png",
        "mutation_rates_by_allele_size.1bp_and_2bp_motifs.png",
        "reference_locus_size_distribution.1bp_motifs.png",
    ] if not exclude_homopolymers else [
        "allele_size_distribution_by_number_of_repeats.x3.only_pure_repeats.png",
        "allele_size_distribution_by_repeat_size_in_base_pairs.x3.only_pure_repeats.png",
        "STR_allele_size_distribution_by_number_of_repeats.x6.only_pure_repeats.png",
        "STR_allele_size_distributions_by_repeat_size_in_base_pairs.x6.only_pure_repeats.png",
    ]) + ([
        "syndip_indel_size_distribution.png",
    ] if dipcall_vcf_input is not None else []):
        if row is not None:
            output_path = os.path.join(output_dir, "figures", f"{row.sample_id}{suffix}.{filename}")
        else:
            output_path = os.path.join(output_dir, "figures", f"combined{suffix}.{filename}")
        plot_step.output(filename, output_path)
        files_to_download[filename.replace(".png", "")] = output_path

    return plot_step, files_to_download


def create_combine_results_step(
        bp, df, suffix, filter_steps,
        annotate_variants_steps, annotate_alleles_steps, output_dir, exclude_homopolymers=False,
        keep_loci_that_have_overlapping_variants=False,
        use_preemptibles=True):
    combine_steps = []
    combined_output_dir = os.path.join(output_dir, "combined")

    for combine_step_type, combine_step_prefix, variants_or_alleles, combine_step_suffix, cpu in [
        ("concat tsvs", "concat", "variants", "tsv", 4),
        ("concat tsvs", "concat", "annotated.variants", "tsv", 8),
        ("concat tsvs", "concat", "alleles", "tsv", 4),
        ("concat tsvs", "concat", "annotated.alleles", "tsv", 8),
        ("join tsvs", "joined", "variants", "tsv", 2),
        ("combine beds", "combined", "variants", "bed", 1),
    ]:
        if variants_or_alleles == "annotated.variants" and annotate_variants_steps is None:
            combine_steps.append(None)
            continue
        if variants_or_alleles == "annotated.alleles" and annotate_alleles_steps is None:
            combine_steps.append(None)
            continue

        combine_step = bp.new_step(
            f"{combine_step_type}: {variants_or_alleles}: {len(df)} samples",
            arg_suffix="combine-step",
            image=FILTER_VCFS_DOCKER_IMAGE,
            cpu=cpu,
            preemptible=use_preemptibles,
            memory="highmem",
            output_dir=combined_output_dir,
        )

        if variants_or_alleles == "annotated.variants":
            combine_step.depends_on(annotate_variants_steps)
        elif variants_or_alleles == "annotated.alleles":
            combine_step.depends_on(annotate_alleles_steps)
        else:
            combine_step.depends_on(filter_steps)

        combine_step.command("set -exuo pipefail")

        input_files = [
            combine_step.input(
                os.path.join(output_dir, f"{row.sample_id}/{row.sample_id}{suffix}.{variants_or_alleles}.{combine_step_suffix}.gz")
            ) for _, row in df.iterrows()
        ]

        if combine_step_type == "concat tsvs":
            concat_tsv_output_filename = f"{combine_step_prefix}.{len(df)}_samples.{variants_or_alleles}.{combine_step_suffix}.gz"
            if not exclude_homopolymers:
                combine_step.storage("20G")
            combine_step.command(
                f"python3 /filter_vcfs/scripts/concat_per_sample_tables.py -o {concat_tsv_output_filename} " +
                " ".join(i.local_path for i in input_files))
            combine_step.output(concat_tsv_output_filename)
        elif combine_step_type == "join tsvs":
            joined_tsv_output_filename = f"joined.{len(df)}_samples.variants.tsv.gz"
            joined_tsv_output_stats_filename = f"joined.{len(df)}_samples.variants.stats.tsv.gz"
            combine_step.command(
                f"python3 /filter_vcfs/scripts/join_per_sample_variant_tables.py "
                f"--output-stats-tsv {joined_tsv_output_stats_filename} "
                f"-o {joined_tsv_output_filename} " +
                " ".join(i.local_path for i in input_files))
            combine_step.output(joined_tsv_output_filename)
            combine_step.output(joined_tsv_output_stats_filename)
        elif combine_step_type == "combine beds":
            combined_bed_output_filename = f"combined.{len(df)}_samples.variants.bed.gz"
            combined_bed_output_stats_filename = f"combined.{len(df)}_bed_files.variants.stats.tsv.gz"
            combine_step.command(
                f"python3 /filter_vcfs/scripts/combine_per_sample_variant_bed_files.py "
                f"--output-stats-tsv {combined_bed_output_stats_filename} "
                f"-o {combined_bed_output_filename} " +
                " ".join(i.local_path for i in input_files))
            combine_step.output(combined_bed_output_filename)
            combine_step.output(combined_bed_output_filename + ".tbi")
            combine_step.output(combined_bed_output_stats_filename)

        combine_steps.append(combine_step)

    _, concat_annotated_variants_tsv_step, _, concat_annotated_alleles_tsv_step, join_tsvs_step, _ = combine_steps

    # create plots
    figures_to_download_dict = {}
    if concat_annotated_alleles_tsv_step is not None:
        combine_plot_step, figures_to_download_dict = create_plot_step(
            bp, suffix, combined_output_dir,
            alleles_tsv_step=concat_annotated_alleles_tsv_step,
            exclude_homopolymers=exclude_homopolymers)

    # convert to catalogs
    combined_variant_catalogs_step = bp.new_step(
        f"combined variant catalogs: {len(df)} samples",
        arg_suffix="combine-variant-catalogs-step",
        cpu=2,
        preemptible=use_preemptibles,
        memory="highmem",
        storage="20G",
        image=FILTER_VCFS_DOCKER_IMAGE)

    hg38_fasta_input, _ = combined_variant_catalogs_step.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    joined_tsv_input, _ = combined_variant_catalogs_step.use_previous_step_outputs_as_inputs(join_tsvs_step)

    combined_variant_catalogs_step.command("set -exuo pipefail")
    combined_variant_catalogs_step.command(
        "python3 -u /str-truth-set/tool_comparison/scripts/convert_truth_set_to_variant_catalogs.py "
        f"--ref-fasta {hg38_fasta_input} "
        f"--expansion-hunter-loci-per-run {EXPANSION_HUNTER_LOCI_PER_RUN} "
        f"--gangstr-loci-per-run {GANGSTR_LOCI_PER_RUN} "
        f"--straglr-loci-per-run {STRAGLR_LOCI_PER_RUN} "
        "--skip-negative-loci "
        "--skip-popstr "
        f"{joined_tsv_input}")

    combined_variant_catalogs_step.command(
        f"find tool_comparison -name '*.json'")
    combined_variant_catalogs_step.command(
        f"find tool_comparison -name '*.bed'")
    combined_variant_catalogs_step.output(
        f"./tool_comparison/variant_catalogs/expansion_hunter/positive_loci.EHv5.001_of_001.json",
        os.path.join(output_dir, "combined", f"combined.{len(df)}_samples.positive_loci.json")
    )

    return figures_to_download_dict


def main():
    bp = pipeline("filter_vcf_to_STRs", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--only-pure-repeats", action="store_true")
    parser.add_argument("--exclude-homopolymers", action="store_true")
    parser.add_argument("--output-negative-loci", action="store_true")
    parser.add_argument("--skip-combine-steps", action="store_true")
    parser.add_argument("--use-nonpreemptibles", action="store_true")
    parser.add_argument("-s", "--sample-id", action="append",
                        help="Process only this sample. Can be specified more than once.")
    parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf")
    args = bp.parse_known_args()

    bp.precache_file_paths("gs://str-truth-set-v2/filter_vcf/**/*.*")

    df = pd.read_table("../dipcall_pipeline/hprc_assemblies.tsv")
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

    filter_steps_keeping_all_loci = []
    filter_steps = []
    annotate_alleles_steps = []
    annotate_variants_steps = []
    figures_to_download = collections.defaultdict(list)
    for row_i, (_, row) in enumerate(df.iterrows()):
        output_dir = os.path.join(args.output_dir, output_dir_suffix, row.sample_id)

        # The default filter options discard detected STR variants if they overlap any other variant in the input VCF
        # (ie. SNVs or non-STR indels). This is because overlapping variants make the true STR genotype ambiguous,
        # and so are not suitable for inclusion in an STR truth set. However, they are still polymorphic STR loci, and
        # so can be included in a set of all polymorphic STRs in the genome.
        # The steps below first compute the set of all STRs without this post-filter (for use as a polymorphic STR catalog)
        # and then compute the post-filtered set of STRs without any overlapping variants (for use as an STR truth set).
        filter_step_keeping_all_loci = create_filter_step(bp, row, f"{suffix}.keeping_loci_that_have_overlapping_variants",
                                            os.path.join(args.output_dir, f"{output_dir_suffix}_keeping_loci_that_have_overlapping_variants", row.sample_id),
                                            exclude_homopolymers=args.exclude_homopolymers,
                                            only_pure_repeats=args.only_pure_repeats,
                                            keep_loci_that_have_overlapping_variants=True,
                                            use_preemptibles=not args.use_nonpreemptibles)

        filter_steps_keeping_all_loci.append(filter_step_keeping_all_loci)


        filter_step = create_filter_step(bp, row, suffix, output_dir,
                                         exclude_homopolymers=args.exclude_homopolymers,
                                         only_pure_repeats=args.only_pure_repeats,
                                         use_preemptibles=not args.use_nonpreemptibles)

        filter_steps.append(filter_step)

        annotate_variants_step, annotate_alleles_step = create_annotate_steps(bp, row, suffix, output_dir,
                                                                              exclude_homopolymers=args.exclude_homopolymers,
                                                                              use_preemptibles=not args.use_nonpreemptibles,
                                                                              negative_loci=False)
        annotate_variants_step.depends_on(filter_step)
        annotate_alleles_step.depends_on(filter_step)
        annotate_variants_steps.append(annotate_variants_step)
        annotate_alleles_steps.append(annotate_alleles_step)

        variant_catalogs_step = create_variant_catalogs_step(bp, row, suffix, output_dir,
                                                             exclude_homopolymers=args.exclude_homopolymers,
                                                             output_negative_loci=args.output_negative_loci,
                                                             use_preemptibles=not args.use_nonpreemptibles)
        variant_catalogs_step.depends_on(filter_step)

        #TODO: also annotate the negative loci
        #if args.output_negative_loci:
        #    annotate_variants_step, annotate_alleles_step = create_annotate_steps(
        #        bp, row, suffix, output_dir,
        #        exclude_homopolymers=args.exclude_homopolymers,
        #        use_preemptibles=not args.use_nonpreemptibles,
        #        negative_loci=True)

        if row.sample_id in ("HG002", "CHM1_CHM13"):
            table_for_tool_comparisons_step = create_table_for_tool_comparisons_step(bp, row, suffix, output_dir,
                                                                            exclude_homopolymers=args.exclude_homopolymers,
                                                                            output_negative_loci=False,  #TODO: args.output_negative_loci,
                                                                            use_preemptibles=not args.use_nonpreemptibles)
            table_for_tool_comparisons_step.depends_on(annotate_alleles_step)

            plot_step, figures_to_download_dict = create_plot_step(bp, suffix, output_dir, row=row,
                                                                   exclude_homopolymers=args.exclude_homopolymers)

            plot_step.depends_on(annotate_alleles_step)

            for label, file_path in figures_to_download_dict.items():
                figures_to_download[label].append(file_path)


    if not args.skip_combine_steps:
        figures_to_download_dict = create_combine_results_step(
            bp, df, suffix,
            filter_steps, annotate_variants_steps, annotate_alleles_steps,
            output_dir=os.path.join(args.output_dir, output_dir_suffix),
            exclude_homopolymers=args.exclude_homopolymers,
            use_preemptibles=not args.use_nonpreemptibles)

        for label, file_path in figures_to_download_dict.items():
            figures_to_download[label].append(file_path)

        create_combine_results_step(
            bp, df, f"{suffix}.keeping_loci_that_have_overlapping_variants",
            filter_steps_keeping_all_loci, None, None,
            output_dir=os.path.join(args.output_dir, f"{output_dir_suffix}_keeping_loci_that_have_overlapping_variants"),
            exclude_homopolymers=args.exclude_homopolymers,
            use_preemptibles=not args.use_nonpreemptibles)

    bp.run()

    # download figures
    for label, file_paths in figures_to_download.items():
        local_dir = "results_without_homopolymers" if args.exclude_homopolymers else "results_with_homopolymers"
        local_dir = f"../{local_dir}/figures/{label}"
        if not os.path.exists(local_dir):
            os.system(f"mkdir -p {local_dir}")
        for file_path in file_paths:
            local_path = os.path.join(local_dir, os.path.basename(file_path))
            if not os.path.exists(local_path):
                print(f"Downloading {file_path} to {local_dir}")
                try:
                    hfs.copy(file_path, local_path)
                except Exception as e:
                    print(e)


if __name__ == "__main__":
    main()


#%%
