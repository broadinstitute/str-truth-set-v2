"""
This pipeline runs TR genotyping tools on the catalog(s) of interest.

Per-sample inputs:
    short read or long read cram/crai file
    filter_vcf output catalogs for the tools of interest
    output directory
"""

import collections
import json
import hailtop.fs as hfs
import os
import pandas as pd
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize
import sys
import tempfile

sys.path.append("../str-truth-set/tool_comparison/hail_batch_pipelines")
from expansion_hunter_pipeline import create_expansion_hunter_steps, DOCKER_IMAGE as EH_DOCKER_IMAGE
from gangstr_pipeline import create_gangstr_steps
from hipstr_pipeline import create_hipstr_steps
from trgt_pipeline import create_trgt_step, DOCKER_IMAGE as TRGT_V5_DOCKER_IMAGE, TRGT_V3_DOCKER_IMAGE
from longtr_pipeline import create_longtr_steps
from inquistr_pipeline import create_inquistr_steps
from vamos_pipeline import create_vamos_step
from atarva_pipeline import create_atarva_step
from ensembletr_pipeline import create_ensembletr_steps

# EnsembleTR is a consensus/merge tool run in two modes (each a separate "tool" in the comparison). It consumes the
# already-computed per-caller outputs for the same (sample, data_type, coverage): the ExpansionHunter json (from
# ENSEMBLETR_EH_SOURCE_TOOL) plus the HipSTR (+GangSTR) native VCFs.
ENSEMBLETR_EH_SOURCE_TOOL = "EHv5-bw2-optimized"
ENSEMBLETR_TOOLS = {
    "EnsembleTR-EH+HipSTR",
    "EnsembleTR-EH+HipSTR+GangSTR",
}

SHORT_READ_TOOLS = {
    "IlluminaEHv5",
    "EHv5",
    "EHv5-bw2-optimized",
    "GangSTR",
    "HipSTR",
} | ENSEMBLETR_TOOLS

LONG_READ_TOOLS = {
    "TRGTv3",
    "TRGTv5",
    "LongTR",
    "inquiSTR",
    "vamos",
    "ATaRVa",
}

# Tools whose output VCF carries an actual allele sequence in REF/ALT, so their sequence accuracy can be scored by
# edit distance against the assembly truth allele sequences (see create_extract_allele_sequences_step). Every other
# tool reports only a repeat count, which the accuracy plots already cover. This is the single place that decides
# which tools get the extra extract step, the extra add-columns invocation, and the extra plot loop.
SEQUENCE_ACCURACY_TOOLS = {
    "TRGTv5",
    "ATaRVa",
    "HipSTR",
}

# The add-columns and plot steps use the /str-truth-set baked into FILTER_VCFS_DOCKER_IMAGE (rebuild that image and
# update its digest to pick up new str-truth-set scripts).

# Motif size bins (min, max) used to stratify the accuracy plots.
MOTIF_SIZE_BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (2, 6), (7, 24), (25, 1000)]

# IlluminaEHv5 (official Illumina EH v5 build) crashes with "numIndels out of range" on large loci, so its
# prefiltered catalog drops any locus whose reference interval is at least this many base pairs wide.
ILLUMINA_EH_MAX_LOCUS_SIZE_BP = 500

SHORT_READ_DATA_TYPES = {
    "illumina",
    "illumina_exome",
    "element",
    "ultima",
    "illumina_rnaseq",
}

LONG_READ_DATA_TYPES = {
    "pacbio",
    "ONT",
    "pacbio_isoseq",
}

# RNA-seq data types are labeled by total bases sequenced (Gbp) rather than genome coverage, since
# genome-wide depth is meaningless for transcript data (reads only cover expressed loci).
RNASEQ_DATA_TYPES = {
    "illumina_rnaseq",
    "pacbio_isoseq",
}

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

FILTER_VCFS_DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:025432ff71ee297d21f4d72dafd898c45ec35bd54512f65b45ea7dc57dc42062"

# Image for run_tools scripts run as Hail Batch steps (built by .github/workflows/build_run_tools_image.yml from
# run_tools/docker/Dockerfile). Used by the per-sample build-catalogs step, which runs
# convert_truth_set_to_variant_catalogs.py baked into the image.
RUN_TOOLS_DOCKER_IMAGE = "weisburd/run-tools@sha256:311ee8747a37235ab13da3024733b44006a80bd99806422766a8dcd07e5dcc10"

DEFAULT_OUTPUT_DIR = "gs://str-truth-set-v2/tool_results"

# --benchmark-resources mode (see run_resource_benchmark): benchmark one tool's runtime/memory/cost vs. catalog
# size, by subsampling the TRExplorer catalog (trexplorer.broadinstitute.org) to the N most-polymorphic loci.
# TREXPLORER_BIGQUERY_TABLE is the live BigQuery snapshot; HPRC256_Stdev is the per-locus standard deviation of allele
# sizes across the 256 HPRC samples, used to rank loci by polymorphism (highest stdev = most polymorphic). Update
# the table id when TRExplorer reloads it (it's the TABLE_ID const in tandem-repeat-explorer/website/index.html).
TREXPLORER_BIGQUERY_TABLE = "cmg-analysis.tandem_repeat_explorer.catalog_20260712_002034"
BENCHMARK_TOOL = "EHv5-bw2-optimized"
# Full range of catalog sizes (# loci) — this is the range behind the resource_metrics.json files the viewer reads.
# Pass a smaller subset via --benchmark-catalog-sizes for a quick test run (the sizes are nested, so a subset of
# these reuses the already-uploaded catalogs).
BENCHMARK_CATALOG_SIZES = [2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000]

# Per-tool resource-benchmark config. "cpu" is the VM's provisioned cores (used for VM sizing + the recorded VM
# metadata); "threads" is how many cores the tool actually uses (Runtime CPU-hours = wall x threads). They can differ
# (a single-threaded tool provisioned with extra cores just for RAM), so runtime never counts idle cores. GangSTR and
# HipSTR run cpu=1 single-threaded on a bed catalog -- HipSTR stays lightweight (~0.5 GB) because the catalog is
# filtered to loci it can handle (see BENCHMARK_MAX_* above); IlluminaEHv5 uses all 16 cores (--threads 16) on a
# prefiltered EH json; the bw2 fork reads the EH json directly at cpu=1. Runtime as CPU-hours keeps them comparable.
BENCHMARK_TOOL_CONFIG = {
    "EHv5-bw2-optimized": {"cpu": 1,  "threads": 1,  "memory": "highmem",  "job_name": "Run EHv5:optimized-streaming", "catalog": "eh"},
    "IlluminaEHv5":       {"cpu": 16, "threads": 16, "memory": "highmem",  "job_name": "Run EHv5:streaming",           "catalog": "illumina_eh"},
    "GangSTR":            {"cpu": 1,  "threads": 1,  "memory": "standard",  "job_name": "Run GangSTR",                  "catalog": "gangstr"},
    "HipSTR":             {"cpu": 1,  "threads": 1,  "memory": "standard",  "job_name": "Run HipSTR",                   "catalog": "hipstr"},
    # Long-read tools (PacBio HiFi / ONT). Each genotypes the SAME 120bp-filtered TRExplorer catalog (converted to the
    # tool's format) on the long-read BAM. cpu/threads match each tool's own default in its hail_batch_pipelines
    # factory; "memory" is only used to report the VM RAM (factories that set no tier get Hail's 3.75 GB/core
    # "standard" default). "catalog" names the converter that feeds the tool (see _convert_eh_catalog_step) -- vamos
    # reads the raw EH json in-step so it needs no convert step. TRGT v3/v5 share one factory, selected by "docker".
    "TRGTv5":             {"cpu": 16, "threads": 16, "memory": "standard",  "job_name": "Run TRGT on",                  "catalog": "trgt",   "docker": TRGT_V5_DOCKER_IMAGE},
    "TRGTv3":             {"cpu": 16, "threads": 16, "memory": "standard",  "job_name": "Run TRGT on",                  "catalog": "trgt",   "docker": TRGT_V3_DOCKER_IMAGE},
    "LongTR":             {"cpu": 1,  "threads": 1,  "memory": "standard",  "job_name": "Run LongTR on",                "catalog": "longtr"},
    "vamos":              {"cpu": 16, "threads": 16, "memory": "standard",  "job_name": "Run Vamos on",                 "catalog": "eh"},
    "ATaRVa":             {"cpu": 4,  "threads": 4,  "memory": "standard",  "job_name": "Run ATaRVa on",                "catalog": "bed"},
    "inquiSTR":           {"cpu": 16, "threads": 16, "memory": "standard",  "job_name": "Run inquiSTR on",              "catalog": "bed"},
}
# Hail Batch RAM per core (GiB) by memory tier, used to report the VM RAM in the recorded metadata.
HAIL_MEM_GB_PER_CORE = {"lowmem": 0.9, "standard": 3.75, "highmem": 6.5}

# Primary assembly contigs (no 'chr' prefix, matching the TRExplorer 'chrom' column). The benchmark catalog is
# restricted to these (like production catalogs) so a locus's +/-flank extension near a contig start can't go
# negative -- the official Illumina EH build crashes on e.g. chrM loci whose +/-1000bp extension is a negative coord.
BENCHMARK_PRIMARY_CONTIGS = [str(i) for i in range(1, 23)] + ["X", "Y"]

# HipSTR can't genotype the largest polymorphic VNTRs (its memory explodes and it segfaults on huge expansions),
# so the benchmark catalog is restricted to loci HipSTR can handle: motif <= 9bp (HipSTR's own motif-size limit)
# and reference span <= 120bp. The SAME filtered catalog is used for every tool so the comparison is on an identical
# locus set. (These caps were chosen so HipSTR finishes -- re-verify HipSTR completes if you loosen them further.)
BENCHMARK_MAX_MOTIF_SIZE_BP = 9
BENCHMARK_MAX_LOCUS_SPAN_BP = 120

def main():
    sample_table_path = "HPRC_all_aligned_short_read_and_long_read_samples.tsv"
    df = pd.read_table(sample_table_path)

    bp = pipeline("run_genotyping_tools", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    # Run on non-preemptible (non-spot) machines when NONPREEMPTIBLE=1, for long runs that must not be interrupted.
    if os.environ.get("NONPREEMPTIBLE", "").lower() in ("1", "true", "yes"):
        bp.default_preemptible(False)

    parser = bp.get_config_arg_parser()
    parser.add_argument("-s", "--sample-id", action="append",
                        help="Process only this sample. Can be specified more than once.")
    parser.add_argument("-t", "--tool", action="append", choices=SHORT_READ_TOOLS|LONG_READ_TOOLS, help="The tool to run.")
    parser.add_argument("--data-type", action="append", choices=SHORT_READ_DATA_TYPES|LONG_READ_DATA_TYPES, help="Which data type(s) to process")
    parser.add_argument("-k", "--filename-keyword", help="If specified, only BAM paths that contain this keyword will be processed", action="append")
    parser.add_argument("--filter-vcf-dir", default="gs://str-truth-set-v2/filter_vcf_v2",
                        help="Base dir under which the per-sample tool catalogs live ({sample_id}/{sample_id}.EHv5* "
                             "etc.). The build-catalogs step writes them here (from the truth set in "
                             "--truth-set-genotypes-dir) and the genotyping steps read them back from here.")
    parser.add_argument("--truth-set-genotypes-dir", default="gs://str-truth-set-v2/filter_vcf_v2",
                        help="Base dir for the filter_vcf_to_tandem_repeats genotype step output "
                             "({sample_id}/{sample_id}.tandem_repeat_genotypes.tsv.gz), used as the truth set "
                             "(it carries the per-allele repeat purity used by the purity-stratified plots)")
    parser.add_argument("--custom-catalog-path", help="If specified, use this catalog instead of the filter_vcf catalogs")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--output-subdir", help="If specified, append this extra subdirectory after the "
                        "{sample}/{data_type}/{tool}/{coverage}_coverage/ output path. Useful to keep a "
                        "--custom-catalog-path run's results separate from the per-sample-catalog results.")

    # --benchmark-resources mode (handled by run_resource_benchmark, which reuses --output-dir but ignores the
    # --tool/--data-type/--sample-id/--custom-catalog-path selection args above).
    parser.add_argument("--benchmark-resources", action="store_true",
                        help="Resource-benchmark mode: subsample the TRExplorer catalog to the N most-polymorphic "
                             "loci for several catalog sizes, run one tool at cpu=1 on each, and record wall-clock "
                             "runtime, peak RSS, and Hail Batch cost into resource_metrics.json (read by "
                             "docs/tool_comparison_viewer.html). Reuses --output-dir but ignores the "
                             "--tool/--data-type/--sample-id selection args. Run with the tool's combine-skip flag "
                             "(--skip-combine-expansion-hunter-step for the EH tools, --skip-combine-gangstr-step for "
                             "GangSTR, --skip-combine-hipstr-step for HipSTR) so the benchmark doesn't pay for the "
                             "combine step (recorded metrics are genotyping-only).")
    parser.add_argument("--benchmark-tool", default=BENCHMARK_TOOL, choices=SHORT_READ_TOOLS|LONG_READ_TOOLS,
                        help="Tool to benchmark in --benchmark-resources mode.")
    parser.add_argument("--benchmark-sample-id", default="HG002", help="Sample to benchmark.")
    parser.add_argument("--benchmark-data-type", default="illumina", help="Sequencing data type to benchmark.")
    parser.add_argument("--benchmark-coverage-keyword", default="downsampled_to_10x",
                        help="Substring selecting which coverage's read-data row to benchmark (matched against "
                             "read_data_path).")
    parser.add_argument("--benchmark-catalog-sizes", default=",".join(str(s) for s in BENCHMARK_CATALOG_SIZES),
                        help="Comma-separated catalog sizes (# loci) to benchmark, smallest first.")
    parser.add_argument("--trexplorer-bigquery-table", default=TREXPLORER_BIGQUERY_TABLE,
                        help="BigQuery table (project.dataset.table) of the TRExplorer catalog to subsample.")
    parser.add_argument("--benchmark-scrape-batch-id", type=int,
                        help="Recovery path for --benchmark-resources: skip submitting any steps and just scrape "
                             "this already-finished Hail Batch id (printed by an earlier run) and (re)write "
                             "resource_metrics.json. Use the same --benchmark-* args as the original run.")
    args = bp.parse_known_args()

    if args.benchmark_resources:
        run_resource_benchmark(bp, args, df)
        return

    if not args.tool:
        args.tool = ["TRGTv5"]
    if not args.data_type:
        args.data_type = ["pacbio"]

    if args.sample_id:
        df = df[df.sample_id.isin(args.sample_id)]

    if args.data_type:
        df = df[df.sequencing_data_type.isin(args.data_type)]

    if not args.sample_id:
        # default to HG002 + CHM1_CHM13 (the standard tool-evaluation samples) unless explicit --sample-id given
        df = df[df.sample_id.isin(["HG002", "CHM1_CHM13"])]

    # A custom catalog's results must not collide with the per-sample-catalog results under the default output dir;
    # require either a non-default --output-dir or an --output-subdir that nests them under a separate path.
    if args.custom_catalog_path and args.output_dir == DEFAULT_OUTPUT_DIR and not args.output_subdir:
        parser.error("--custom-catalog-path is set without also setting --output-dir or --output-subdir")

    # Precache only the output subtrees actually selected, not the whole bucket. A single
    # precache of "{output_dir}/**/*.*" forces gcloud to list ALL ~90k objects (the ~1200 svgs per combo
    # dominate) and filter client-side -- 5-9 min, the slowest part of every run. Instead precache one narrow
    # "{sample}/{data_type}/{tool}/**/*.tsv.gz" prefix per selected (sample, data_type, tool): each lists only
    # that tool's handful of coverage dirs (seconds). *.tsv.gz covers the add-columns for_comparison table and the
    # combine-step .tsv.gz outputs, so skip detection works for those steps. It does NOT cover ExpansionHunter's
    # per-shard genotyping outputs, which are written as json/*.json (create_expansion_hunter_steps step1), so a
    # plain re-run RE-GENOTYPES EHv5/EHv5-bw2-optimized/IlluminaEHv5 instead of skipping. To replot those without
    # re-genotyping, pass --skip-run-expansion-hunter-step --skip-combine-expansion-hunter-step --skip-add-columns-step
    # together with --force-plot-accuracy-step.
    precache_tools = args.tool
    for precache_sample in df.sample_id.unique():
        for precache_data_type in df.loc[df.sample_id == precache_sample, "sequencing_data_type"].unique():
            for precache_tool in precache_tools:
                bp.precache_file_paths(os.path.join(
                    args.output_dir, precache_sample, precache_data_type, precache_tool, "**/*.tsv.gz"))
                # also precache the per-shard genotyping json so existing shards are skipped without re-genotyping
                # (lists only the json dirs -- ~20 files per combo, fast). This makes a no-force re-run surgical:
                # only shards whose json is absent are re-run, the rest are reused. The *.json* glob matches both
                # the uncompressed .json and the .json.gz the bw2 fork writes (-z is always on for EHv5/EHv5-bw2-optimized).
                bp.precache_file_paths(os.path.join(
                    args.output_dir, precache_sample, precache_data_type, precache_tool, "**/json/*.json*"))
    # precache the prefiltered IlluminaEHv5 catalog(s) so the prefilter step is skipped once it already exists
    bp.precache_file_paths(os.path.join(args.filter_vcf_dir, "**/*.for_illumina_eh.json"))


    download_to_dir = "../results"

    # IlluminaEHv5 prefilter steps, keyed by source EHv5 catalog path so the catalog is filtered once and the step
    # is reused across all IlluminaEHv5 data types/coverages; values are (step, filtered_catalog_path) tuples.
    illumina_eh_prefilter_steps = {}

    # Per-sample build-catalogs steps, keyed by sample_id so each sample's catalogs are built once and reused
    # across all its data types/coverages/tools; values are (step, catalog_paths_by_tool) tuples.
    variant_catalog_steps = {}

    for row_i, (_, row) in enumerate(df.iterrows()):
        if args.filename_keyword:
            if not any(keyword in row.read_data_path for keyword in args.filename_keyword):
                continue

        coverage = int(round(float(row.depth_of_coverage)))
        # RNA-seq rows carry total bases sequenced (Gbp) rather than genome coverage, so label them
        # "{N}G" (e.g. "24G") instead of "{N}x" in the output dir, plots, and the Coverage column.
        coverage_label = f"{coverage}G" if row.sequencing_data_type in RNASEQ_DATA_TYPES else f"{coverage}x"
        for tool in args.tool:
            if tool in SHORT_READ_TOOLS and row.sequencing_data_type not in SHORT_READ_DATA_TYPES:
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                      f"doesn't support {row.sequencing_data_type} data")
                continue
            if tool in LONG_READ_TOOLS and row.sequencing_data_type not in LONG_READ_DATA_TYPES:
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} since {tool} "
                      f"doesn't support {row.sequencing_data_type} data")
                continue

            # GangSTR and IlluminaEHv5 are only run on illumina and illumina_exome data (GangSTR never completes on
            # ultima, and the original Illumina ExpansionHunter build is only meaningful on Illumina WGS/exome data).
            if tool in ("GangSTR", "IlluminaEHv5") and row.sequencing_data_type not in ("illumina", "illumina_exome"):
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} "
                      f"({tool} is only run on illumina and illumina_exome data)")
                continue

            # The EH+HipSTR+GangSTR EnsembleTR mode needs GangSTR, which is only run on illumina/illumina_exome.
            if tool == "EnsembleTR-EH+HipSTR+GangSTR" and row.sequencing_data_type not in ("illumina", "illumina_exome"):
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} "
                      f"(the EH+HipSTR+GangSTR mode is only run on illumina and illumina_exome data)")
                continue

            # vamos is not run on pacbio_isoseq data.
            if tool == "vamos" and row.sequencing_data_type == "pacbio_isoseq":
                print(f"WARNING: Skipping {tool} for {row.sample_id} {row.sequencing_data_type} "
                      f"(vamos is not run on pacbio_isoseq data)")
                continue

            # Build this sample's tool catalogs once (reused across its data types/coverages/tools) from the truth
            # set, unless a --custom-catalog-path was given. Cached by sample_id like illumina_eh_prefilter_steps.
            build_catalogs_step = None
            catalog_paths_by_tool = None
            if not args.custom_catalog_path:
                if row.sample_id not in variant_catalog_steps:
                    variant_catalog_steps[row.sample_id] = create_variant_catalogs_step(
                        bp,
                        sample_id=row.sample_id,
                        genotypes_tsv_path=os.path.join(args.truth_set_genotypes_dir, row.sample_id,
                            f"{row.sample_id}.tandem_repeat_genotypes.tsv.gz"),
                        output_dir=os.path.join(args.filter_vcf_dir, row.sample_id))
                build_catalogs_step, catalog_paths_by_tool = variant_catalog_steps[row.sample_id]

            # Resolve the catalog path(s) this tool reads. With --custom-catalog-path, hfs.ls it as before.
            # Otherwise reference the build-catalogs step's known single-shard outputs (created above); hfs.ls-ing
            # them here would find nothing until that step runs, so the genotyping steps depend on it instead.
            if args.custom_catalog_path:
                repeat_catalog_paths = [x.path for x in hfs.ls(args.custom_catalog_path)]
            elif tool in ("inquiSTR", "ATaRVa"):
                # inquiSTR and ATaRVa both genotype from the plain {sample_id}.bed.gz loci catalog
                # (chrom, start0, end, motif, motif_length; bgzipped + tabix-indexed)
                repeat_catalog_paths = [catalog_paths_by_tool["inquiSTR"]]
            elif tool == "vamos" or tool in ("EHv5", "EHv5-bw2-optimized", "IlluminaEHv5") or tool in ENSEMBLETR_TOOLS:
                # vamos, all three ExpansionHunter v5 variants, and both EnsembleTR modes read the single unsharded
                # EHv5 catalog json (IlluminaEHv5 prefilters it first; vamos/EnsembleTR derive from it in-step)
                repeat_catalog_paths = [catalog_paths_by_tool["EHv5"]]
            elif tool in ("TRGTv3", "TRGTv5"):
                # both TRGT versions read the same TRGT BED catalog
                repeat_catalog_paths = [catalog_paths_by_tool["TRGT"]]
            elif tool in ("GangSTR", "HipSTR", "LongTR"):
                repeat_catalog_paths = [catalog_paths_by_tool[tool]]

            print(f"Catalogs for {tool}: {repeat_catalog_paths}")
            output_dir = os.path.join(args.output_dir, row.sample_id, row.sequencing_data_type, tool, f"{coverage_label}_coverage")
            if args.output_subdir:
                output_dir = os.path.join(output_dir, args.output_subdir)
            if tool in ("EHv5", "EHv5-bw2-optimized", "IlluminaEHv5"):
                # Three ExpansionHunter v5 variants, all genotyped with create_expansion_hunter_steps:
                #   EHv5               - bw2 fork, --analysis-mode low-mem-streaming
                #   EHv5-bw2-optimized - bw2 fork, --analysis-mode optimized-streaming --improved-genotyping
                #   IlluminaEHv5       - original Illumina build, --analysis-mode streaming
                use_illumina_expansion_hunter = (tool == "IlluminaEHv5")
                if tool == "EHv5":
                    analysis_mode = "low-mem-streaming"
                elif tool == "EHv5-bw2-optimized":
                    # optimized-streaming implies --improved-genotyping inside create_expansion_hunter_steps.
                    # EHV5_ANALYSIS_MODE overrides this (e.g. "streaming") to benchmark other bw2-fork modes
                    # with the same binary; cpu/threads/memory then come from EHV5_STREAMING_* as usual.
                    analysis_mode = os.environ.get("EHV5_ANALYSIS_MODE", "optimized-streaming")
                else:  # IlluminaEHv5
                    analysis_mode = "streaming"

                catalog_prefilter_step = None
                if args.custom_catalog_path:
                    variant_catalog_file_paths = repeat_catalog_paths
                elif tool == "IlluminaEHv5":
                    # The official Illumina ExpansionHunter build aborts on the first catalog locus with >5 Ns in its
                    # +/-1000bp flanks ("Flanks can contain at most 5 characters N but found x Ns"); the bw2 fork
                    # tolerates them. So prefilter the single unsharded EHv5 catalog to drop those loci before
                    # genotyping. The prefilter step is built once per source catalog (cached in
                    # illumina_eh_prefilter_steps) and reused across all IlluminaEHv5 data types/coverages.
                    source_eh_catalog_path = next(p for p in repeat_catalog_paths if p.endswith(".EHv5.001_of_001.json"))
                    if source_eh_catalog_path not in illumina_eh_prefilter_steps:
                        prefilter_step = create_illumina_eh_catalog_prefilter_step(
                            bp,
                            eh_catalog_path=source_eh_catalog_path,
                            reference_fasta=REFERENCE_FASTA_PATH,
                            reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                            output_dir=os.path.join(args.filter_vcf_dir, row.sample_id))
                        # the prefilter reads the build-catalogs step's EHv5 json, so it must wait for that step
                        if build_catalogs_step is not None:
                            prefilter_step[0].depends_on(build_catalogs_step)
                        illumina_eh_prefilter_steps[source_eh_catalog_path] = prefilter_step
                    catalog_prefilter_step, filtered_catalog_path = illumina_eh_prefilter_steps[source_eh_catalog_path]
                    variant_catalog_file_paths = [filtered_catalog_path]
                else:
                    # the streaming EHv5 / EHv5-bw2-optimized variants use the single unsharded catalog; match the
                    # exact .EHv5.001_of_001.json so the IlluminaEHv5 prefilter outputs that share the
                    # EHv5.001_of_001 stem (.for_illumina_eh.json, .without_loci_with_flanking_Ns.json/.filtered_loci.txt)
                    # are excluded
                    variant_catalog_file_paths = [p for p in repeat_catalog_paths if p.endswith(".EHv5.001_of_001.json")]
                    # genotyping reads the build-catalogs step's EHv5 json, so gate it on that step
                    catalog_prefilter_step = build_catalogs_step

                current_step = create_expansion_hunter_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    variant_catalog_file_paths=variant_catalog_file_paths,
                    output_dir=output_dir,
                    output_prefix= f"{row.sample_id}.{tool}",
                    analysis_mode=analysis_mode,
                    loci_to_exclude=None,
                    min_locus_coverage=None,
                    use_illumina_expansion_hunter=use_illumina_expansion_hunter,
                    # IlluminaEHv5 genotypes the single prefiltered catalog and must wait for the prefilter step
                    catalog_prefilter_step=catalog_prefilter_step,
                    # EHv5/EHv5-bw2-optimized stream single-threaded. Always run UNSHARDED (one job over the whole
                    # catalog) -- both for the in-run per-sample build (whose catalog doesn't exist at construction,
                    # so it can't be sharded) and for a --custom-catalog-path (sharding it would re-download +
                    # re-scan the full CRAM once per shard). The single job is sized cpu=2 / --threads 4 / highmem,
                    # the June cost benchmark's balanced optimum: the 4 threads parallelize the htslib CRAM
                    # decompression scan (~half the cost of cpu=4/threads=8 for ~25-35% more wall). Baked in code,
                    # not env vars. No-op for IlluminaEHv5 (hardcoded 16/highmem).
                    num_shards=1,
                    streaming_cpu=2,
                    streaming_threads=4,
                    streaming_memory="highmem")
            elif tool == "GangSTR":
                current_step = create_gangstr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    repeat_spec_file_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}",
                    # wait for the build-catalogs step that produces this sample's GangSTR bed (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
            elif tool == "HipSTR":
                current_step = create_hipstr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    regions_bed_file_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}",
                    # wait for the build-catalogs step that produces this sample's HipSTR bed (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
            elif tool in ("TRGTv3", "TRGTv5"):
                current_step = create_trgt_step(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    trgt_catalog_bed_paths=repeat_catalog_paths,
                    parse_reference_region_from_locus_id=True,
                    output_dir=output_dir,
                    output_prefix= f"{row.sample_id}.{tool}",
                    docker_image=TRGT_V3_DOCKER_IMAGE if tool == "TRGTv3" else TRGT_V5_DOCKER_IMAGE,
                    # wait for the build-catalogs step that produces this sample's TRGT bed (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
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
                    output_prefix=f"{row.sample_id}.{tool}",
                    # ONT base qualities sit below LongTR's default --min-mean-qual 30, so without a lower
                    # threshold every ONT read is filtered out and LongTR emits an empty VCF.
                    min_mean_qual=10 if row.sequencing_data_type == "ONT" else None,
                    # wait for the build-catalogs step that produces this sample's LongTR bed (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
            elif tool == "inquiSTR":
                current_step = create_inquistr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    inquistr_catalog_bed_paths=repeat_catalog_paths,
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}",
                    # wait for the build-catalogs step that produces this sample's loci bed (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
            elif tool == "ATaRVa":
                # ATaRVa reads the same bgzipped + tabix-indexed {sample}.bed.gz loci catalog as inquiSTR
                current_step = create_atarva_step(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    male_or_female=row.male_or_female,
                    regions_bed_path=repeat_catalog_paths[0],
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}",
                    # wait for the build-catalogs step that produces this sample's loci bed (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
            elif tool == "vamos":
                # use the unsharded ExpansionHunter catalog json; the vamos step converts it to a vamos catalog
                current_step = create_vamos_step(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    input_bam=row.read_data_path,
                    input_bai=row.read_data_index_path,
                    # vamos needs the single unsharded EHv5 catalog; match the exact .EHv5.001_of_001.json so the
                    # IlluminaEHv5 prefilter outputs sharing the EHv5.001_of_001 stem (.for_illumina_eh.json,
                    # .without_loci_with_flanking_Ns.json/.filtered_loci.txt) don't make this a multi-catalog list and
                    # crash the vamos step. A --custom-catalog-path is passed through unfiltered (its filename won't
                    # contain that token, so filtering would leave an empty list and crash the vamos step).
                    expansion_hunter_catalog_paths=(
                        repeat_catalog_paths if args.custom_catalog_path
                        else [p for p in repeat_catalog_paths if p.endswith(".EHv5.001_of_001.json")]),
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}",
                    # wait for the build-catalogs step that produces this sample's EHv5 json (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
            elif tool in ENSEMBLETR_TOOLS:
                # EnsembleTR merges the already-computed per-caller outputs for this (sample, data_type, coverage):
                # the ExpansionHunter json (from ENSEMBLETR_EH_SOURCE_TOOL) and the HipSTR (+GangSTR) native VCFs.
                tool_results_base = os.path.join(args.output_dir, row.sample_id, row.sequencing_data_type)
                eh_json_glob = os.path.join(
                    tool_results_base, ENSEMBLETR_EH_SOURCE_TOOL, f"{coverage_label}_coverage", "json", "*.json*")
                hipstr_vcf_glob = os.path.join(
                    tool_results_base, "HipSTR", f"{coverage_label}_coverage", "vcf", "*.vcf.gz")
                eh_json_paths = sorted(x.path for x in hfs.ls(eh_json_glob))
                hipstr_vcf_paths = sorted(x.path for x in hfs.ls(hipstr_vcf_glob))
                if not eh_json_paths:
                    raise ValueError(f"No ExpansionHunter json files found for EnsembleTR at {eh_json_glob}")
                if not hipstr_vcf_paths:
                    raise ValueError(f"No HipSTR vcf files found for EnsembleTR at {hipstr_vcf_glob}")
                gangstr_vcf_paths = None
                if tool == "EnsembleTR-EH+HipSTR+GangSTR":
                    gangstr_vcf_glob = os.path.join(
                        tool_results_base, "GangSTR", f"{coverage_label}_coverage", "vcf", "*.vcf.gz")
                    gangstr_vcf_paths = sorted(x.path for x in hfs.ls(gangstr_vcf_glob))
                    if not gangstr_vcf_paths:
                        raise ValueError(f"No GangSTR vcf files found for EnsembleTR at {gangstr_vcf_glob}")
                current_step = create_ensembletr_steps(
                    bp,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                    eh_json_paths=eh_json_paths,
                    hipstr_vcf_paths=hipstr_vcf_paths,
                    gangstr_vcf_paths=gangstr_vcf_paths,
                    variant_catalog_path=(
                        repeat_catalog_paths[0] if args.custom_catalog_path
                        else next(p for p in repeat_catalog_paths if p.endswith(".EHv5.001_of_001.json"))),
                    output_dir=output_dir,
                    output_prefix=f"{row.sample_id}.{tool}",
                    sample_id=row.sample_id,
                    male_or_female=row.male_or_female,
                    # wait for the build-catalogs step that produces this sample's EHv5 json (None for --custom-catalog-path)
                    catalog_step=build_catalogs_step)
            else:
                raise ValueError(f"Unknown tool: {tool}")


            # TRGT, ATaRVa and HipSTR report an allele sequence in their VCF, so extract those sequences from the
            # VCF the genotyping step already wrote. No re-genotyping is involved -- this just re-reads that VCF.
            allele_sequences_step = allele_sequences_path = None
            if tool in SEQUENCE_ACCURACY_TOOLS:
                if tool == "HipSTR":
                    # HipSTR writes one vcf per catalog shard under {output_dir}/vcf/, named after the shard's bed
                    # file. Derive the expected path(s) from this run's own catalog shards -- globbing
                    # {output_dir}/vcf/ would also pick up stale shard vcfs left over from an older, differently-
                    # sharded catalog (older runs used a 6-way split) and silently mix them with this run's output.
                    tool_vcf_paths = [
                        os.path.join(output_dir, "vcf",
                                     re.sub(r"\.bed(\.gz)?$", "", os.path.basename(catalog_path)) + ".vcf.gz")
                        for catalog_path in repeat_catalog_paths]
                else:
                    tool_vcf_paths = [os.path.join(output_dir, f"{row.sample_id}.{tool}.vcf.gz")]

                allele_sequences_step, allele_sequences_path = create_extract_allele_sequences_step(
                    bp,
                    current_step,
                    tool=tool,
                    vcf_paths=tool_vcf_paths,
                    output_prefix=f"{row.sample_id}.{tool}",
                    output_dir=output_dir,
                    reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH)

            # EHv5, EHv5-bw2-optimized, and IlluminaEHv5 each keep their own label downstream (all three are
            # registered in add_tool_results_columns.py / add_concordance_columns.py / plot_tool_accuracy_by_allele_size.py).
            add_columns_step = add_tool_comparison_columns_step(
                bp,
                current_step,
                tool=tool,
                coverage_label=coverage_label,
                sample_id=row.sample_id,
                output_dir=output_dir,
                truth_set_genotypes_path=os.path.join(
                    args.truth_set_genotypes_dir, row.sample_id, f"{row.sample_id}.tandem_repeat_genotypes.tsv.gz"),
                tool2="Truth",
                allele_sequences_step=allele_sequences_step,
                allele_sequences_path=allele_sequences_path,
                download_to_dir=download_to_dir)

            plot_tool_accuracy_step = create_plot_tool_accuracy_steps(
                bp,
                add_columns_step,
                tool=tool,
                coverage_label=coverage_label,
                sequencing_data_type=row.sequencing_data_type,
                sample_id=row.sample_id,
                output_dir=output_dir)
    bp.run()


def run_resource_benchmark(bp, args, df):
    """Benchmark one tool's runtime / memory / cost vs. catalog size, for the resource viewer.

    Subsamples the TRExplorer catalog (--trexplorer-bigquery-table) to the N most-polymorphic loci (ranked by the
    per-locus HPRC256_Stdev annotation, descending) for each requested catalog size, builds the tool's catalog
    (EH json for the EH variants, prefiltered EH json for IlluminaEHv5, converted bed for GangSTR/HipSTR), and runs
    --benchmark-tool once per catalog on the selected sample/coverage. After the batch finishes, it scrapes each
    genotyping job's wall-clock runtime and peak RSS from the /usr/bin/time output in the job log and its cost from
    the Hail Batch client, and merges resource_metrics.json (the schema read by docs/tool_comparison_viewer.html)
    into the tool_results coverage dir.

    Runtime is recorded as CPU-hours (wall-clock x the VM's cpu count) so tools on different-sized VMs are
    comparable: GangSTR/HipSTR/EHv5-bw2-optimized run at cpu=1, IlluminaEHv5 at cpu=16 (see BENCHMARK_TOOL_CONFIG).
    Each tool's VM (cpu, RAM, preemptible) and full genotyping command line (filenames only) are recorded alongside
    the metrics. The metrics are for the genotyping step ONLY; run with --skip-combine-{expansion-hunter,gangstr,
    hipstr}-step so the benchmark doesn't pay for the downstream combine steps.

    The sizes are nested subsets (each is the top-N most-polymorphic loci); note the actual genotyped locus count is
    smaller than N for tools that filter (HipSTR drops motifs >9bp; IlluminaEHv5 drops flanking-N and >=500bp loci).

    Args:
        bp: the step_pipeline pipeline object.
        args: parsed args (uses the --benchmark-* and --trexplorer-bigquery-table options).
        df: the sample table (HPRC_all_aligned_short_read_and_long_read_samples.tsv), used to resolve the sample's
            read-data path, index, sex, and coverage.
    """
    tool = args.benchmark_tool
    if tool not in BENCHMARK_TOOL_CONFIG:
        raise ValueError(f"--benchmark-resources supports {sorted(BENCHMARK_TOOL_CONFIG)}, got '{tool}'")

    sizes = sorted(int(s) for s in args.benchmark_catalog_sizes.split(","))

    # Resolve the one read-data row for the requested sample / data type / coverage.
    rows = df[(df.sample_id == args.benchmark_sample_id)
              & (df.sequencing_data_type == args.benchmark_data_type)
              & (df.read_data_path.str.contains(args.benchmark_coverage_keyword, regex=False))]
    if len(rows) != 1:
        raise ValueError(f"Expected exactly 1 read-data row for {args.benchmark_sample_id} "
                         f"{args.benchmark_data_type} matching '{args.benchmark_coverage_keyword}', found {len(rows)}")
    row = rows.iloc[0]
    cov = int(round(float(row.depth_of_coverage)))
    # RNA-seq rows are labeled by total bases sequenced (Gbp) not fold-coverage, matching main() and the viewer's
    # SAMPLES tok (e.g. "24G"); using "x" here would write to a coverage dir the viewer never requests.
    coverage_label = f"{cov}G" if args.benchmark_data_type in RNASEQ_DATA_TYPES else f"{cov}x"

    # resource_metrics.json goes into the coverage dir the viewer reads
    # ({output_dir}/{sample}/{data_type}/{coverage}_coverage/resource_metrics.json); the subsampled catalogs and
    # per-size EH outputs go under a resource_benchmark subtree so they don't collide with the accuracy results.
    coverage_dir = os.path.join(args.output_dir, args.benchmark_sample_id, args.benchmark_data_type,
                                f"{coverage_label}_coverage")
    benchmark_base = os.path.join(args.output_dir, args.benchmark_sample_id, args.benchmark_data_type,
                                  "resource_benchmark", tool, f"{coverage_label}_coverage")
    billing_project = getattr(args, "batch_billing_project", None) or "tgg-rare-disease"

    # Recovery path: with --benchmark-scrape-batch-id, skip resubmitting any steps and just (re)scrape a finished
    # batch and write resource_metrics.json. Used when a blocking run's scrape/write failed after the (paid) batch
    # already ran, or to finish a --no-wait submission. The other --benchmark-* args must match the original run.
    if args.benchmark_scrape_batch_id:
        _scrape_and_write(args.benchmark_scrape_batch_id, billing_project, sizes, tool, coverage_dir)
        return

    # 1. Query TRExplorer for the top max(sizes) loci by HPRC256_Stdev, build nested EH catalogs, upload to GCS.
    #    (GangSTR/HipSTR convert these to bed; IlluminaEHv5 prefilters them; the bw2 EH fork reads them directly.)
    #    The 120bp-filtered nested catalogs depend only on (bigquery_table, sizes, filter params) -- NOT on the
    #    sample / coverage / tool -- so every combo shares ONE catalog dir (keyed by the filter params so a filter
    #    change can't silently reuse stale catalogs), and _build_trexplorer_catalogs skips any size already uploaded.
    shared_catalogs_dir = os.path.join(args.output_dir, "resource_benchmark",
                                       f"shared_catalogs_span{BENCHMARK_MAX_LOCUS_SPAN_BP}_motif{BENCHMARK_MAX_MOTIF_SIZE_BP}")
    eh_catalog_paths = _build_trexplorer_catalogs(args.trexplorer_bigquery_table, sizes, shared_catalogs_dir)

    # 2. Per catalog size: build the tool's catalog (convert/prefilter as needed) + one genotyping step on the tool's
    #    VM. The per-size catalog basename carries the trexplorer_top_{size} token, so the genotyping job name does too,
    #    which is how _scrape_benchmark_metrics keys the catalog size.
    for size in sizes:
        eh_catalog = eh_catalog_paths[size]
        size_dir = os.path.join(benchmark_base, f"size_{size}")
        output_prefix = f"{args.benchmark_sample_id}.{tool}.trexplorer_top_{size}"

        if tool in ("EHv5-bw2-optimized", "IlluminaEHv5"):
            use_illumina = (tool == "IlluminaEHv5")
            if use_illumina:
                # IlluminaEHv5 needs the flanking-N / >=500bp-locus prefilter; cpu/memory are forced to 16/highmem
                # inside the factory (the streaming_* args below are ignored for it).
                prefilter_step, catalog_path = create_illumina_eh_catalog_prefilter_step(
                    bp, eh_catalog_path=eh_catalog, reference_fasta=REFERENCE_FASTA_PATH,
                    reference_fasta_fai=REFERENCE_FASTA_FAI_PATH, output_dir=os.path.join(size_dir, "catalog"))
                analysis_mode = "streaming"
            else:
                prefilter_step, catalog_path, analysis_mode = None, eh_catalog, "optimized-streaming"
            create_expansion_hunter_steps(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                male_or_female=row.male_or_female,
                variant_catalog_file_paths=[catalog_path],
                output_dir=size_dir,
                output_prefix=output_prefix,
                analysis_mode=analysis_mode,
                loci_to_exclude=None,
                min_locus_coverage=None,
                use_illumina_expansion_hunter=use_illumina,
                catalog_prefilter_step=prefilter_step,
                num_shards=1,
                streaming_cpu=1,
                streaming_threads=1,
                streaming_memory="highmem")
        elif tool == "GangSTR":
            convert_step, bed_path = _convert_eh_catalog_step(bp, eh_catalog, "GangSTR", size,
                                                              os.path.join(size_dir, "catalog"))
            create_gangstr_steps(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                male_or_female=row.male_or_female,
                repeat_spec_file_paths=[bed_path],
                output_dir=size_dir,
                output_prefix=output_prefix,
                catalog_step=convert_step)
        elif tool == "HipSTR":
            convert_step, bed_path = _convert_eh_catalog_step(bp, eh_catalog, "HipSTR", size,
                                                              os.path.join(size_dir, "catalog"))
            create_hipstr_steps(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                male_or_female=row.male_or_female,
                regions_bed_file_paths=[bed_path],
                output_dir=size_dir,
                output_prefix=output_prefix,
                catalog_step=convert_step,
                cpu=BENCHMARK_TOOL_CONFIG["HipSTR"]["cpu"],
                memory=BENCHMARK_TOOL_CONFIG["HipSTR"]["memory"])
        elif tool in ("TRGTv5", "TRGTv3"):
            convert_step, bed_path = _convert_eh_catalog_step(
                bp, eh_catalog, "TRGT", size, os.path.join(size_dir, "catalog"),
                reference_fasta=REFERENCE_FASTA_PATH, reference_fasta_fai=REFERENCE_FASTA_FAI_PATH)
            create_trgt_step(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                male_or_female=row.male_or_female,
                trgt_catalog_bed_paths=[bed_path],
                output_dir=size_dir,
                output_prefix=output_prefix,
                # the downstream json conversion derives the ReferenceRegion from the TRGT VCF's own coords (the
                # single-motif catalog takes the --parse-genotype-from-AL-field path), so it never has to parse the
                # catalog's locus-id -- keeping the genotyping step from failing on an unexpected TRID format.
                parse_reference_region_from_locus_id=False,
                cpu=BENCHMARK_TOOL_CONFIG[tool]["cpu"],
                docker_image=BENCHMARK_TOOL_CONFIG[tool]["docker"],
                catalog_step=convert_step)
        elif tool == "LongTR":
            convert_step, bed_path = _convert_eh_catalog_step(bp, eh_catalog, "LongTR", size,
                                                              os.path.join(size_dir, "catalog"))
            create_longtr_steps(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                male_or_female=row.male_or_female,
                regions_bed_paths=[bed_path],
                output_dir=size_dir,
                output_prefix=output_prefix,
                # ONT reads are noisier; match the main pipeline's ONT quality floor
                min_mean_qual=(10 if args.benchmark_data_type == "ONT" else None),
                catalog_step=convert_step)
        elif tool == "vamos":
            # vamos converts the EH json to its own motif catalog IN-STEP, so it reads the raw EH catalog directly
            # (no separate convert step); cpu drives both the CRAM/BAM scan threads and vamos -t.
            create_vamos_step(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                expansion_hunter_catalog_paths=[eh_catalog],
                output_dir=size_dir,
                output_prefix=output_prefix,
                cpu=BENCHMARK_TOOL_CONFIG["vamos"]["cpu"])
        elif tool == "ATaRVa":
            convert_step, bed_path = _convert_eh_catalog_step(bp, eh_catalog, "ATaRVa", size,
                                                              os.path.join(size_dir, "catalog"))
            create_atarva_step(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                male_or_female=row.male_or_female,
                regions_bed_path=bed_path,   # ATaRVa takes a single bed path (not a list)
                output_dir=size_dir,
                output_prefix=output_prefix,
                cpu=BENCHMARK_TOOL_CONFIG["ATaRVa"]["cpu"],
                catalog_step=convert_step)
        elif tool == "inquiSTR":
            convert_step, bed_path = _convert_eh_catalog_step(bp, eh_catalog, "inquiSTR", size,
                                                              os.path.join(size_dir, "catalog"))
            create_inquistr_steps(
                bp,
                reference_fasta=REFERENCE_FASTA_PATH,
                reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
                input_bam=row.read_data_path,
                input_bai=row.read_data_index_path,
                male_or_female=row.male_or_female,
                inquistr_catalog_bed_paths=[bed_path],
                output_dir=size_dir,
                output_prefix=output_prefix,
                cpu=BENCHMARK_TOOL_CONFIG["inquiSTR"]["cpu"],
                catalog_step=convert_step)

    # 3. Submit + wait. bp.run() blocks unless --no-wait, and returns the hail batch handle (result.id).
    result = bp.run()
    if result is None:
        print("No steps were run (all skipped?); not writing resource_metrics.json.")
        return
    batch_id = getattr(result, "id", None)
    print(f"Benchmark batch id: {batch_id}")
    if getattr(args, "no_wait", False) or getattr(args, "dry_run", False):
        print(f"--no-wait/--dry-run does not scrape metrics. Once the batch finishes, write resource_metrics.json "
              f"by re-running with --benchmark-scrape-batch-id {batch_id} (plus the same --benchmark-* args).")
        return

    # 4. Scrape wall-clock + peak RSS (from the /usr/bin/time output in each job log) and cost (from the Hail Batch
    #    client), keyed by catalog size, then write resource_metrics.json. The batch already ran (and was paid for),
    #    so on failure surface the recovery command instead of silently losing the results.
    try:
        _scrape_and_write(batch_id, billing_project, sizes, tool, coverage_dir)
    except Exception:
        print(f"ERROR scraping/writing metrics for finished batch {batch_id}. Retry without resubmitting the batch "
              f"by re-running with --benchmark-scrape-batch-id {batch_id} (plus the same --benchmark-* args).")
        raise


def _scrape_and_write(batch_id, billing_project, sizes, tool, coverage_dir):
    """Scrape a finished benchmark batch and merge its metrics + provenance + VM metadata into resource_metrics.json."""
    print(f"Scraping metrics from finished batch {batch_id} ...")
    metrics_by_size, vm = _scrape_benchmark_metrics(batch_id, billing_project, sizes, tool)
    _write_resource_metrics_json(metrics_by_size, sizes, tool, coverage_dir, vm or _benchmark_vm(tool))


def _benchmark_vm(tool):
    """Table/env fallback VM metadata for this tool: {cpu, mem_gb, preemptible} (used when no job was scraped)."""
    cfg = BENCHMARK_TOOL_CONFIG[tool]
    preemptible = os.environ.get("NONPREEMPTIBLE", "").lower() not in ("1", "true", "yes")
    return {"cpu": cfg["cpu"], "mem_gb": round(cfg["cpu"] * HAIL_MEM_GB_PER_CORE[cfg["memory"]], 2),
            "preemptible": preemptible}


def _read_job_vm(status, tool):
    """Read the VM (cpu, mem_gb, preemptible) from a finished job's Hail Batch status, with the table as fallback.

    cpu is the actual allocated cores (msec_mcpu / duration_ms / 1000, Hail's own accounting); preemptible is read
    from the cost_breakdown resource strings (which carry a preemptible/nonpreemptible token); RAM is the provisioned
    tier RAM for that cpu count (cpu x HAIL_MEM_GB_PER_CORE[tier]). Falls back to BENCHMARK_TOOL_CONFIG + env when a
    field can't be read from the job.
    """
    cfg = BENCHMARK_TOOL_CONFIG[tool]
    cpu = cfg["cpu"]
    msec_mcpu, duration = status.get("msec_mcpu"), status.get("duration")
    if msec_mcpu and duration:
        cpu = max(1, round(msec_mcpu / duration / 1000))
    resource_strs = " ".join((cb or {}).get("resource", "") for cb in (status.get("cost_breakdown") or []))
    if "nonpreemptible" in resource_strs:      # check nonpreemptible first (it contains "preemptible")
        preemptible = False
    elif "preemptible" in resource_strs:
        preemptible = True
    else:
        preemptible = os.environ.get("NONPREEMPTIBLE", "").lower() not in ("1", "true", "yes")
    return {"cpu": cpu, "mem_gb": round(cpu * HAIL_MEM_GB_PER_CORE[cfg["memory"]], 2), "preemptible": preemptible}


def _sanitize_command(cmd):
    """Collapse whitespace and replace every path-like token with just its basename (record filenames, not paths)."""
    return " ".join(t.rsplit("/", 1)[-1] if "/" in t else t for t in cmd.split())


def _convert_eh_catalog_step(bp, eh_catalog_path, tool, size, output_dir, reference_fasta=None, reference_fasta_fai=None):
    """Convert an ExpansionHunter JSON catalog to the tool's catalog format via the str_analysis converters baked into
    RUN_TOOLS_DOCKER_IMAGE. The output basename carries the trexplorer_top_{size} token so the downstream genotyping
    job name does too (that's how _scrape_benchmark_metrics keys the catalog size). Returns (step, gs:// catalog path).

    tool -> converter (all baked into RUN_TOOLS_DOCKER_IMAGE):
      GangSTR         -> convert_expansion_hunter_catalog_to_gangstr_spec   (plain .bed)
      HipSTR          -> convert_expansion_hunter_catalog_to_hipstr_format  (plain .bed)
      TRGT            -> convert_expansion_hunter_catalog_to_trgt_catalog   (plain .bed; requires -R reference_fasta)
      LongTR          -> convert_expansion_hunter_catalog_to_longtr_format  (plain .bed)
      ATaRVa/inquiSTR -> convert_expansion_hunter_catalog_to_bed            (bgzipped + tabix-indexed .bed.gz; ATaRVa
                                                                             localizes the .tbi, inquiSTR ignores it;
                                                                             ATaRVa's col5 is rewritten to the motif
                                                                             length, which it requires)
    """
    module = {"GangSTR": "convert_expansion_hunter_catalog_to_gangstr_spec",
              "HipSTR": "convert_expansion_hunter_catalog_to_hipstr_format",
              "TRGT": "convert_expansion_hunter_catalog_to_trgt_catalog",
              "LongTR": "convert_expansion_hunter_catalog_to_longtr_format",
              "ATaRVa": "convert_expansion_hunter_catalog_to_bed",
              "inquiSTR": "convert_expansion_hunter_catalog_to_bed"}[tool]
    # the plain-bed converter bgzips + tabix-indexes when its -o path ends in .bed.gz (ATaRVa needs the .tbi)
    bgzipped = tool in ("ATaRVa", "inquiSTR")
    out_name = f"trexplorer_top_{size}.{tool}.bed.gz" if bgzipped else f"trexplorer_top_{size}.{tool}.bed"
    step = bp.new_step(
        name=f"Convert catalog to {tool} format (trexplorer_top_{size})",
        arg_suffix="convert-benchmark-catalog-step",
        image=RUN_TOOLS_DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        storage="10Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)
    step.command("set -ex")
    local_json = step.input(eh_catalog_path)
    ref_arg = ""
    if tool == "TRGT":
        # the TRGT converter reads each locus's reference sequence, so it requires the reference fasta (+ .fai)
        local_fasta = step.input(reference_fasta)
        step.input(reference_fasta_fai or f"{reference_fasta}.fai")
        ref_arg = f"-R {local_fasta} "
    step.command(f"python3 -m str_analysis.{module} {ref_arg}{local_json} -o {out_name}")
    if tool == "ATaRVa":
        # ATaRVa strictly requires 5 columns: chrom, start, end, motif, motif_length -- it aborts if column 5 isn't
        # the motif length. The plain-bed converter writes "." in column 5 (which inquiSTR tolerates), so rewrite
        # column 5 to len(motif) and re-index. The benchmark catalog is single-motif, so column 4 holds one motif.
        step.command(f"zcat {out_name} | awk 'BEGIN{{OFS=\"\\t\"}} {{$5=length($4); print}}' | bgzip > {out_name}.fixed")
        step.command(f"mv {out_name}.fixed {out_name}")
        step.command(f"tabix -f {out_name}")
    step.command("ls -lhrt")
    step.output(out_name)
    if bgzipped:
        step.output(f"{out_name}.tbi")
    return step, os.path.join(output_dir, out_name)


def _build_trexplorer_catalogs(bigquery_table, sizes, catalogs_gcs_dir):
    """Query the TRExplorer BigQuery catalog for the most-polymorphic loci and write one EH catalog per size.

    Selects the top max(sizes) loci ordered by HPRC256_Stdev descending (LocusId as a deterministic tiebreak),
    restricted to the primary contigs (chr1-22,X,Y; matching production catalogs, and avoiding the chrM near-start
    loci whose flank extension goes negative and crashes IlluminaEHv5), and
    excluding loci with >5 Ns in their flanks (matching TRExplorer's own ExpansionHunter catalog export — the
    official Illumina EH build rejects these; the bw2 fork tolerates them, but excluding them keeps the benchmark
    catalog identical to what TRExplorer distributes). Each size takes the top-N of that sigma-ranked list, so the
    sizes are nested subsets (same loci); the loci WITHIN each written catalog are then ordered by canonical motif
    to match how production EH catalogs are laid out (convert_truth_set_to_variant_catalogs.py sorts by canonical
    motif to improve the optimized-streaming cache hit rate), so the benchmarked runtime is representative. Each
    catalog is written as an ExpansionHunter variant catalog JSON and uploaded to
    {catalogs_gcs_dir}/trexplorer_top_{size}.EH.json.

    Args:
        bigquery_table: "project.dataset.table" of the TRExplorer catalog.
        sizes: list of catalog sizes (# loci).
        catalogs_gcs_dir: gs:// dir to upload the catalogs into.

    Returns:
        dict mapping size -> gs:// path of that size's EH catalog.
    """
    from google.cloud import bigquery   # lazy import: only needed in --benchmark-resources mode

    # The catalogs are shared across all (sample, coverage, tool) combos, so skip any size already uploaded and only
    # query/build the missing ones. When every requested size is present, skip the BigQuery query entirely.
    gcs_paths = {size: os.path.join(catalogs_gcs_dir, f"trexplorer_top_{size}.EH.json") for size in sizes}
    missing = [size for size in sizes if not hfs.exists(gcs_paths[size])]
    if not missing:
        print(f"All {len(sizes)} catalog(s) already present in {catalogs_gcs_dir}; skipping BigQuery query + rebuild.")
        return gcs_paths
    print(f"{len(sizes) - len(missing)}/{len(sizes)} catalog(s) already present in {catalogs_gcs_dir}; "
          f"building missing size(s) {missing}.")

    max_size = max(missing)
    project = bigquery_table.split(".")[0]
    contigs_in = ", ".join(f"'{c}'" for c in BENCHMARK_PRIMARY_CONTIGS)
    sql = f"""
        SELECT ReferenceRegion, CONCAT('(', ReferenceMotif, ')*') AS LocusStructure, LocusId, CanonicalMotif
        FROM `{bigquery_table}`
        WHERE HPRC256_Stdev IS NOT NULL AND NsInFlanks <= 5 AND chrom IN ({contigs_in})
          AND MotifSize <= {BENCHMARK_MAX_MOTIF_SIZE_BP}
          AND (end_1based - start_0based) <= {BENCHMARK_MAX_LOCUS_SPAN_BP}
        ORDER BY HPRC256_Stdev DESC, LocusId ASC
        LIMIT {max_size}
    """
    print(f"Querying {bigquery_table} for the top {max_size} loci by HPRC256_Stdev ...")
    rows = list(bigquery.Client(project=project).query(sql).result())
    if len(rows) < max_size:
        print(f"WARNING: only {len(rows)} loci have a non-null HPRC256_Stdev (< requested {max_size}); larger "
              f"catalog sizes will be truncated to {len(rows)}.")
    # kept in sigma-descending order (defines which loci each size includes); CanonicalMotif is carried so each
    # written catalog can be motif-sorted below, then dropped from the EH records.
    catalog_records = [
        {"LocusId": r["LocusId"], "LocusStructure": r["LocusStructure"],
         "ReferenceRegion": r["ReferenceRegion"], "VariantType": "Repeat", "_CanonicalMotif": r["CanonicalMotif"]}
        for r in rows]

    tmp_dir = tempfile.mkdtemp()
    for size in missing:
        # top-N most-polymorphic loci (prefix of the sigma-ranked list), ordered by canonical motif for a
        # production-representative layout; drop the helper _CanonicalMotif key from the written EH records.
        subset = sorted(catalog_records[:size], key=lambda r: (r["_CanonicalMotif"] or "", r["LocusId"]))
        records = [{k: v for k, v in r.items() if k != "_CanonicalMotif"} for r in subset]
        local_path = os.path.join(tmp_dir, f"trexplorer_top_{size}.EH.json")
        with open(local_path, "w") as f:
            json.dump(records, f)
        print(f"Uploading {min(size, len(catalog_records))}-locus catalog to {gcs_paths[size]}")
        if os.system(f"gcloud storage cp '{local_path}' '{gcs_paths[size]}'") != 0:
            raise RuntimeError(f"Failed to upload catalog to {gcs_paths[size]}")
    return gcs_paths


def _scrape_benchmark_metrics(batch_id, billing_project, sizes, tool):
    """Scrape per-catalog-size runtime / peak RSS / cost + provenance from a finished Hail Batch.

    Matches this tool's genotyping jobs by BENCHMARK_TOOL_CONFIG[tool]["job_name"], parses the catalog size from the
    'trexplorer_top_{size}' token in the job name, the wall-clock runtime and peak RSS from the /usr/bin/time --verbose
    output in the job log, and the cost from the Hail Batch client. Runtime is reported as CPU-hours (wall-clock x
    the job's actual allocated cpu) so different-VM tools are comparable; memory as peak RSS in GB. For each size it
    also records the batch id, job id, and the sanitized genotyping command line (from that job's `set -x` echo) so
    every data point is traceable back to the job that produced it.

    Args:
        batch_id: the finished batch's id.
        billing_project: hail batch billing project to open the client with.
        sizes: catalog sizes expected (used to warn about any missing job).
        tool: which tool's jobs to scrape (selects the job-name matcher).

    Returns:
        (metrics_by_size, vm): metrics_by_size maps size -> {"runtime": cpu_hours, "memory": gb, "cost": usd,
        "batch_id": int, "job_id": int, "command": str} (any unmeasured field is None); vm is the {cpu, mem_gb,
        preemptible} read from the first matched job, or None if no job matched.
    """
    import hailtop.batch_client.client as hb_client   # lazy import: only needed in --benchmark-resources mode

    job_name_match = BENCHMARK_TOOL_CONFIG[tool]["job_name"]
    batch = hb_client.BatchClient(billing_project).get_batch(batch_id)
    metrics = {}
    vm = None
    for job_info in batch.jobs():
        name = job_info.get("name") or ""
        if job_name_match not in name:
            continue
        # tolerant of both the current "trexplorer_top_{size}" naming and the legacy "trex_top_{size}" (batches
        # submitted before the trex -> trexplorer rename still carry trex_top_ in their job names).
        size_match = re.search(r"trex(?:plorer)?_top_(\d+)", name)
        if not size_match:
            continue
        size = int(size_match.group(1))
        # only scrape the requested sizes -- lets a scrape target a subset and skips fetching logs for sizes we
        # intentionally exclude (e.g. GangSTR's huge multi-MB 100k/200k logs that time out and are capped out anyway).
        if size not in sizes:
            continue

        job = batch.get_job(job_info["job_id"])
        status = job.status()
        cost = status.get("cost")
        if isinstance(cost, str):
            cost = float(cost.lstrip("$")) if cost.strip() else None

        # runtime CPU-hours = wall-clock x the cores the tool actually uses (its thread count), NOT the provisioned
        # cpu -- e.g. HipSTR is single-threaded on a cpu=2 VM (2nd core only for RAM), so its idle core isn't counted.
        threads = BENCHMARK_TOOL_CONFIG[tool]["threads"]
        if vm is None:
            vm = _read_job_vm(status, tool)   # VM record still uses the provisioned cpu (read from the job)

        # a slow/unavailable job log must not abort the whole scrape (one Hail log-fetch TimeoutError previously
        # crashed the run after other jobs were already fetched) -- treat an unfetchable log as "no metrics for this
        # size" (cost is still recorded from the job status above), and move on.
        try:
            log = job.log() or {}
        except Exception as e:
            print(f"WARNING: could not fetch log for job {job_info['job_id']} (catalog size {size}): {e}")
            log = {}
        main_log = log.get("main", "") if isinstance(log, dict) else str(log)

        runtime_cpu_hours = None
        wall_match = re.search(r"Elapsed \(wall clock\) time.*?\): ([0-9:.]+)", main_log)
        if wall_match:
            # e.g. "7:53.62" (m:ss) or "1:07:53.6" (h:mm:ss): sum tokens right-to-left with 60**i weights, then
            # x threads to get CPU-hours (== wall-hours for the single-threaded tools; x16 for IlluminaEHv5).
            parts = [float(x) for x in wall_match.group(1).split(":")][::-1]
            runtime_cpu_hours = sum(value * 60**i for i, value in enumerate(parts)) / 3600.0 * threads

        memory_gb = None
        rss_match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", main_log)
        if rss_match:
            memory_gb = int(rss_match.group(1)) / 1e6

        # the tool's own command (strip the /usr/bin/time --verbose measurement wrapper), filenames only
        command = None
        cmd_match = re.search(r"/usr/bin/time --verbose (\S[^\n]*)", main_log)
        if cmd_match:
            command = _sanitize_command(cmd_match.group(1))

        if runtime_cpu_hours is None or memory_gb is None:
            print(f"WARNING: could not parse runtime/RSS for catalog size {size} from job {job_info['job_id']}")
        metrics[size] = {"runtime": runtime_cpu_hours, "memory": memory_gb, "cost": cost,
                         "batch_id": batch_id, "job_id": job_info["job_id"], "command": command}

    for size in sizes:
        if size not in metrics:
            print(f"WARNING: no {tool} genotyping job found for catalog size {size}")
            metrics[size] = {k: None for k in ("runtime", "memory", "cost", "batch_id", "job_id", "command")}
    return metrics, vm


def _write_resource_metrics_json(metrics_by_size, sizes, tool, coverage_dir, vm):
    """Merge this run's metrics + provenance + metadata into resource_metrics.json (viewer schema) at {coverage_dir}.

    Schema read by docs/tool_comparison_viewer.html:
        {"catalog_sizes": [...], "tools": {tool: {"runtime": [...], "memory": [...], "cost": [...],
                                                  "batch_id": [...], "job_id": [...], "command": [...],
                                                  "vm": {"cpu", "mem_gb", "preemptible"}}}}
    Every per-size array (runtime/memory/cost + batch_id/job_id/command) is aligned index-wise to catalog_sizes
    (null where not measured), so each data point is traceable to the batch+job that produced it; vm is per-tool.

    A single benchmark run only covers one tool and the --benchmark-catalog-sizes it was given, but the viewer's
    "all" mode overlays every tool and the full-range plan re-runs with more sizes; so this reads any existing file
    and MERGES (union of catalog_sizes, this tool's arrays + metadata re-aligned and added/updated) rather than
    overwriting, which would discard previously recorded sizes/tools.
    """
    gcs_path = os.path.join(coverage_dir, "resource_metrics.json")
    per_size = ("runtime", "memory", "cost", "batch_id", "job_id", "command")

    # Unpack any existing file into {tool: {size: {field}}} + {tool: {vm}}, then overlay this run. The per-field
    # read is defensive (missing / non-list / short arrays -> None) so a schema change can't corrupt the merge.
    by_tool, meta_by_tool = {}, {}
    existing = _read_existing_resource_metrics(gcs_path)
    if existing:
        existing_sizes = existing.get("catalog_sizes", [])
        for t, entry in existing.get("tools", {}).items():
            for i, s in enumerate(existing_sizes):
                by_tool.setdefault(t, {})[s] = {
                    f: (entry[f][i] if isinstance(entry.get(f), list) and i < len(entry[f]) else None)
                    for f in per_size}
            if entry.get("vm") is not None:
                meta_by_tool[t] = {"vm": entry["vm"]}
    # Overlay this run per-field, but never overwrite a previously-recorded value with a None: on a re-run with a
    # wider size range, step_pipeline skips the already-completed sizes, so _scrape_benchmark_metrics has no job for
    # them and reports all-None -- writing those Nones would wipe the real earlier measurements + provenance.
    for s in sizes:
        slot = by_tool.setdefault(tool, {}).setdefault(s, {f: None for f in per_size})
        for f in per_size:
            if metrics_by_size[s][f] is not None:
                slot[f] = metrics_by_size[s][f]
    if vm is not None:
        meta_by_tool.setdefault(tool, {})["vm"] = vm

    all_sizes = sorted({s for size_map in by_tool.values() for s in size_map})
    data = {"catalog_sizes": all_sizes, "tools": {}}
    for t, size_map in by_tool.items():
        entry = {f: [(size_map.get(s) or {}).get(f) for s in all_sizes] for f in per_size}
        entry.update(meta_by_tool.get(t, {}))
        data["tools"][t] = entry

    local_path = os.path.join(tempfile.mkdtemp(), "resource_metrics.json")
    with open(local_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Writing resource metrics to {gcs_path}:\n{json.dumps(data, indent=2)}")
    if os.system(f"gcloud storage cp '{local_path}' '{gcs_path}'") != 0:
        raise RuntimeError(f"Failed to upload resource_metrics.json to {gcs_path}")


def _read_existing_resource_metrics(gcs_path):
    """Return the parsed resource_metrics.json at gcs_path, or None if it doesn't exist yet."""
    if not hfs.exists(gcs_path):
        return None
    with hfs.open(gcs_path, "rb") as f:
        return json.loads(f.read())


def create_variant_catalogs_step(bp, *, sample_id, genotypes_tsv_path, output_dir):
    """Build the per-tool repeat catalogs for one sample from its filter_vcf_to_tandem_repeats truth set.

    Runs convert_truth_set_to_variant_catalogs.py (baked into RUN_TOOLS_DOCKER_IMAGE) on the sample's
    {sample_id}.tandem_repeat_genotypes.tsv.gz, writing a single unsharded catalog per tool so the genotyping
    steps can reference exact 001_of_001 paths (and depend on this step) instead of hfs.ls-ing files that don't
    exist until this step runs. --gangstr-loci-per-run is set huge so GangSTR/HipSTR are single-shard like the rest.

    Args:
        bp: the step_pipeline pipeline object.
        sample_id: sample id, used as the catalog filename prefix and to name the step.
        genotypes_tsv_path: gs:// path of {sample_id}.tandem_repeat_genotypes.tsv.gz (the truth set).
        output_dir: directory to write the catalogs into (the --filter-vcf-dir/{sample_id} dir, reused across runs).

    Returns:
        A (step, catalog_paths_by_tool) tuple. catalog_paths_by_tool maps
        "EHv5"/"GangSTR"/"HipSTR"/"LongTR"/"TRGT"/"inquiSTR" to the gs:// path of that tool's catalog; the
        genotyping steps read those paths and must depend on the returned step.
    """
    catalog_filenames = {
        "EHv5": f"{sample_id}.EHv5.001_of_001.json",
        "GangSTR": f"{sample_id}.GangSTR.001_of_001.bed",
        "HipSTR": f"{sample_id}.HipSTR.001_of_001.bed",
        "LongTR": f"{sample_id}.LongTR.001_of_001.bed",
        "TRGT": f"{sample_id}.TRGT_repeat_catalog.bed",
        "inquiSTR": f"{sample_id}.bed.gz",
    }

    step = bp.new_step(
        name=f"Build tool catalogs for {sample_id}",
        arg_suffix="build-variant-catalogs-step",
        image=RUN_TOOLS_DOCKER_IMAGE,
        cpu=2,
        memory="highmem",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)
    step.command("set -ex")
    local_tsv = step.input(genotypes_tsv_path)
    # single unsharded catalog per tool (--gangstr-loci-per-run huge -> GangSTR/HipSTR are 001_of_001 too)
    step.command(
        f"python3 /run_tools/convert_truth_set_to_variant_catalogs.py "
        f"--output-filename-prefix {sample_id} "
        f"--gangstr-loci-per-run 1000000000 "
        f"--output-dir . "
        f"{local_tsv}")
    step.command("ls -lhrt")
    for fname in catalog_filenames.values():
        step.output(fname)
    # inquiSTR bed is bgzipped + tabixed by the converter, so also delocalize its index
    step.output(f"{sample_id}.bed.gz.tbi")

    catalog_paths_by_tool = {tool: os.path.join(output_dir, fname) for tool, fname in catalog_filenames.items()}
    return step, catalog_paths_by_tool


def create_illumina_eh_catalog_prefilter_step(bp, *, eh_catalog_path, reference_fasta, reference_fasta_fai, output_dir):
    """Build a step that prefilters an ExpansionHunter catalog for the official Illumina EH v5 build (IlluminaEHv5).

    The official build chokes on two classes of loci that the bw2 fork tolerates, so two filters are applied in
    sequence to produce the IlluminaEHv5 catalog:
      1. Loci with >5 Ns in their +/-1000bp flanks, which make the official build abort with
         "Flanks can contain at most 5 characters N but found x Ns" (str_analysis.filter_out_loci_with_Ns_in_flanks,
         using the EH default --region-extension-length of 1000bp).
      2. Loci whose reference interval is >= ILLUMINA_EH_MAX_LOCUS_SIZE_BP wide, which make the official build crash
         with "numIndels out of range" on large VNTRs (e.g. a 12kb/4kb-motif locus). For multi-region loci the
         interval is measured from the first region's start to the last region's end.

    The official build can't read a gzipped catalog, so the output is written as plain .json.

    Args:
        bp: the step_pipeline pipeline object.
        eh_catalog_path: gs:// path of the source ExpansionHunter (EHv5) variant catalog json.
        reference_fasta: gs:// path of the reference fasta (used to read each locus's flanking sequence).
        reference_fasta_fai: gs:// path of the reference fasta .fai index.
        output_dir: directory to write the filtered catalog into (the filter_vcf dir, so it is reused across runs).

    Returns:
        A (step, filtered_catalog_path) tuple. filtered_catalog_path is the plain-json catalog the IlluminaEHv5
        genotyping step should read; the genotyping step must depend on the returned step.
    """
    base = re.sub(r"\.json(\.gz)?$", "", os.path.basename(eh_catalog_path))
    n_filtered_filename = f"{base}.without_loci_with_flanking_Ns.json"
    n_filtered_loci_filename = f"{base}.without_loci_with_flanking_Ns.filtered_loci.txt"
    # final IlluminaEHv5 catalog: N-flank-filtered AND with large loci removed
    filtered_catalog_filename = f"{base}.for_illumina_eh.json"
    filtered_catalog_path = os.path.join(output_dir, filtered_catalog_filename)

    # the EH image (weisburd/str-analysis-with-expansion-hunter, pinned by digest) bakes in str_analysis, so the
    # filter runs reproducibly without a runtime pip install
    step = bp.new_step(
        name=f"Prefilter EHv5 catalog (drop flanking-N and large loci) for IlluminaEHv5: {os.path.basename(eh_catalog_path)}",
        arg_suffix="prefilter-illumina-eh-catalog-step",
        image=EH_DOCKER_IMAGE,
        cpu=2,
        memory="standard",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)
    step.command("set -ex")
    local_fasta = step.input(reference_fasta)
    step.input(reference_fasta_fai)
    local_catalog = step.input(eh_catalog_path)
    # filter 1: drop loci with >5 Ns in their +/-1000bp flanks
    step.command(
        f"python3 -m str_analysis.filter_out_loci_with_Ns_in_flanks "
        f"-R {local_fasta} --region-extension-length 1000 "
        f"-o {n_filtered_filename} "
        f"-f {n_filtered_loci_filename} "
        f"{local_catalog}")
    # filter 2: drop loci whose reference interval is >= ILLUMINA_EH_MAX_LOCUS_SIZE_BP wide (official build crashes
    # with "numIndels out of range" on large VNTRs). Quoted heredoc so the shell leaves the regex/anchors alone; the
    # f-string only fills in the file names and the size cutoff.
    step.command(
        f"""python3 <<'PYEOF'
import json, re
cat = json.load(open("{n_filtered_filename}"))
kept = []
dropped = 0
for rec in cat:
    rr = rec["ReferenceRegion"]
    regions = rr if isinstance(rr, list) else [rr]
    starts = []
    ends = []
    for reg in regions:
        m = re.match(r"^(.+):([0-9]+)-([0-9]+)$", reg)
        starts.append(int(m.group(2)))
        ends.append(int(m.group(3)))
    if max(ends) - min(starts) >= {ILLUMINA_EH_MAX_LOCUS_SIZE_BP}:
        dropped += 1
        continue
    kept.append(rec)
json.dump(kept, open("{filtered_catalog_filename}", "wt"), indent=4)
print("Dropped %d loci >= {ILLUMINA_EH_MAX_LOCUS_SIZE_BP}bp wide; kept %d" % (dropped, len(kept)))
PYEOF""")
    step.command("ls -lhrt")
    step.output(n_filtered_filename)
    step.output(n_filtered_loci_filename)
    step.output(filtered_catalog_filename)
    return step, filtered_catalog_path


def create_extract_allele_sequences_step(bp, tool_results_step, *, tool, vcf_paths, output_prefix, output_dir,
                                         reference_fasta, reference_fasta_fai):
    """Build a step that extracts the tool's per-allele sequences from the VCF it already wrote.

    This is what makes the sequence-accuracy benchmark cheap: every tool VCF is already on GCS, so nothing has to be
    re-genotyped. The step just re-reads that VCF at cpu=1, resolves each record's GT to the REF/ALT sequences, trims
    them to the truth locus interval, and writes a small side-car table that add_sequence_accuracy_columns.py joins
    onto the comparison table.

    Args:
        bp: the step_pipeline pipeline object.
        tool_results_step: the tool's genotyping/combine step, depended on so the VCF exists before this step runs.
        tool: "TRGTv5", "ATaRVa", or "HipSTR" (one of SEQUENCE_ACCURACY_TOOLS).
        vcf_paths: gs:// path(s) of the tool's output VCF(s). HipSTR writes one per catalog shard.
        output_prefix: filename prefix for the output table ("{sample_id}.{tool}").
        output_dir: the tool's {coverage}_coverage output dir.
        reference_fasta: gs:// path of the reference fasta, used to confirm each record's REF matches the reference
            genome over the locus interval (a locus that doesn't match is reported rather than silently mis-scored).
        reference_fasta_fai: gs:// path of the reference fasta .fai index.

    Returns:
        A (step, allele_sequences_path) tuple; allele_sequences_path is the gs:// path of the output table.
    """
    allele_sequences_filename = f"{output_prefix}.allele_sequences.tsv.gz"

    step = bp.new_step(
        name=f"Extract {tool} allele sequences for {os.path.basename(output_dir)}",
        arg_suffix="extract-allele-sequences-step",
        image=FILTER_VCFS_DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)

    step.depends_on(tool_results_step)

    # set -ex before the inputs, matching the other step factories here: the GSUTIL_COPY localization commands are
    # emitted at .input() time, so this makes a failed localization abort the job instead of falling through to the
    # extractor with a missing file.
    step.command("set -ex")
    local_fasta = step.input(reference_fasta)
    step.input(reference_fasta_fai)
    local_vcfs = [step.input(vcf_path) for vcf_path in vcf_paths]

    step.command(
        f"python3 -u -m str_analysis.extract_allele_sequences_from_vcf "
        f"--tool {tool} "
        f"-R {local_fasta} "
        f"-o {allele_sequences_filename} "
        + " ".join(str(local_vcf) for local_vcf in local_vcfs))
    step.command("ls -lhrt")

    step.output(allele_sequences_filename)

    return step, os.path.join(output_dir, allele_sequences_filename)


def add_tool_comparison_columns_step(bp, tool_results_step, *, tool, coverage_label, sample_id, output_dir, truth_set_genotypes_path, tool2="Truth", allele_sequences_step=None, allele_sequences_path=None, download_to_dir=None):
    tool_results_path = None
    for output_spec in tool_results_step.get_outputs():
        if output_spec.output_path.endswith(".variants.tsv.gz"):
            tool_results_path = output_spec.output_path

    if tool_results_path is None:
        raise ValueError(f"Couldn't find a .variants.tsv.gz table among the output files for {tool}: "
                         f"{[s.output_path for s in tool_results_step.get_outputs()]}")

    add_columns_step = bp.new_step(
        name=f"Add {sample_id} {tool} results columns to {tool2} table for {os.path.basename(output_dir)}",
        arg_suffix=f"add-columns-step",
        image=FILTER_VCFS_DOCKER_IMAGE,
        # cpu=2/highmem (~13GB) OOM-kills add_concordance_columns.py on the 1.6M-locus combined catalog
        # (~762k-row pandas table); cpu=4/highmem (~26GB) gives headroom.
        cpu=4,
        memory="highmem",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)

    add_columns_step.depends_on(tool_results_step)

    add_columns_step.command("set -ex")

    local_tool_results_input = add_columns_step.input(tool_results_path)
    # the truth set is the v2 genotype table; compute_truth_set_tsv_for_comparisons.py normalizes its columns and
    # carries the per-allele repeat purity (RepeatPurity: Allele 1/2) used by the purity-stratified plots
    local_truth_set_genotypes = add_columns_step.input(truth_set_genotypes_path)

    # for the tools that report allele sequences, also localize the extract step's side-car table so the
    # sequence-accuracy columns can be added between add_tool_results_columns and add_concordance_columns
    local_allele_sequences = None
    if allele_sequences_step is not None:
        add_columns_step.depends_on(allele_sequences_step)
        local_allele_sequences = add_columns_step.input(allele_sequences_path)

    add_columns_step.command(f"""python3 <<EOF
import pandas as pd
print("Adding columns to {local_tool_results_input}")
try:
    df = pd.read_table("{local_tool_results_input}", dtype=str)
except pd.errors.EmptyDataError:
    # The tool genotyped nothing (e.g. HipSTR/EH on single-end ultima reads, or LongTR on ONT before the
    # --min-mean-qual fix), so the combined variants table is an empty gzip. Keep going with an empty table:
    # add_tool_results_columns.py marks every truth-set locus as a No Call for this tool.
    df = pd.DataFrame()
df.loc[:, "Coverage"] = "{coverage_label}"
print(f"Writing {{len(df):,d}} rows to {local_tool_results_input}")
df.to_csv("{local_tool_results_input}", sep="\\t", index=False, header=True)
EOF
""")

    # matches the output name compute_truth_set_tsv_for_comparisons.py derives from the input basename
    for_comparison_filename = os.path.basename(truth_set_genotypes_path).replace(".tsv", ".for_comparison.tsv")
    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py "
               f"--output-dir . "
               f"{local_truth_set_genotypes}")

    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/add_tool_results_columns.py "
               f"--tool {tool} "
               f"{local_tool_results_input} "
               f"{for_comparison_filename} ")

    add_columns_step.command("ls -lhrt")

    local_tsv_file_path = for_comparison_filename.replace(".tsv.gz", "") + f".with_{tool}_results.tsv.gz"
    output_filename = for_comparison_filename.replace(".tsv.gz", "") + f".with_{tool}_vs_{tool2}_columns.tsv.gz"

    # add the per-allele edit distances in place, before add_concordance_columns.py melts the table into the alleles
    # table (which is what the plots read). The sequences themselves are dropped inside this script, so they never
    # reach the variants / alleles tables.
    if local_allele_sequences is not None:
        add_columns_step.command(
            f"python3 -u /str-truth-set/tool_comparison/scripts/add_sequence_accuracy_columns.py "
            f"--tool {tool} "
            f"--truth-set-genotypes {local_truth_set_genotypes} "
            f"--allele-sequences {local_allele_sequences} "
            f"{local_tsv_file_path}")
        add_columns_step.command("ls -lhrt")

    add_columns_step.command(f"python3 -u /str-truth-set/tool_comparison/scripts/add_concordance_columns.py "
               f"--tool {tool} "
               f"--compare-to {tool2} "
               f"--output-tsv {output_filename} "
               f"{local_tsv_file_path}")

    add_columns_step.command("ls -lhrt")

    add_columns_step.output(output_filename, download_to_dir=download_to_dir)
    add_columns_step.output(output_filename.replace(".tsv", ".alleles.tsv"))

    return add_columns_step


def create_plot_tool_accuracy_steps(bp, add_columns_step, *, tool, coverage_label, sequencing_data_type, sample_id, output_dir):
    plot_tool_accuracy_step = bp.new_step(
        name=f"Plot {sample_id} {tool} accuracy for {os.path.basename(output_dir)}",
        arg_suffix=f"plot-accuracy-step",
        image=FILTER_VCFS_DOCKER_IMAGE,
        cpu=1,
        memory="highmem",
        storage="20Gi",
        output_dir=output_dir)

    local_variants_tsv, local_alleles_tsv = plot_tool_accuracy_step.use_previous_step_outputs_as_inputs(add_columns_step)

    plot_tool_accuracy_step.command("set -ex")

    # The plot script stratifies internally by purity bin, IsPureRepeat, and no-call loci, so each invocation produces
    # many svg files. They're all captured below with a wildcard output.
    for min_motif_size, max_motif_size in MOTIF_SIZE_BINS:
        plot_tool_accuracy_step.command(
            f"python3 -u /str-truth-set/figures_and_tables/plot_tool_accuracy_by_allele_size.py "
            "--verbose "
            f"--tool {tool} "
            f"--coverage {coverage_label} "
            f"--sequencing-data-type {sequencing_data_type} "
            "--q-threshold 0 "
            f"--min-motif-size {min_motif_size} "
            f"--max-motif-size {max_motif_size} "
            "--image-type svg "
            "--show-title "
            f"{local_alleles_tsv} ")
        plot_tool_accuracy_step.command("ls -lhrt")

    # Also generate the unstratified "all motif sizes" plot (".all_motifs", surfaced as the "all" bin in the viewer)
    # by running the script in its default mode with no --min/--max-motif-size. --all-motifs-only restricts that mode to
    # the all_motifs plot so it no longer emits the unused .2bp_motifs/.3to6bp_motifs/.25to50bp_motifs coarse-bin svgs.
    plot_tool_accuracy_step.command(
        f"python3 -u /str-truth-set/figures_and_tables/plot_tool_accuracy_by_allele_size.py "
        "--verbose "
        f"--tool {tool} "
        f"--coverage {coverage_label} "
        f"--sequencing-data-type {sequencing_data_type} "
        "--q-threshold 0 "
        "--all-motifs-only "
        "--image-type svg "
        "--show-title "
        f"{local_alleles_tsv} ")
    plot_tool_accuracy_step.command("ls -lhrt")

    # gzip each svg in place (keeping the .svg name) and upload it with the content headers set inline on the cp via
    # step.output's content_encoding/content_type, so the browser fetches the .svg URL, receives gzipped bytes (it
    # sends Accept-Encoding: gzip), and renders the svg directly. Uploaded with a wildcard via gcloud storage cp
    # (Delocalize.GSUTIL_COPY), since Delocalize.COPY needs an explicit filename per output. Setting the headers at
    # upload time avoids a separate "gcloud storage objects update" step, which intermittently hit "HTTPError 409 ...
    # edited during the operation" when updating freshly-created objects.
    # For the tools that report allele sequences, also plot how close those sequences are to the truth sequences.
    # add_sequence_accuracy_columns.py put the per-allele edit distances into the same alleles table, so this reuses
    # the input that's already localized. One invocation per motif bin (each emits 2 metrics x 3 genotype subsets),
    # matching the motif loop above, plus the unstratified all-motifs invocation.
    if tool in SEQUENCE_ACCURACY_TOOLS:
        for min_motif_size, max_motif_size in MOTIF_SIZE_BINS + [(None, None)]:
            plot_tool_accuracy_step.command(
                f"python3 -u /str-truth-set/figures_and_tables/plot_tool_sequence_accuracy.py "
                f"--tool {tool} "
                f"--coverage {coverage_label} "
                f"--sequencing-data-type {sequencing_data_type} "
                + (f"--min-motif-size {min_motif_size} --max-motif-size {max_motif_size} "
                   if min_motif_size is not None else "")
                + "--image-type svg "
                "--show-title "
                f"{local_alleles_tsv} ")
            plot_tool_accuracy_step.command("ls -lhrt")

    plot_tool_accuracy_step.command('for f in tool_accuracy_by_true_allele_size.*.svg; do gzip "$f"; mv "$f.gz" "$f"; done')
    plot_tool_accuracy_step.output(
        "tool_accuracy_by_true_allele_size.*.svg",
        delocalize_by=Delocalize.GSUTIL_COPY,
        content_encoding="gzip",
        content_type="image/svg+xml")

    if tool in SEQUENCE_ACCURACY_TOOLS:
        plot_tool_accuracy_step.command(
            'for f in tool_sequence_accuracy_by_true_allele_size.*.svg; do gzip "$f"; mv "$f.gz" "$f"; done')
        plot_tool_accuracy_step.output(
            "tool_sequence_accuracy_by_true_allele_size.*.svg",
            delocalize_by=Delocalize.GSUTIL_COPY,
            content_encoding="gzip",
            content_type="image/svg+xml")

    #plot_tool_accuracy_step.command(f"python3 -u /str-truth-set/figures_and_tables/plot_tool_accuracy_vs_Q.py "
    #                                "--verbose "
    #                                f"--coverage {coverage_label} "
    #                                f"--min-motif-size {min_motif_size} "
    #                                f"--max-motif-size {max_motif_size} "
    #                                "--genotype all "
    #                                "--show-no-call-loci "
    #                                "--image-type svg "
    #                                "--show-title "
    #                                f"{local_alleles_tsv} ")

    plot_tool_accuracy_step.command("ls -lhrt")

    return plot_tool_accuracy_step

#
#    python3 -u plot_tool_accuracy_percent_exactly_right.py --show-title --output-dir ${output_dir}
#    python3 -u plot_tool_accuracy_percent_exactly_right.py --show-title --exclude-hipstr-no-call-loci --output-dir ${output_dir}
#    python3 -u plot_tool_accuracy_by_motif_size.py --output-dir ${output_dir}
#
#    # generate tool accuracy vs Q and tool accuracy by num repeats plots
#    python3 figures_pipeline.py --force --batch-size 25
#
#    gsutil -m cp -r gs://str-truth-set/hg38/figures/accuracy_vs_Q .
#    gsutil -m cp -r gs://str-truth-set/hg38/figures/accuracy_by_allele_size  .
#
#    #python3 -u plot_tool_accuracy_vs_Q.py --verbose
#    #python3 -u plot_tool_accuracy_by_allele_size.py --verbose
#
#    ./generate_figure_panels_for_paper.sh
#
#    TODO plot fraction of loci that are polymorphic by size in reference

if __name__ == "__main__":
    main()
