"""Generate the rankings.json files that power the per-bin tool-ranking row in the tool-comparison viewer.

For every bin combination the viewer exposes (motif size x genotype x coverage x purity x chromosome class x
"Hide No Call loci"), this recomputes each tool's exact-match accuracy -- the "% green" of the stacked
histograms, i.e. the "...got N out of M alleles (P%) exactly right" number in the plot titles -- and writes
one JSON per (sample, sequencing data type) that the front-end loads to render the ranking row.

Zero-drift by construction: the bin/filter/exact-match logic is IMPORTED from
str-truth-set/figures_and_tables/plot_tool_accuracy_by_allele_size.py (define_hue_column, PURITY_BINS,
CHROM_CLASSES, chrom_class_of_locus_id, the NO_CALL/Same labels), not reimplemented. Only the per-plot loop
and the JSON emission are new here, so each pct equals the value the plot script bakes into the SVG title.

Input  (already produced by run_genotyping_tools.py):
    gs://.../tool_results/{sample}/{dtype}/{tool}/{cov}_coverage/*.with_{tool}_vs_Truth_columns.alleles.tsv.gz
Output (one per sample x dtype):
    gs://.../tool_results/{sample}/{dtype}/rankings.json

JSON shape -- keyed by the SVG filename "middle" (everything between the
"tool_accuracy_by_true_allele_size." prefix and the ".{tool}.svg" suffix), which is exactly the key the
viewer builds from its controls (see buildBinKey() in docs/tool_comparison_viewer.html). Each value is the
tools sorted by descending exact-match percent (tie-broken by tool name for determinism):

    {"<bin key>": [{"tool": "TRGTv3", "pct": 96.2, "exact": 12345, "total": 12834}, ...], ...}

NOTE: the front-end fetch()es these cross-origin from storage.googleapis.com, so the bucket needs a CORS
policy that allows GET from the viewer's origin (the SVGs load via <img> and don't need it, but fetch() does).
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

# Headless matplotlib: importing the plot module pulls in matplotlib/seaborn, which we never draw with here.
os.environ.setdefault("MPLBACKEND", "Agg")
import pandas as pd

# Reuse the plot script's bin-defining logic verbatim so the ranking matches the plotted "% exactly right".
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "str-truth-set", "figures_and_tables"))
import plot_tool_accuracy_by_allele_size as P

DEFAULT_TOOL_RESULTS_ROOT = "gs://str-truth-set-v2/tool_results"

# Motif-size bins surfaced by the viewer == run_genotyping_tools.py MOTIF_SIZE_BINS (each emitted by a
# --min/--max plot invocation) plus the "all_motifs" default-mode bin (2 <= MotifSize <= 50). The token is the
# ".{...}_motifs" filename segment; (min, max) None means the unstratified all_motifs filter.
MOTIF_BINS = [
    ("all_motifs", None, None),
    ("1to1bp_motifs", 1, 1),
    ("2to2bp_motifs", 2, 2),
    ("3to3bp_motifs", 3, 3),
    ("4to4bp_motifs", 4, 4),
    ("5to5bp_motifs", 5, 5),
    ("6to6bp_motifs", 6, 6),
    ("2to6bp_motifs", 2, 6),
    ("7to24bp_motifs", 7, 24),
    ("25to1000bp_motifs", 25, 1000),
]

# Genotype subset == the plot script's genotype_subset loop; token is the filename segment.
GENOTYPE_BINS = [("all_genotypes", "all"), ("HET", "HET"), ("HOM", "HOM")]


def build_bin_key(motif_token, genotype_token, coverage_label, purity_name, chrom_class, hide_no_call):
    """Returns the SVG filename between the "tool_accuracy_by_true_allele_size." prefix and ".{tool}.svg".

    Mirrors plot_tool_accuracy_by_allele_size.py output_image_filename (minus the tool/extension) and the
    viewer's buildBinKey(), so the front-end looks each combination up directly.
    """
    key = f"{motif_token}.{genotype_token}.{coverage_label}"
    if purity_name and purity_name != "all":
        key += f".purity_{purity_name}"
    if chrom_class and chrom_class != "all":
        key += f".{chrom_class}"
    if hide_no_call:
        key += ".exclude_no_call_loci"
    return key


def compute_tool_accuracy(tsv_path, tool, coverage_label):
    """Recomputes exact-match accuracy for one tool over every viewer bin combination.

    Replicates the filter/skip chain in plot_tool_accuracy_by_allele_size.py for the way the pipeline invokes
    it (q-threshold forced to 0 -> no Filtered loci; both no-call-kept and no-call-excluded variants), so each
    returned pct equals the plotted "% exactly right".

    Args:
        tsv_path (str): local path to a *.with_{tool}_vs_Truth_columns.alleles.tsv.gz table.
        tool (str): the tool name (also the suffix of its per-allele columns).
        coverage_label (str): coverage token (e.g. "30x", "24G"), used as the Coverage value and bin key.

    Returns:
        dict: {bin_key: {"exact": int, "total": int, "pct": float}}
    """
    df = pd.read_table(tsv_path)
    df = df[~df["DiffFromRefRepeats: Allele: Truth"].isna()]
    if "IsFoundInReference" in df.columns:
        df = df[df["IsFoundInReference"]]
    if "PositiveOrNegative" in df.columns:
        df = df[df["PositiveOrNegative"] == "positive"]
    if "Coverage" not in df.columns:
        df["Coverage"] = coverage_label  # single-coverage table; matches the plot script's fallback

    P.define_hue_column(df, tool)  # adds the "DiffRepeats: Allele: {tool} - Truth (bin)" hue column
    bin_column = f"DiffRepeats: Allele: {tool} - Truth (bin)"

    if P.PURITY_COLUMN in df.columns:
        df.loc[:, P.PURITY_COLUMN] = pd.to_numeric(df[P.PURITY_COLUMN], errors="coerce")
        purity_bins = [("all", None, None)] + P.PURITY_BINS
    else:
        purity_bins = [("all", None, None)]
    df["ChromClass"] = df["LocusId"].apply(P.chrom_class_of_locus_id)

    result = {}
    for chrom_class in P.CHROM_CLASSES:
        df_chrom = df if chrom_class == "all" else df[df["ChromClass"] == chrom_class]
        for purity_name, purity_low, purity_high in purity_bins:
            if purity_low is None:
                df_purity = df_chrom
            else:
                df_purity = df_chrom[(df_chrom[P.PURITY_COLUMN] >= purity_low) & (df_chrom[P.PURITY_COLUMN] < purity_high)]
            for motif_token, min_motif, max_motif in MOTIF_BINS:
                # HipSTR doesn't support motifs larger than 9bp (plot script skip_condition1): all_motifs and any
                # bin whose whole range is >9bp produce a placeholder image with no percent, so omit them.
                if tool == "HipSTR" and (min_motif is None or min_motif > 9) and (max_motif is None or max_motif > 9):
                    continue
                if min_motif is None:
                    df_motif = df_purity[(2 <= df_purity["MotifSize"]) & (df_purity["MotifSize"] <= 50)]
                else:
                    df_motif = df_purity[(min_motif <= df_purity["MotifSize"]) & (df_purity["MotifSize"] <= max_motif)]
                for genotype_token, genotype in GENOTYPE_BINS:
                    if genotype == "all":
                        df_geno = df_motif
                    elif genotype == "HET":
                        df_geno = df_motif[df_motif["SummaryString"].str.contains(":HET")]
                    else:
                        df_geno = df_motif[df_motif["SummaryString"].str.contains(":HOM")]
                    for hide_no_call in (False, True):
                        df_plot = df_geno if not hide_no_call else df_geno[df_geno[bin_column] != P.NO_CALL_LABEL]
                        if len(df_plot) < 10:  # plot script skip_condition2: not enough alleles to draw a histogram
                            continue
                        exact = int((df_plot[bin_column] == "0").sum() + (df_plot[bin_column] == P.NO_DIFFERENCE_LABEL).sum())
                        result[build_bin_key(motif_token, genotype_token, coverage_label, purity_name, chrom_class, hide_no_call)] = {
                            "exact": exact, "total": len(df_plot), "pct": round(100.0 * exact / len(df_plot), 1)}
    return result


def gcs_ls(gcs_dir):
    """Lists immediate children of a gs:// "directory"; returns full gs:// URIs (empty on a missing path)."""
    result = subprocess.run(["gcloud", "storage", "ls", gcs_dir.rstrip("/") + "/"], capture_output=True, text=True)
    if result.returncode != 0:
        return []
    return [line.strip() for line in result.stdout.splitlines() if line.strip().startswith("gs://")]


def basename(uri):
    """Last path segment of a gs:// URI, ignoring any trailing slash."""
    return uri.rstrip("/").split("/")[-1]


def generate_for_sample_data_type(sample, dtype_uri, tmp_dir, dry_run):
    """Builds and uploads rankings.json for one (sample, data type) by scanning every tool/coverage under it."""
    dtype = basename(dtype_uri)
    rankings = {}  # bin_key -> {tool -> {exact, total, pct}}
    for tool_uri in [u for u in gcs_ls(dtype_uri) if u.endswith("/")]:
        tool = basename(tool_uri)
        for cov_uri in [u for u in gcs_ls(tool_uri) if u.rstrip("/").endswith("_coverage")]:
            coverage_label = basename(cov_uri)[:-len("_coverage")]
            alleles_uris = [u for u in gcs_ls(cov_uri) if u.endswith(f".with_{tool}_vs_Truth_columns.alleles.tsv.gz")]
            if not alleles_uris:
                print(f"  {sample}/{dtype}/{tool}/{coverage_label}: no alleles tsv; skipping")
                continue
            local_tsv = os.path.join(tmp_dir, f"{sample}.{dtype}.{tool}.{coverage_label}.alleles.tsv.gz")
            try:
                subprocess.run(["gcloud", "storage", "cp", alleles_uris[0], local_tsv], check=True, capture_output=True, text=True)
                per_bin = compute_tool_accuracy(local_tsv, tool, coverage_label)
            except Exception as e:
                print(f"  {sample}/{dtype}/{tool}/{coverage_label}: ERROR ({e}); skipping")
                continue
            finally:
                if os.path.exists(local_tsv):
                    os.remove(local_tsv)
            for key, stats in per_bin.items():
                rankings.setdefault(key, {})[tool] = stats
            print(f"  {sample}/{dtype}/{tool}/{coverage_label}: {len(per_bin):,d} bins")

    if not rankings:
        print(f"{sample}/{dtype}: no data; not writing rankings.json")
        return

    output = {
        key: [{"tool": t, "pct": s["pct"], "exact": s["exact"], "total": s["total"]}
              for t, s in sorted(tool_stats.items(), key=lambda kv: (-kv[1]["pct"], kv[0]))]
        for key, tool_stats in rankings.items()
    }
    local_json = os.path.join(tmp_dir, f"{sample}.{dtype}.rankings.json")
    with open(local_json, "w") as f:
        json.dump(output, f, separators=(",", ":"), sort_keys=True)
    dest = f"{dtype_uri.rstrip('/')}/rankings.json"
    if dry_run:
        print(f"{sample}/{dtype}: [dry-run] {len(output):,d} bins -> {dest} ({os.path.getsize(local_json):,d} bytes)")
    else:
        # no-cache so the public storage.googleapis.com edge revalidates every fetch; otherwise the default
        # "public, max-age=3600" makes a regenerated rankings.json invisible to the viewer for up to an hour.
        subprocess.run(["gcloud", "storage", "cp", "--content-type", "application/json",
                        "--cache-control", "no-cache, max-age=0", local_json, dest], check=True)
        print(f"{sample}/{dtype}: wrote {len(output):,d} bins -> {dest}")
    os.remove(local_json)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tool-results-root", default=DEFAULT_TOOL_RESULTS_ROOT,
                   help="gs:// root the genotyping pipeline writes to (default: %(default)s)")
    p.add_argument("--sample", action="append", help="Sample id to process (repeatable; default: HG002)")
    p.add_argument("--data-type", action="append",
                   help="Restrict to this sequencing data type (repeatable; default: all found under the sample)")
    p.add_argument("--dry-run", action="store_true", help="Compute and report sizes but don't upload rankings.json")
    args = p.parse_args()

    samples = args.sample or ["HG002"]
    root = args.tool_results_root.rstrip("/")
    with tempfile.TemporaryDirectory() as tmp_dir:
        for sample in samples:
            for dtype_uri in [u for u in gcs_ls(f"{root}/{sample}") if u.endswith("/")]:
                if args.data_type and basename(dtype_uri) not in args.data_type:
                    continue
                generate_for_sample_data_type(sample, dtype_uri, tmp_dir, args.dry_run)


if __name__ == "__main__":
    main()
