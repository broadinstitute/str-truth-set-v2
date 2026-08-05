"""Generate the seq_rankings.json files that power the Accuracy Summary row in the tool-comparison viewer's
Sequence Accuracy section.

For every bin combination that section exposes (motif size x genotype x coverage), this recomputes each
sequence-reporting tool's exact-match rate -- the "% green" of the stacked histograms, i.e. the "...got N out
of M allele sequences (P%) exactly right" number in the plot titles -- and writes one JSON per (sample,
sequencing data type) that the front-end loads to render the row.

Sibling of generate_tool_ranking_json.py (same idea, for the repeat-count Accuracy section above it). This one
has a narrower bin space -- no purity / chromosome-class / hide-no-call facets, since the sequence-accuracy
plots don't expose those -- and no metric dimension, since exact-match (edit distance == 0) is identical
whether it's read off the raw or the normalized distance column.

Input  (already produced by run_genotyping_tools.py):
    gs://.../tool_results/{sample}/{dtype}/{tool}/{cov}_coverage/*.with_{tool}_vs_Truth_columns.alleles.tsv.gz
Output (one per sample x dtype):
    gs://.../tool_results/{sample}/{dtype}/seq_rankings.json

JSON shape -- keyed by "{motif_token}.{genotype_token}.{coverage_label}", with a trailing
".exclude_no_call_loci" on the Hide "No Call" loci variant (see buildSeqRankingKey() in
docs/tool_comparison_viewer.html). Each value is the tools sorted by descending exact-match percent (tie-
broken by tool name for determinism):

    {"<bin key>": [{"tool": "TRGTv5", "pct": 88.2, "exact": 12345, "total": 13985}, ...], ...}
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "str-truth-set", "figures_and_tables"))
from plot_tool_sequence_accuracy import SEQUENCE_ACCURACY_TOOLS

DEFAULT_TOOL_RESULTS_ROOT = "gs://str-truth-set-v2/tool_results"

# Motif-size bins == run_genotyping_tools.py MOTIF_SIZE_BINS, plus the unstratified "all_motifs" bin (mirrors
# generate_tool_ranking_json.py's MOTIF_BINS; token is the plot script's ".{...}_motifs" filename segment).
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

# Genotype subset == plot_tool_sequence_accuracy.py's genotype_subset loop; token is the filename segment.
GENOTYPE_BINS = [("all_genotypes", "all"), ("HET", "HET"), ("HOM", "HOM")]


def build_bin_key(motif_token, genotype_token, coverage_label, exclude_no_call_loci=False):
    """Returns the key the viewer looks up: "{motif_token}.{genotype_token}.{coverage_label}", plus a trailing
    ".exclude_no_call_loci" for the Hide "No Call" loci variant. See buildSeqRankingKey() in the viewer."""
    return (f"{motif_token}.{genotype_token}.{coverage_label}"
            + (".exclude_no_call_loci" if exclude_no_call_loci else ""))


def compute_tool_accuracy(tsv_path, tool, coverage_label):
    """Recomputes exact-match accuracy for one tool over every viewer bin combination.

    Replicates the filter chain in plot_tool_sequence_accuracy.py's main()/generate_all_plots() (drop rows
    with no true allele size, restrict to this coverage, motif/genotype filters, HipSTR's >9bp skip, the
    <10-alleles skip) so each returned pct equals the plotted "% exactly right". Exact-match is
    `distance == 0`, which NaN (no-call / unscored) never satisfies -- the same classification
    bin_edit_distance() gives it via NO_CALL_LABEL, just without the row-wise apply.

    Args:
        tsv_path (str): local path to a *.with_{tool}_vs_Truth_columns.alleles.tsv.gz table.
        tool (str): the tool name (also the suffix of its per-allele columns).
        coverage_label (str): coverage token (e.g. "30x", "37x"), used as the Coverage value and bin key.

    Returns:
        dict: {bin_key: {"exact": int, "total": int, "pct": float}}
    """
    df = pd.read_table(tsv_path)
    df = df[~df["DiffFromRefRepeats: Allele: Truth"].isna()]
    if "Coverage" not in df.columns:
        df["Coverage"] = coverage_label
    df = df[df["Coverage"] == coverage_label]

    distance_column = f"SequenceEditDistance: Allele: {tool}"
    if distance_column not in df.columns or len(df) == 0:
        return {}

    result = {}
    for motif_token, min_motif, max_motif in MOTIF_BINS:
        # HipSTR doesn't support motifs larger than 9bp (plot script skip condition): all_motifs and any bin
        # whose whole range is >9bp produce a placeholder image with no percent, so omit them.
        if tool == "HipSTR" and min_motif is not None and min_motif > 9 and (max_motif is None or max_motif > 9):
            continue
        if min_motif is None:
            df_motif = df[(2 <= df["MotifSize"]) & (df["MotifSize"] <= 50)]
        else:
            df_motif = df[(min_motif <= df["MotifSize"]) & (df["MotifSize"] <= max_motif)]
        for genotype_token, genotype in GENOTYPE_BINS:
            if genotype == "all":
                df_geno = df_motif
            elif genotype == "HET":
                df_geno = df_motif[df_motif["SummaryString"].str.contains(":HET")]
            else:
                df_geno = df_motif[df_motif["SummaryString"].str.contains(":HOM")]
            # Both no-call variants, matching the plot script's exclude_no_call_loci loop. Dropping the
            # unscored (NaN) alleles shrinks the denominator but never the numerator, so the two keys report
            # different percentages for the same locus set.
            for exclude_no_call_loci in (False, True):
                df_plot = df_geno[~df_geno[distance_column].isna()] if exclude_no_call_loci else df_geno
                if len(df_plot) < 10:  # plot script skip_condition2: not enough alleles to draw a histogram
                    continue
                exact = int((df_plot[distance_column] == 0).sum())
                result[build_bin_key(motif_token, genotype_token, coverage_label, exclude_no_call_loci)] = {
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
    """Builds and uploads seq_rankings.json for one (sample, data type) by scanning its SEQUENCE_ACCURACY_TOOLS."""
    dtype = basename(dtype_uri)
    rankings = {}  # bin_key -> {tool -> {exact, total, pct}}
    for tool_uri in [u for u in gcs_ls(dtype_uri) if u.endswith("/")]:
        tool = basename(tool_uri)
        if tool not in SEQUENCE_ACCURACY_TOOLS:
            continue
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
        print(f"{sample}/{dtype}: no sequence-accuracy data; not writing seq_rankings.json")
        return

    output = {
        key: [{"tool": t, "pct": s["pct"], "exact": s["exact"], "total": s["total"]}
              for t, s in sorted(tool_stats.items(), key=lambda kv: (-kv[1]["pct"], kv[0]))]
        for key, tool_stats in rankings.items()
    }
    local_json = os.path.join(tmp_dir, f"{sample}.{dtype}.seq_rankings.json")
    with open(local_json, "w") as f:
        json.dump(output, f, separators=(",", ":"), sort_keys=True)
    dest = f"{dtype_uri.rstrip('/')}/seq_rankings.json"
    if dry_run:
        print(f"{sample}/{dtype}: [dry-run] {len(output):,d} bins -> {dest} ({os.path.getsize(local_json):,d} bytes)")
    else:
        # no-cache so the public storage.googleapis.com edge revalidates every fetch; otherwise the default
        # "public, max-age=3600" makes a regenerated seq_rankings.json invisible to the viewer for up to an hour.
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
    p.add_argument("--dry-run", action="store_true", help="Compute and report sizes but don't upload seq_rankings.json")
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
