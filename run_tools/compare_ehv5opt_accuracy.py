"""Stratified before/after accuracy comparison for the EHv5-bw2-optimized tool.

Reads the per-(sample,data_type) rankings.json files from two directories (baseline vs new) and, for every
stratification bin, extracts the EHv5-bw2-optimized exact-match percentage. Reports the headline per-dataset
numbers, the largest improvements/regressions, and breakdowns by stratification dimension so it is easy to see
where the latest bw2/ExpansionHunter optimized-streaming commits changed accuracy.

Usage:
  python3 compare_ehv5opt_accuracy.py --old-dir /tmp/ehv5opt_rankings_before --new-dir /tmp/ehv5opt_rankings_after
"""
import argparse
import glob
import json
import os

TOOL = "EHv5-bw2-optimized"


def load_tool_bins(path):
    """Returns {bin_key -> {pct,exact,total}} holding only the EHv5-bw2-optimized entry of each bin."""
    out = {}
    with open(path) as f:
        data = json.load(f)
    for bin_key, ranked in data.items():
        for entry in ranked:
            if entry["tool"] == TOOL:
                out[bin_key] = {"pct": entry["pct"], "exact": entry["exact"], "total": entry["total"]}
                break
    return out


def parse_key(bin_key):
    """Splits a bin key into its stratification dimensions.

    Key layout: {motif}.{genotype}.{coverage}[.purity_{...}][.{chrom}][.exclude_no_call_loci]
    Returns dict with motif, genotype, coverage, purity, chrom, no_call.
    """
    parts = bin_key.split(".")
    motif, genotype, coverage = parts[0], parts[1], parts[2]
    rest = parts[3:]
    no_call = "exclude_no_call_loci" in rest
    rest = [p for p in rest if p != "exclude_no_call_loci"]
    purity = "all"
    chrom = "all"
    i = 0
    while i < len(rest):
        if rest[i] == "purity_0":  # purity token is split by '.' into purity_0 + 65_to_0 + 75 ; rejoin
            purity = ".".join(rest[i:i + 3])
            i += 3
        elif rest[i] in ("autosomes", "chrX", "chrY"):
            chrom = rest[i]
            i += 1
        else:
            i += 1
    return {"motif": motif, "genotype": genotype, "coverage": coverage,
            "purity": purity, "chrom": chrom, "no_call": "exclude" if no_call else "include"}


def fmt_delta(d):
    return f"{d:+.1f}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--old-dir", required=True)
    p.add_argument("--new-dir", required=True)
    p.add_argument("--min-total", type=int, default=50,
                   help="Ignore bins with fewer than this many alleles when ranking changes (noise filter).")
    p.add_argument("--top", type=int, default=25)
    args = p.parse_args()

    datasets = []  # (sample, dtype, old_bins, new_bins)
    for old_path in sorted(glob.glob(os.path.join(args.old_dir, "*.rankings.json"))):
        name = os.path.basename(old_path)[:-len(".rankings.json")]
        new_path = os.path.join(args.new_dir, os.path.basename(old_path))
        if not os.path.exists(new_path):
            print(f"!! no new rankings.json for {name}; skipping")
            continue
        datasets.append((name, load_tool_bins(old_path), load_tool_bins(new_path)))

    # ---- Headline: unstratified per-coverage accuracy (all_motifs.all_genotypes, purity all, chrom all, include no-call)
    print("=" * 100)
    print("HEADLINE  EHv5-bw2-optimized exact-match %  (all motifs, all genotypes, all chrom, include no-call loci)")
    print("=" * 100)
    print(f"{'dataset':<34} {'coverage':<8} {'old%':>7} {'new%':>7} {'Δ':>7} {'old exact/total':>22} {'new exact/total':>22}")
    for name, old, new in datasets:
        covs = sorted({parse_key(k)["coverage"] for k in set(old) | set(new)
                       if parse_key(k)["motif"] == "all_motifs" and parse_key(k)["genotype"] == "all_genotypes"
                       and parse_key(k)["purity"] == "all" and parse_key(k)["chrom"] == "all"
                       and parse_key(k)["no_call"] == "include"})
        for cov in covs:
            key = f"all_motifs.all_genotypes.{cov}"
            o, n = old.get(key), new.get(key)
            if not o and not n:
                continue
            op = o["pct"] if o else float("nan")
            np_ = n["pct"] if n else float("nan")
            d = (np_ - op) if (o and n) else float("nan")
            oet = f"{o['exact']}/{o['total']}" if o else "-"
            net = f"{n['exact']}/{n['total']}" if n else "-"
            print(f"{name:<34} {cov:<8} {op:>7.1f} {np_:>7.1f} {fmt_delta(d) if o and n else '   n/a':>7} {oet:>22} {net:>22}")

    # ---- Overall change tally + biggest movers (matched bins only, total>=min-total in BOTH)
    print("\n" + "=" * 100)
    print(f"CHANGE TALLY across all matched bins (both sides total >= {args.min_total})")
    print("=" * 100)
    all_rows = []  # (delta, name, bin_key, o, n)
    for name, old, new in datasets:
        for key in set(old) & set(new):
            o, n = old[key], new[key]
            if o["total"] >= args.min_total and n["total"] >= args.min_total:
                all_rows.append((round(n["pct"] - o["pct"], 1), name, key, o, n))
    improved = [r for r in all_rows if r[0] > 0]
    worsened = [r for r in all_rows if r[0] < 0]
    same = [r for r in all_rows if r[0] == 0]
    print(f"matched bins: {len(all_rows)}   improved: {len(improved)}   worsened: {len(worsened)}   unchanged: {len(same)}")
    if all_rows:
        mean_d = sum(r[0] for r in all_rows) / len(all_rows)
        print(f"mean Δ (pct points, unweighted over bins): {mean_d:+.2f}")

    def show(rows, title):
        print(f"\n--- {title} ---")
        print(f"{'Δ':>7} {'old%':>6} {'new%':>6} {'dataset':<30} {'bin'}")
        for d, name, key, o, n in rows:
            print(f"{fmt_delta(d):>7} {o['pct']:>6.1f} {n['pct']:>6.1f} {name:<30} {key}")

    show(sorted(worsened, key=lambda r: r[0])[:args.top], f"TOP {args.top} REGRESSIONS")
    show(sorted(improved, key=lambda r: -r[0])[:args.top], f"TOP {args.top} IMPROVEMENTS")

    # ---- By-dimension breakdown: mean Δ grouped by one dimension at a time (only the "all" slices of the others
    #      to keep strata non-overlapping where possible). Reports unweighted mean Δ and counts.
    print("\n" + "=" * 100)
    print("MEAN Δ BY DIMENSION (matched bins, both sides total >= min-total)")
    print("=" * 100)
    for dim in ("coverage", "chrom", "genotype", "motif", "purity", "no_call"):
        buckets = {}
        for d, name, key, o, n in all_rows:
            v = parse_key(key)[dim]
            buckets.setdefault(v, []).append(d)
        print(f"\n[{dim}]")
        for v, ds in sorted(buckets.items(), key=lambda kv: sum(kv[1]) / len(kv[1])):
            print(f"  {v:<28} n={len(ds):>5}  meanΔ={sum(ds)/len(ds):+6.2f}  worsened={sum(1 for x in ds if x<0):>5}  improved={sum(1 for x in ds if x>0):>5}")


if __name__ == "__main__":
    main()
