"""Generate per-tool repeat catalogs from a filter_vcf_to_tandem_repeats genotype output.

This is a str-truth-set-v2 port of
str-truth-set/tool_comparison/scripts/convert_truth_set_to_variant_catalogs.py, adapted to
read the new-format truth set table (<sample>.tandem_repeat_genotypes.tsv.gz) directly.

It extracts the positive (variant) loci -- those where the genotyped short or long allele
repeat count differs from the reference -- and writes EHv5, GangSTR, HipSTR, TRGT, and LongTR
catalogs plus a plain positive_loci.bed.gz (consumed by inquiSTR). Output filenames match the
globs that run_genotyping_tools.py expects under --filter-vcf-dir/<sample>/:

    <prefix>.EHv5.{NNN}_of_{NNN}.json
    <prefix>.GangSTR.{NNN}_of_{NNN}.bed
    <prefix>.HipSTR.{NNN}_of_{NNN}.bed
    <prefix>.LongTR.001_of_001.bed
    <prefix>.TRGT_repeat_catalog.bed
    <prefix>.bed.gz                       (inquiSTR)

Positive (variant) loci are the default and carry no label; negative loci, if ever added,
would be labeled explicitly. Negative loci are NOT generated here.
"""

import argparse
import json
import os
import pandas as pd

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", default=".", help="Directory where to write the catalogs")
    p.add_argument("--output-filename-prefix", help="Output filename prefix (defaults to the input "
                   "filename up to the first '.', e.g. 'HG002')")
    p.add_argument("--only", choices=["eh", "gangstr", "hipstr", "trgt", "longtr"], action="append",
                   help="Only generate catalogs for the specified tool(s). Can be repeated.")
    p.add_argument("--gangstr-loci-per-run", type=int, default=100000, help="GangSTR/HipSTR shard size. "
                   "The positive loci are split into bed files of this size.")
    p.add_argument("--expansion-hunter-loci-per-run", type=int, default=1000, help="ExpansionHunter (IlluminaEHv5) "
                   "shard size. Positive loci are split into catalog json files of this size; a single unsharded "
                   "catalog (001_of_001.json) is also written when there is more than one shard.")
    p.add_argument("genotypes_tsv_path", help="Path of a local <sample>.tandem_repeat_genotypes.tsv(.gz) file "
                   "produced by 'filter_vcf_to_tandem_repeats genotype'")
    return p.parse_args()


def generate_set_of_positive_loci(df):
    """Return the set of (chrom, start_0based, end_1based, motif) tuples for variant loci.

    A locus is positive (variant) when its genotyped short or long allele repeat count differs
    from the reference repeat count. Rows whose repeat counts can't be parsed as numbers, or whose
    genotype is missing (no-call: empty repeat-count cells, read by pandas as NaN), are skipped.

    Args:
        df (pandas.DataFrame): rows from a *.tandem_repeat_genotypes.tsv(.gz) table.

    Returns:
        set: tuples of (chrom, start_0based, end_1based, motif).
    """
    positive_loci = set()
    skipped = 0
    for _, row in df.iterrows():
        try:
            ref = float(row.NumRepeatsInReference)
            short = float(row.NumRepeatsShortAllele)
            long = float(row.NumRepeatsLongAllele)
        except (ValueError, TypeError):
            skipped += 1
            continue
        if pd.isna(ref) or pd.isna(short) or pd.isna(long):
            skipped += 1  # no-call / missing genotype (empty cells parse to NaN, which float() doesn't reject)
            continue
        if short == ref and long == ref:
            continue  # hom-ref -> not a positive locus
        # Note: unlike the v1 truth set, the new format allows loci that span a non-integer number
        # of repeats (partial repeats in VNTRs), so no integer-multiple check is applied here.
        positive_loci.add((row.Chrom, int(row.Start0Based), int(row.End), row.Motif))

    if skipped:
        print(f"Skipped {skipped:,d} rows with unparseable repeat counts")
    return positive_loci


def write_expansion_hunter_variant_catalogs(locus_set, output_path_prefix, loci_per_run):
    """Write EHv5 variant catalog json(s), sorted by canonical motif.

    Writes sharded catalogs `<prefix>.{NNN}_of_{MMM}.json` of up to loci_per_run loci each, used by the
    unoptimized IlluminaEHv5 build (which needs parallel shards to avoid running out of memory). When there is
    more than one shard, also writes the single unsharded `<prefix>.001_of_001.json` used by the streaming
    EHv5 / EHv5-bw2-optimized variants and by vamos. With only one shard the lone `001_of_001` file serves both.
    """
    variant_catalog = []
    for unmodified_chrom, start_0based, end_1based, motif in sorted(
            locus_set, key=lambda x: compute_canonical_motif(x[3], include_reverse_complement=True)):
        chrom = unmodified_chrom.replace("chr", "")
        variant_catalog.append({
            "LocusId": f"{chrom}-{start_0based}-{end_1based}-{motif}",
            "ReferenceRegion": f"{unmodified_chrom}:{start_0based}-{end_1based}",
            "LocusStructure": f"({motif})*",
            "VariantType": "Repeat",
        })

    batches = [variant_catalog[i:i+loci_per_run] for i in range(0, len(variant_catalog), loci_per_run)] or [[]]
    for batch_i, current_variant_catalog in enumerate(batches):
        with open(f"{output_path_prefix}.{batch_i+1:03d}_of_{len(batches):03d}.json", "wt") as f:
            json.dump(current_variant_catalog, f, indent=3)
    if len(batches) > 1:
        with open(f"{output_path_prefix}.001_of_001.json", "wt") as f:
            json.dump(variant_catalog, f, indent=3)
    print(f"Wrote {len(batches):,d} ExpansionHunter variant catalog shard(s) "
          f"({len(variant_catalog):,d} loci) to {output_path_prefix}*.json")


def write_gangstr_hipstr_or_longtr_repeat_specs(locus_set, output_path_prefix, tool, loci_per_run=None):
    """Write GangSTR/HipSTR/LongTR repeat-spec bed files (sharded into loci_per_run-sized batches)."""
    if tool not in ("gangstr", "hipstr", "longtr"):
        raise ValueError(f"Invalid tool arg: '{tool}'. Must be 'gangstr', 'hipstr', or 'longtr'")

    locus_list = list(sorted(locus_set))
    batches = [locus_list] if loci_per_run is None else [
        locus_list[i:i+loci_per_run] for i in range(0, len(locus_list), loci_per_run)]

    for batch_i, current_repeat_specs in enumerate(batches):
        with open(f"{output_path_prefix}.{batch_i+1:03d}_of_{len(batches):03d}.bed", "wt") as f:
            for chrom, start_0based, end_1based, motif in current_repeat_specs:
                if (end_1based - start_0based) / len(motif) <= 1:
                    # only 1 repeat in the reference -> HipSTR/LongTR error out on these
                    continue
                if tool == "gangstr":
                    output_fields = [chrom, start_0based + 1, end_1based, len(motif), motif]
                elif tool == "hipstr":
                    if len(motif) > 9:
                        continue  # HipSTR doesn't support motifs longer than 9bp
                    output_fields = [chrom, start_0based + 1, end_1based, len(motif),
                                     int((end_1based - start_0based)/len(motif)),
                                     f"{chrom}-{start_0based}-{end_1based}-{motif}"]
                else:  # longtr
                    output_fields = [chrom, start_0based + 1, end_1based, len(motif),
                                     int((end_1based - start_0based)/len(motif)),
                                     f"{chrom}-{start_0based}-{end_1based}-{motif}"]
                f.write("\t".join(map(str, output_fields)) + "\n")

    print(f"Wrote {len(batches):,d} {tool} repeat spec bed file(s) to {output_path_prefix}*.bed")


def write_trgt_catalog(locus_set, output_path):
    """Write a TRGT catalog bed (ID=..;MOTIFS=..;STRUC=(motif)n)."""
    with open(output_path, "wt") as f:
        for unmodified_chrom, start_0based, end_1based, motif in sorted(locus_set):
            chrom = unmodified_chrom.replace("chr", "")
            column4 = f"ID={chrom}-{start_0based}-{end_1based}-{motif};MOTIFS={motif};STRUC=({motif})n"
            f.write("\t".join(map(str, [unmodified_chrom, start_0based, end_1based, column4])) + "\n")
    print(f"Wrote {len(locus_set):,d} loci to {output_path}")


def write_positive_loci_bed(locus_set, output_path):
    """Write the plain positive_loci bed (chrom, start0, end1, motif, motif_len), then bgzip+tabix.

    This file is consumed directly by inquiSTR.
    """
    with open(output_path, "wt") as f:
        for chrom, start_0based, end_1based, motif in sorted(locus_set):
            f.write("\t".join(map(str, [chrom, start_0based, end_1based, motif, len(motif)])) + "\n")
    os.system(f"bgzip -f {output_path}")
    os.system(f"tabix -f {output_path}.gz")
    print(f"Wrote {len(locus_set):,d} loci to {output_path}.gz")


def main():
    args = parse_args()
    prefix = args.output_filename_prefix or os.path.basename(args.genotypes_tsv_path).split(".")[0]

    df = pd.read_table(args.genotypes_tsv_path, dtype={"Chrom": str, "Motif": str})
    print(f"Parsed {len(df):,d} rows from {args.genotypes_tsv_path}")

    positive_loci = generate_set_of_positive_loci(df)
    print(f"Generated {len(positive_loci):,d} positive (variant) loci")

    os.makedirs(args.output_dir, exist_ok=True)
    out = lambda name: os.path.join(args.output_dir, f"{prefix}.{name}")

    if not args.only or "eh" in args.only:
        write_expansion_hunter_variant_catalogs(positive_loci, out("EHv5"),
                                                loci_per_run=args.expansion_hunter_loci_per_run)
    if not args.only or "gangstr" in args.only:
        write_gangstr_hipstr_or_longtr_repeat_specs(positive_loci, out("GangSTR"), "gangstr",
                                                    loci_per_run=args.gangstr_loci_per_run)
    if not args.only or "hipstr" in args.only:
        write_gangstr_hipstr_or_longtr_repeat_specs(positive_loci, out("HipSTR"), "hipstr",
                                                    loci_per_run=args.gangstr_loci_per_run)
    if not args.only or "longtr" in args.only:
        write_gangstr_hipstr_or_longtr_repeat_specs(positive_loci, out("LongTR"), "longtr", loci_per_run=None)
    if not args.only or "trgt" in args.only:
        write_trgt_catalog(positive_loci, out("TRGT_repeat_catalog.bed"))

    # always write the plain loci bed (used by inquiSTR)
    write_positive_loci_bed(positive_loci, os.path.join(args.output_dir, f"{prefix}.bed"))

    print("Done")


if __name__ == "__main__":
    main()
