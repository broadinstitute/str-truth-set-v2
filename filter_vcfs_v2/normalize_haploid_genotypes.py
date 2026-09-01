"""Rewrite dipcall's half-missing chrX/chrY genotypes as haploid genotypes.

dipcall fills in both haplotype slots on every record, so the regions it genotyped from a single assembly haplotype
come out as ".|1" (male chrX outside the PAR) or "1|." (male chrY). Benchmarking tools cannot match a "." haplotype,
so they score every one of those records as a false negative. Measured with Aardvark v1.0.0 against GIAB T2T-Q100
v1.1, the HG002 truth set gets recall of exactly 0 on both chromosomes, which accounts for 55% of all of its apparent
false negatives genome-wide. GIAB moved its own benchmark off this encoding for the same reason. This script rewrites
those genotypes to the haploid form ("1") that benchmarking tools expect.

Which records are haploid follows dipcall's own construction of the confident regions (run-dipcall v0.3), which for a
male sample takes chrX outside the PAR from hap2 only, all of chrY from hap1 only, and chrX inside the PAR from the
intersection of both haplotypes. chrX PAR records are therefore left alone, as is any record where both haplotypes
were called, so a genotype is never collapsed where dipcall actually called two alleles.

Reads a VCF on stdin, writes it to stdout, and changes nothing but the GT subfield of the records described above.

Usage:
    bedtools intersect ... | python3 normalize_haploid_genotypes.py --sex male | bgzip > out.vcf.gz
"""

import argparse
import sys

# hg38 chrX pseudoautosomal regions, copied from dipcall's data/hs38.PAR.bed, so each pair is a BED start and end.
# dipcall calls a record pseudoautosomal when start <= POS and POS + len(REF) <= end, and reproducing that comparison
# on these exact bounds is what keeps this script's notion of a haploid record identical to dipcall's.
CHRX_PAR_INTERVALS = ((0, 2781479), (155701383, 156030895))

PROVENANCE_HEADER_LINE = (
    '##normalizeHaploidGenotypes=<sex={sex},Description="chrX outside the PAR and all of chrY were genotyped from a '
    'single assembly haplotype, so their genotypes are written as haploid (1) rather than as a diploid genotype with '
    'one missing haplotype (.|1)">\n')


def is_haploid_region(chrom, pos_1based, ref, sex):
    """Whether dipcall genotyped this record from a single assembly haplotype.

    Args:
        chrom (str): the record's CHROM, with or without a "chr" prefix.
        pos_1based (int): the record's POS.
        ref (str): the record's REF allele.
        sex (str): "male" or "female".

    Returns:
        bool: True if the record falls in a region dipcall treats as haploid.
    """
    if sex != "male":
        return False

    chrom = chrom[3:] if chrom.startswith("chr") else chrom
    if chrom == "Y":
        # dipcall subtracts the PAR from chrX but not from chrY, so a male's whole chrY is hap1-only.
        return True
    if chrom != "X":
        return False

    return not any(start <= pos_1based and pos_1based + len(ref) <= end for start, end in CHRX_PAR_INTERVALS)


def genotype(fields):
    """The record's GT value.

    Args:
        fields (list): the record's tab-delimited fields.

    Returns:
        str: the GT value, or None if the record's FORMAT column doesn't declare a GT.
    """
    format_keys = fields[8].split(":")
    if "GT" not in format_keys:
        return None

    return fields[9].split(":")[format_keys.index("GT")]


def haplotypes(gt):
    """The genotype's per-haplotype allele indices.

    Args:
        gt (str): the record's GT, e.g. ".|1", "0/1" or "1".

    Returns:
        list: one string per haplotype, e.g. [".", "1"]. An uncalled haplotype is ".".
    """
    return gt.replace("/", "|").split("|")


def haploid_genotype(gt):
    """The haploid form of a two-haplotype genotype that has at most one called haplotype.

    Args:
        gt (str): the record's GT, e.g. ".|1", "1|.", "0|1" or "1".

    Returns:
        str: the haploid genotype, or None if the genotype should be left as it is because it is already haploid or
            because both of its haplotypes were called.
    """
    if len(haplotypes(gt)) != 2:
        return None

    called = [haplotype for haplotype in haplotypes(gt) if haplotype != "."]
    if len(called) == 1:
        return called[0]
    if not called:
        return "."
    return None


def normalize_line(line, sex):
    """Rewrite one VCF line's genotype as haploid if dipcall genotyped that record from a single haplotype.

    Args:
        line (str): one line of a VCF, including its trailing newline. Header lines are returned unchanged.
        sex (str): "male" or "female". A female sample has no haploid regions, so every line is returned unchanged.

    Returns:
        tuple: (line, changed), where changed is True only if the genotype was rewritten.
    """
    if line.startswith("#"):
        return line, False

    fields = line.rstrip("\n").split("\t")
    if len(fields) != 10:
        raise ValueError(f"Expected 10 tab-delimited fields in a single-sample VCF, found {len(fields)}: {line!r}")

    if not is_haploid_region(fields[0], int(fields[1]), fields[3], sex) or genotype(fields) is None:
        return line, False

    new_gt = haploid_genotype(genotype(fields))
    if new_gt is None:
        return line, False

    format_keys = fields[8].split(":")
    sample_values = fields[9].split(":")
    sample_values[format_keys.index("GT")] = new_gt
    fields[9] = ":".join(sample_values)

    return "\t".join(fields) + "\n", True


def normalize_vcf(input_lines, output_file, sex):
    """Write the VCF with the genotypes of dipcall's haploid regions rewritten as haploid, and count what changed.

    Args:
        input_lines (iterable): the input VCF's lines.
        output_file (file): where to write the output VCF.
        sex (str): "male" or "female".

    Returns:
        dict: counts of records whose genotype was "rewritten", and of records that have exactly one missing haplotype
            and were left alone ("half_missing_not_rewritten").
    """
    counts = {"rewritten": 0, "half_missing_not_rewritten": 0}

    for line in input_lines:
        if line.startswith("#CHROM"):
            output_file.write(PROVENANCE_HEADER_LINE.format(sex=sex))
        if line.startswith("#"):
            output_file.write(line)
            continue

        output_line, changed = normalize_line(line, sex)
        output_file.write(output_line)

        if changed:
            counts["rewritten"] += 1
            continue

        gt = genotype(line.rstrip("\n").split("\t"))
        if gt is not None and haplotypes(gt).count(".") == 1:
            counts["half_missing_not_rewritten"] += 1

    return counts


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sex", choices=["male", "female"], required=True,
                   help="The sample's sex. Only a male sample has haploid chrX/chrY regions, so a female sample's VCF "
                        "is passed through unchanged and only checked.")
    args = p.parse_args()

    counts = normalize_vcf(sys.stdin, sys.stdout, args.sex)

    print(f"Rewrote {counts['rewritten']:,d} genotypes as haploid", file=sys.stderr)
    print(f"Left {counts['half_missing_not_rewritten']:,d} records with a missing haplotype alone", file=sys.stderr)

    # A female sample has no haploid regions, so dipcall had both assembly haplotypes everywhere it called anything.
    # A missing haplotype means this sample's recorded sex disagrees with the VCF dipcall produced from it, and
    # continuing would write a truth set that still gets 0 recall on chrX/chrY in every benchmarking tool.
    if args.sex == "female" and counts["half_missing_not_rewritten"] > 0:
        raise ValueError(
            f"{counts['half_missing_not_rewritten']:,d} records have a missing haplotype, which dipcall only writes "
            f"where it genotyped from a single assembly haplotype. That should not happen for a sample recorded as "
            f"female. Check this sample's sex in the metadata table.")


if __name__ == "__main__":
    main()
