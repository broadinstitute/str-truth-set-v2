"""Uppercase the REF and ALT alleles in a dipcall VCF.

dipcall copies alleles straight out of the reference, so it carries hg38's soft-masking through into REF and ALT. The
tandem repeat catalog lives in repeat regions, which is exactly where hg38 is soft-masked, so most of the records that
matter here come out lowercase. Measured on HG002, 2,783,855 of 4,772,787 records (58%) carry a lowercase base.

Benchmarks are written in uppercase, and the VCF specification treats the alleles as case-insensitive, so tools are
free to compare them either way. A tool that compares them as written scores every soft-masked record as a mismatch.
Aardvark loses about 0.0004 of basepair F1 on the HG002 truth set for this reason alone, and a stricter tool would
lose more.

This does not change what any step of this pipeline computes: str_analysis.filter_vcf_to_tandem_repeats uppercases
REF, ALT and the reference sequence as it reads them, so the catalog is identical either way. What changes is the
published high_confidence_regions VCF, which is what anyone benchmarking against this truth set actually reads.

An allele is only rewritten if it is a plain base sequence. Symbolic alleles ("<DEL>"), the missing allele ("*") and
breakend notation ("N[chr1:123[") are left exactly as they are, since uppercasing a breakend would corrupt the
contig name inside it. dipcall emits none of these, so this is a guard rather than a transformation that fires.

Reads a VCF on stdin, writes it to stdout, and changes nothing but the REF and ALT columns.

Usage:
    bedtools intersect ... | python3 uppercase_ref_and_alt.py | bgzip > out.vcf.gz
"""

import argparse
import sys

BASES = set("ACGTNacgtn")

PROVENANCE_HEADER_LINE = (
    '##uppercaseRefAndAlt=<Description="REF and ALT alleles were uppercased. dipcall carries the reference\'s '
    'soft-masking through into the alleles, and a tool that compares them as written scores every soft-masked '
    'record as a mismatch.">\n')


def uppercase_allele(allele):
    """Uppercase one allele, leaving anything that is not a plain base sequence alone.

    Args:
        allele (str): a single allele, e.g. "acgt", "<DEL>", "*" or "N[chr1:123[".

    Returns:
        str: the uppercased allele, or the allele unchanged if it is not a plain base sequence.
    """
    if allele and all(base in BASES for base in allele):
        return allele.upper()

    return allele


def uppercase_alt_field(alt):
    """Uppercase every allele in an ALT column, which may list several separated by commas.

    Args:
        alt (str): the record's ALT column.

    Returns:
        str: the ALT column with each plain base sequence uppercased.
    """
    return ",".join(uppercase_allele(allele) for allele in alt.split(","))


def uppercase_line(line):
    """Uppercase one VCF line's REF and ALT alleles.

    Args:
        line (str): one line of a VCF, including its trailing newline. Header lines are returned unchanged.

    Returns:
        tuple: (line, changed), where changed is True only if an allele was rewritten.
    """
    if line.startswith("#"):
        return line, False

    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        raise ValueError(f"Expected at least 8 tab-delimited fields in a VCF record, found {len(fields)}: {line!r}")

    ref, alt = uppercase_allele(fields[3]), uppercase_alt_field(fields[4])
    if (ref, alt) == (fields[3], fields[4]):
        return line, False

    fields[3], fields[4] = ref, alt

    return "\t".join(fields) + "\n", True


def uppercase_vcf(input_lines, output_file):
    """Write the VCF with its REF and ALT alleles uppercased, and count what changed.

    Args:
        input_lines (iterable): the input VCF's lines.
        output_file (file): where to write the output VCF.

    Returns:
        int: the number of records whose REF or ALT was rewritten.
    """
    changed_records = 0

    for line in input_lines:
        if line.startswith("#CHROM"):
            output_file.write(PROVENANCE_HEADER_LINE)

        output_line, changed = uppercase_line(line)
        output_file.write(output_line)
        changed_records += changed

    return changed_records


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()

    print(f"Uppercased the alleles of {uppercase_vcf(sys.stdin, sys.stdout):,d} records", file=sys.stderr)


if __name__ == "__main__":
    main()
