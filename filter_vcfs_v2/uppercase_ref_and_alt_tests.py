"""Tests for uppercase_ref_and_alt.py"""

import io
import unittest

from uppercase_ref_and_alt import (
    PROVENANCE_HEADER_LINE, uppercase_allele, uppercase_alt_field, uppercase_line, uppercase_vcf)


def record(ref="a", alt="T", chrom="chr1", pos=1_000_000):
    """A single-sample VCF data line built from the fields a test cares about."""
    return "\t".join([chrom, str(pos), ".", ref, alt, "30", ".", ".", "GT:AD", "1|1:0,2"]) + "\n"


class UppercaseAlleleTests(unittest.TestCase):

    def test_a_lowercase_allele_is_uppercased(self):
        self.assertEqual(uppercase_allele("acgt"), "ACGT")

    def test_an_uppercase_allele_is_unchanged(self):
        self.assertEqual(uppercase_allele("ACGT"), "ACGT")

    def test_a_mixed_case_allele_is_uppercased(self):
        self.assertEqual(uppercase_allele("aCgT"), "ACGT")

    def test_n_is_uppercased_like_any_other_base(self):
        self.assertEqual(uppercase_allele("acgtn"), "ACGTN")

    def test_the_missing_allele_is_left_alone(self):
        self.assertEqual(uppercase_allele("*"), "*")

    def test_a_symbolic_allele_is_left_alone(self):
        self.assertEqual(uppercase_allele("<DEL>"), "<DEL>")

    def test_a_breakend_is_left_alone(self):
        """Uppercasing a breakend would turn the contig name inside it into CHR1 and corrupt the record."""
        self.assertEqual(uppercase_allele("N[chr1:123["), "N[chr1:123[")

    def test_the_empty_allele_is_left_alone(self):
        self.assertEqual(uppercase_allele(""), "")


class UppercaseAltFieldTests(unittest.TestCase):

    def test_a_single_allele_is_uppercased(self):
        self.assertEqual(uppercase_alt_field("acgt"), "ACGT")

    def test_every_allele_in_a_multiallelic_record_is_uppercased(self):
        self.assertEqual(uppercase_alt_field("acgt,tt"), "ACGT,TT")

    def test_a_star_alongside_a_base_allele_keeps_the_star(self):
        self.assertEqual(uppercase_alt_field("acgt,*"), "ACGT,*")

    def test_a_symbolic_allele_alongside_a_base_allele_is_left_alone(self):
        self.assertEqual(uppercase_alt_field("acgt,<DEL>"), "ACGT,<DEL>")


class UppercaseLineTests(unittest.TestCase):

    def test_a_header_line_is_unchanged(self):
        self.assertEqual(uppercase_line("##contig=<ID=chr1>\n"), ("##contig=<ID=chr1>\n", False))

    def test_a_lowercase_ref_is_uppercased(self):
        line, changed = uppercase_line(record(ref="a", alt="T"))
        self.assertEqual(line.split("\t")[3], "A")
        self.assertTrue(changed)

    def test_a_lowercase_alt_is_uppercased(self):
        line, changed = uppercase_line(record(ref="A", alt="caa"))
        self.assertEqual(line.split("\t")[4], "CAA")
        self.assertTrue(changed)

    def test_an_already_uppercase_record_is_reported_as_unchanged(self):
        line, changed = uppercase_line(record(ref="A", alt="T"))
        self.assertEqual(line, record(ref="A", alt="T"))
        self.assertFalse(changed)

    def test_no_other_column_is_touched(self):
        original = record(ref="a", alt="c")
        line, _ = uppercase_line(original)
        self.assertEqual(line.split("\t")[:3], original.split("\t")[:3])
        self.assertEqual(line.split("\t")[5:], original.split("\t")[5:])

    def test_a_truncated_record_raises(self):
        with self.assertRaises(ValueError):
            uppercase_line("chr1\t100\t.\ta\tT\n")

    def test_a_sites_only_record_with_eight_columns_is_accepted(self):
        line, changed = uppercase_line("chr1\t100\t.\ta\tT\t30\t.\t.\n")
        self.assertEqual(line.split("\t")[3], "A")
        self.assertTrue(changed)


class UppercaseVcfTests(unittest.TestCase):

    def test_the_provenance_header_is_written_before_the_chrom_line(self):
        output = io.StringIO()
        uppercase_vcf(["##fileformat=VCFv4.2\n", "#CHROM\tPOS\n"], output)
        self.assertEqual(output.getvalue(),
                         "##fileformat=VCFv4.2\n" + PROVENANCE_HEADER_LINE + "#CHROM\tPOS\n")

    def test_changed_records_are_counted(self):
        output = io.StringIO()
        changed = uppercase_vcf([record(ref="a"), record(ref="A"), record(ref="c", alt="tt")], output)
        self.assertEqual(changed, 2)

    def test_every_record_is_written_whether_or_not_it_changed(self):
        output = io.StringIO()
        uppercase_vcf([record(ref="a"), record(ref="A")], output)
        self.assertEqual(len(output.getvalue().splitlines()), 2)

    def test_running_it_twice_changes_nothing_the_second_time(self):
        """The pipeline may be re-run over its own output, so the transform has to be idempotent."""
        first = io.StringIO()
        uppercase_vcf([record(ref="a", alt="cc")], first)
        second = io.StringIO()
        self.assertEqual(uppercase_vcf(first.getvalue().splitlines(keepends=True), second), 0)


if __name__ == "__main__":
    unittest.main()
