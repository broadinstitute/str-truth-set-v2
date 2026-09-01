"""Tests for normalize_haploid_genotypes.py"""

import io
import unittest

from normalize_haploid_genotypes import (
    PROVENANCE_HEADER_LINE, genotype, haploid_genotype, haplotypes, is_haploid_region, normalize_line, normalize_vcf)


def record(chrom="chrX", pos=3_000_000, ref="A", alt="T", format_keys="GT:AD", sample_values=".|1:0,1"):
    """A single-sample VCF data line built from the fields a test cares about."""
    return "\t".join([chrom, str(pos), ".", ref, alt, "30", ".", ".", format_keys, sample_values]) + "\n"


class IsHaploidRegionTests(unittest.TestCase):

    def test_male_chrx_outside_par_is_haploid(self):
        self.assertTrue(is_haploid_region("chrX", 3_000_000, "A", "male"))

    def test_male_chrx_inside_par1_is_diploid(self):
        self.assertFalse(is_haploid_region("chrX", 1_000_000, "A", "male"))

    def test_male_chrx_inside_par2_is_diploid(self):
        self.assertFalse(is_haploid_region("chrX", 155_800_000, "A", "male"))

    def test_male_chrx_between_the_two_pars_is_haploid(self):
        self.assertTrue(is_haploid_region("chrX", 100_000_000, "A", "male"))

    def test_male_chrx_record_ending_exactly_on_the_par1_boundary_is_diploid(self):
        # dipcall requires POS + len(REF) <= 2781479, so the last diploid record starts one base before that
        self.assertFalse(is_haploid_region("chrX", 2_781_478, "A", "male"))

    def test_male_chrx_record_crossing_the_par1_boundary_is_haploid(self):
        # a deletion anchored inside the PAR but reaching past its end is not fully contained, so dipcall calls it
        # haploid even though its POS is inside the PAR
        self.assertTrue(is_haploid_region("chrX", 2_781_470, "A" * 20, "male"))

    def test_male_chry_par_is_haploid(self):
        # dipcall's PAR file only lists chrX, so a male's whole chrY comes from a single haplotype
        self.assertTrue(is_haploid_region("chrY", 1_000_000, "A", "male"))

    def test_male_chry_outside_par_is_haploid(self):
        self.assertTrue(is_haploid_region("chrY", 20_000_000, "A", "male"))

    def test_male_autosome_is_diploid(self):
        self.assertFalse(is_haploid_region("chr1", 1_000_000, "A", "male"))

    def test_female_chrx_is_diploid(self):
        self.assertFalse(is_haploid_region("chrX", 3_000_000, "A", "female"))

    def test_female_chry_is_diploid(self):
        self.assertFalse(is_haploid_region("chrY", 3_000_000, "A", "female"))

    def test_contig_names_without_a_chr_prefix(self):
        self.assertTrue(is_haploid_region("X", 3_000_000, "A", "male"))
        self.assertTrue(is_haploid_region("Y", 3_000_000, "A", "male"))
        self.assertFalse(is_haploid_region("X", 1_000_000, "A", "male"))
        self.assertFalse(is_haploid_region("1", 3_000_000, "A", "male"))


class GenotypeParsingTests(unittest.TestCase):

    def test_haplotypes_splits_phased_and_unphased_genotypes(self):
        self.assertEqual(haplotypes(".|1"), [".", "1"])
        self.assertEqual(haplotypes("0/1"), ["0", "1"])
        self.assertEqual(haplotypes("1"), ["1"])

    def test_genotype_reads_the_gt_subfield_by_name(self):
        self.assertEqual(genotype(record(format_keys="AD:GT", sample_values="0,1:.|1").rstrip("\n").split("\t")), ".|1")

    def test_genotype_is_none_when_format_declares_no_gt(self):
        self.assertIsNone(genotype(record(format_keys="AD", sample_values="0,1").rstrip("\n").split("\t")))

    def test_missing_first_haplotype(self):
        self.assertEqual(haploid_genotype(".|1"), "1")

    def test_missing_second_haplotype(self):
        self.assertEqual(haploid_genotype("1|."), "1")

    def test_missing_haplotype_with_a_non_first_alt_allele(self):
        self.assertEqual(haploid_genotype(".|2"), "2")

    def test_missing_haplotype_carrying_the_reference_allele(self):
        self.assertEqual(haploid_genotype("0|."), "0")

    def test_unphased_genotype_with_a_missing_haplotype(self):
        self.assertEqual(haploid_genotype("./1"), "1")

    def test_both_haplotypes_missing_becomes_a_haploid_no_call(self):
        self.assertEqual(haploid_genotype(".|."), ".")

    def test_both_haplotypes_called_is_left_alone(self):
        self.assertIsNone(haploid_genotype("0|1"))
        self.assertIsNone(haploid_genotype("1|1"))
        self.assertIsNone(haploid_genotype("2|1"))

    def test_already_haploid_genotype_is_left_alone(self):
        self.assertIsNone(haploid_genotype("1"))
        self.assertIsNone(haploid_genotype("."))


class NormalizeLineTests(unittest.TestCase):

    def test_male_chrx_genotype_is_rewritten_and_the_other_fields_are_untouched(self):
        line, changed = normalize_line(record(chrom="chrX", sample_values=".|1:0,1"), "male")
        self.assertTrue(changed)
        self.assertEqual(line, record(chrom="chrX", sample_values="1:0,1"))

    def test_male_chry_genotype_is_rewritten(self):
        line, changed = normalize_line(record(chrom="chrY", sample_values="1|.:1,0"), "male")
        self.assertTrue(changed)
        self.assertEqual(line, record(chrom="chrY", sample_values="1:1,0"))

    def test_male_chrx_par_genotype_is_left_alone(self):
        line, changed = normalize_line(record(chrom="chrX", pos=1_000_000, sample_values="0|1:1,1"), "male")
        self.assertFalse(changed)
        self.assertEqual(line, record(chrom="chrX", pos=1_000_000, sample_values="0|1:1,1"))

    def test_male_chrx_record_with_both_haplotypes_called_is_left_alone(self):
        line, changed = normalize_line(record(chrom="chrX", sample_values="0|1:1,1"), "male")
        self.assertFalse(changed)
        self.assertEqual(line, record(chrom="chrX", sample_values="0|1:1,1"))

    def test_male_autosome_is_left_alone(self):
        line, changed = normalize_line(record(chrom="chr1", sample_values="0|1:1,1"), "male")
        self.assertFalse(changed)
        self.assertEqual(line, record(chrom="chr1", sample_values="0|1:1,1"))

    def test_female_chrx_is_left_alone(self):
        line, changed = normalize_line(record(chrom="chrX", sample_values=".|1:0,1"), "female")
        self.assertFalse(changed)
        self.assertEqual(line, record(chrom="chrX", sample_values=".|1:0,1"))

    def test_gt_is_rewritten_in_place_when_it_is_not_the_first_format_field(self):
        line, changed = normalize_line(record(format_keys="AD:GT", sample_values="0,1:.|1"), "male")
        self.assertTrue(changed)
        self.assertEqual(line, record(format_keys="AD:GT", sample_values="0,1:1"))

    def test_record_without_a_gt_is_left_alone(self):
        line, changed = normalize_line(record(format_keys="AD", sample_values="0,1"), "male")
        self.assertFalse(changed)
        self.assertEqual(line, record(format_keys="AD", sample_values="0,1"))

    def test_header_lines_are_left_alone(self):
        for header_line in "##fileformat=VCFv4.2\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsyndip\n":
            line, changed = normalize_line(header_line, "male")
            self.assertFalse(changed)
            self.assertEqual(line, header_line)

    def test_multi_sample_vcf_is_rejected(self):
        with self.assertRaises(ValueError):
            normalize_line(record().rstrip("\n") + "\t0|1:1,1\n", "male")


class NormalizeVcfTests(unittest.TestCase):

    HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsyndip\n"

    def normalize(self, lines, sex):
        output_file = io.StringIO()
        counts = normalize_vcf(io.StringIO("".join(lines)), output_file, sex)
        return output_file.getvalue(), counts

    def test_provenance_line_is_written_immediately_before_the_column_header(self):
        output, _ = self.normalize([self.HEADER], "male")
        self.assertEqual(
            output,
            "##fileformat=VCFv4.2\n" + PROVENANCE_HEADER_LINE.format(sex="male") + self.HEADER.splitlines(True)[1])

    def test_counts_rewritten_records(self):
        _, counts = self.normalize(
            [self.HEADER,
             record(chrom="chrX", sample_values=".|1:0,1"),
             record(chrom="chrY", sample_values="1|.:1,0"),
             record(chrom="chr1", sample_values="0|1:1,1")],
            "male")
        self.assertEqual(counts, {"rewritten": 2, "half_missing_not_rewritten": 0})

    def test_counts_a_male_chrx_par_record_that_keeps_a_missing_haplotype(self):
        _, counts = self.normalize(
            [self.HEADER, record(chrom="chrX", pos=1_000_000, sample_values="0|.:1,0")], "male")
        self.assertEqual(counts, {"rewritten": 0, "half_missing_not_rewritten": 1})

    def test_a_female_sample_with_a_missing_haplotype_is_counted_so_the_sex_label_can_be_checked(self):
        output, counts = self.normalize([self.HEADER, record(chrom="chrX", sample_values=".|1:0,1")], "female")
        self.assertEqual(counts, {"rewritten": 0, "half_missing_not_rewritten": 1})
        self.assertIn(".|1:0,1", output)

    def test_a_no_call_on_both_haplotypes_is_not_counted_as_half_missing(self):
        _, counts = self.normalize([self.HEADER, record(chrom="chr1", sample_values=".|.:0,0")], "female")
        self.assertEqual(counts, {"rewritten": 0, "half_missing_not_rewritten": 0})


if __name__ == "__main__":
    unittest.main()
