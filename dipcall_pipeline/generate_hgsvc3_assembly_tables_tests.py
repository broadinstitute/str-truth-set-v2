import unittest

from generate_hgsvc3_assembly_tables import choose_paternal_and_maternal_haplotypes


class ChoosePaternalAndMaternalHaplotypesTests(unittest.TestCase):

    def test_male_with_chry_in_hap1(self):
        percents = {"hap1": {"chrX": 0.3, "chrY": 79.3}, "hap2": {"chrX": 99.0, "chrY": 3.1}}
        self.assertEqual(choose_paternal_and_maternal_haplotypes("HG00096", "male", percents), ("hap1", "hap2"))

    def test_male_with_chry_in_hap2(self):
        percents = {"hap1": {"chrX": 99.8, "chrY": 0.0}, "hap2": {"chrX": 0.0, "chrY": 86.2}}
        self.assertEqual(choose_paternal_and_maternal_haplotypes("HG00731", "male", percents), ("hap2", "hap1"))

    def test_male_with_chrx_and_chry_in_both_haplotypes(self):
        # the HGSVC2 layout, where each haplotype also includes the unassigned chrX and chrY reads
        percents = {"hap1": {"chrX": 98.0, "chrY": 80.0}, "hap2": {"chrX": 98.0, "chrY": 80.0}}
        with self.assertRaisesRegex(ValueError, "aren't in separate haplotypes"):
            choose_paternal_and_maternal_haplotypes("HG00096", "male", percents)

    def test_male_with_too_little_chry(self):
        percents = {"hap1": {"chrX": 0.0, "chrY": 10.0}, "hap2": {"chrX": 99.0, "chrY": 0.0}}
        with self.assertRaisesRegex(ValueError, "aren't in separate haplotypes"):
            choose_paternal_and_maternal_haplotypes("HG00096", "male", percents)

    def test_female(self):
        percents = {"hap1": {"chrX": 99.5, "chrY": 0.0}, "hap2": {"chrX": 99.9, "chrY": 0.2}}
        self.assertEqual(choose_paternal_and_maternal_haplotypes("HG00514", "female", percents), ("hap1", "hap2"))

    def test_female_with_chry(self):
        percents = {"hap1": {"chrX": 0.0, "chrY": 85.0}, "hap2": {"chrX": 99.9, "chrY": 0.0}}
        with self.assertRaisesRegex(ValueError, "don't each have one chrX and no chrY"):
            choose_paternal_and_maternal_haplotypes("HG06807", "female", percents)

    def test_female_with_chrx_in_only_one_haplotype(self):
        percents = {"hap1": {"chrX": 0.0, "chrY": 0.0}, "hap2": {"chrX": 102.7, "chrY": 0.0}}
        with self.assertRaisesRegex(ValueError, "don't each have one chrX and no chrY"):
            choose_paternal_and_maternal_haplotypes("HG00732", "female", percents)
        self.assertEqual(choose_paternal_and_maternal_haplotypes(
            "HG00732", "female", percents, allow_chrx_in_only_one_haplotype=True), ("hap1", "hap2"))

    def test_female_allowed_one_chrx_but_has_none(self):
        percents = {"hap1": {"chrX": 0.0, "chrY": 0.0}, "hap2": {"chrX": 10.0, "chrY": 0.0}}
        with self.assertRaisesRegex(ValueError, "don't each have one chrX and no chrY"):
            choose_paternal_and_maternal_haplotypes("HG00732", "female", percents, allow_chrx_in_only_one_haplotype=True)

    def test_unexpected_sex(self):
        percents = {"hap1": {"chrX": 99.0, "chrY": 0.0}, "hap2": {"chrX": 99.0, "chrY": 0.0}}
        with self.assertRaisesRegex(ValueError, "unexpected sex"):
            choose_paternal_and_maternal_haplotypes("HG00514", None, percents)


if __name__ == "__main__":
    unittest.main()
