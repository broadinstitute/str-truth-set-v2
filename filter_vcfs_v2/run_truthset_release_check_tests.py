"""Tests for run_truthset_release_check.py."""

import unittest

from run_truthset_release_check import (cat_command, compare_to_baseline, failures, format_baseline_tsv,
                                        format_results_table, guarded_metrics, metric_sort_key, parse_baseline_tsv)


def summary_row(comparison, region_label, variant_type="ALL", recall="0.99", precision="0.99", f1="0.99",
                truth_total="1000", truth_fn_gt=""):
    """One Aardvark summary.tsv row, with only the columns this check reads.

    Args:
        comparison (str): the row's comparison type.
        region_label (str): the row's stratification label.
        variant_type (str): the row's variant type.
        recall (str): metric_recall, as Aardvark writes it.
        precision (str): metric_precision.
        f1 (str): metric_f1.
        truth_total (str): truth_total.
        truth_fn_gt (str): truth_fn_gt. Aardvark leaves this empty on non-GT rows.

    Returns:
        dict: the row.
    """
    return {"comparison": comparison, "region_label": region_label, "variant_type": variant_type,
            "metric_recall": recall, "metric_precision": precision, "metric_f1": f1, "truth_total": truth_total,
            "truth_fn_gt": truth_fn_gt}


class GuardedMetricsTests(unittest.TestCase):

    def test_genome_wide_rows_give_recall_precision_and_f1(self):
        metrics = guarded_metrics([summary_row("BASEPAIR", "ALL", recall="0.9939", precision="0.9925", f1="0.9932")])
        self.assertEqual(metrics, {"BASEPAIR/ALL/recall": 0.9939, "BASEPAIR/ALL/precision": 0.9925,
                                   "BASEPAIR/ALL/f1": 0.9932})

    def test_chromosome_rows_give_recall_only(self):
        metrics = guarded_metrics([summary_row("GT", "chrX", recall="0.9873")])
        self.assertEqual(metrics, {"GT/chrX/recall": 0.9873})

    def test_zygosity_agreement_comes_from_the_genome_wide_gt_row(self):
        metrics = guarded_metrics([summary_row("GT", "ALL", truth_total="1113764", truth_fn_gt="683")])
        self.assertAlmostEqual(metrics["GT/ALL/zygosity_agreement"], 1 - 683 / 1113764)

    def test_zygosity_agreement_falls_when_genotypes_are_half_missing(self):
        """The chrX/chrY encoding drove truth_fn_gt from 683 to 27,332, which this metric has to make visible."""
        fixed = guarded_metrics([summary_row("GT", "ALL", truth_total="1113764", truth_fn_gt="683")])
        broken = guarded_metrics([summary_row("GT", "ALL", truth_total="1113764", truth_fn_gt="27332")])
        self.assertGreater(fixed["GT/ALL/zygosity_agreement"] - broken["GT/ALL/zygosity_agreement"], 0.02)

    def test_zygosity_agreement_is_skipped_when_nothing_was_assessed(self):
        metrics = guarded_metrics([summary_row("GT", "ALL", truth_total="0", truth_fn_gt="0")])
        self.assertNotIn("GT/ALL/zygosity_agreement", metrics)

    def test_basepair_rows_have_no_zygosity_agreement(self):
        metrics = guarded_metrics([summary_row("BASEPAIR", "ALL", truth_fn_gt="")])
        self.assertNotIn("BASEPAIR/ALL/zygosity_agreement", metrics)

    def test_per_variant_type_rows_are_ignored(self):
        """Per-type rows are contaminated by cross-type representation swaps, so only the ALL rows are guarded."""
        self.assertEqual(guarded_metrics([summary_row("BASEPAIR", "ALL", variant_type="Snv")]), {})

    def test_other_comparison_types_are_ignored(self):
        self.assertEqual(guarded_metrics([summary_row("RECORD_BP", "ALL")]), {})

    def test_unrecognized_region_labels_are_ignored(self):
        self.assertEqual(guarded_metrics([summary_row("GT", "chrM")]), {})


class CatCommandTests(unittest.TestCase):

    def test_an_uncompressed_file_is_read_with_cat(self):
        self.assertEqual(cat_command("regions.bed"), "cat regions.bed")

    def test_a_gzipped_file_is_decompressed(self):
        self.assertRegex(cat_command("regions.bed.gz"), r"^(gzcat|zcat) regions\.bed\.gz$")


class CompareToBaselineTests(unittest.TestCase):

    def test_a_metric_at_its_baseline_passes(self):
        results = compare_to_baseline({"GT/ALL/f1": 0.98}, {"GT/ALL/f1": 0.98}, tolerance=0.002)
        self.assertEqual([result["status"] for result in results], ["ok"])

    def test_an_improvement_passes(self):
        results = compare_to_baseline({"GT/ALL/f1": 0.99}, {"GT/ALL/f1": 0.98}, tolerance=0.002)
        self.assertEqual(results[0]["status"], "ok")
        self.assertAlmostEqual(results[0]["delta"], 0.01)

    def test_a_drop_within_tolerance_passes(self):
        results = compare_to_baseline({"GT/ALL/f1": 0.9785}, {"GT/ALL/f1": 0.98}, tolerance=0.002)
        self.assertEqual(results[0]["status"], "ok")

    def test_a_drop_beyond_tolerance_regresses(self):
        results = compare_to_baseline({"GT/ALL/f1": 0.975}, {"GT/ALL/f1": 0.98}, tolerance=0.002)
        self.assertEqual(results[0]["status"], "regressed")

    def test_a_baseline_metric_the_run_did_not_produce_is_a_failure(self):
        results = compare_to_baseline({}, {"GT/chrY/recall": 0.93}, tolerance=0.002)
        self.assertEqual(results[0]["status"], "missing")
        self.assertIsNone(results[0]["observed"])
        self.assertEqual(failures(results), results)

    def test_a_metric_with_no_baseline_is_reported_but_does_not_fail(self):
        results = compare_to_baseline({"GT/chr1/recall": 0.98}, {}, tolerance=0.002)
        self.assertEqual(results[0]["status"], "new")
        self.assertEqual(failures(results), [])

    def test_the_sex_chromosome_encoding_regression_is_caught(self):
        """The failure this check exists for: GT recall on chrX and chrY goes to exactly zero."""
        baseline = {"GT/ALL/recall": 0.980049, "GT/chr1/recall": 0.980781, "GT/chrX/recall": 0.987256,
                    "GT/chrY/recall": 0.932407}
        observed = {"GT/ALL/recall": 0.956122, "GT/chr1/recall": 0.980781, "GT/chrX/recall": 0.0,
                    "GT/chrY/recall": 0.0}
        regressed = [result["metric"] for result in failures(compare_to_baseline(observed, baseline, 0.002))]
        self.assertEqual(sorted(regressed), ["GT/ALL/recall", "GT/chrX/recall", "GT/chrY/recall"])

    def test_results_are_ordered_with_genome_wide_metrics_first(self):
        observed = {"GT/chr2/recall": 0.98, "GT/ALL/f1": 0.98, "GT/chr10/recall": 0.98, "GT/chrX/recall": 0.98}
        results = compare_to_baseline(observed, observed, tolerance=0.002)
        self.assertEqual([result["metric"] for result in results],
                         ["GT/ALL/f1", "GT/chr2/recall", "GT/chr10/recall", "GT/chrX/recall"])


class MetricSortKeyTests(unittest.TestCase):

    def test_chromosomes_sort_numerically_rather_than_lexically(self):
        self.assertLess(metric_sort_key("GT/chr2/recall"), metric_sort_key("GT/chr10/recall"))

    def test_the_sex_chromosomes_sort_last(self):
        self.assertLess(metric_sort_key("GT/chr22/recall"), metric_sort_key("GT/chrX/recall"))
        self.assertLess(metric_sort_key("GT/chrX/recall"), metric_sort_key("GT/chrY/recall"))


class BaselineFileTests(unittest.TestCase):

    def test_comments_and_blank_lines_are_ignored(self):
        baseline = parse_baseline_tsv(["# a comment\n", "\n", "GT/ALL/f1\t0.979359\n"])
        self.assertEqual(baseline, {"GT/ALL/f1": 0.979359})

    def test_a_malformed_line_raises(self):
        with self.assertRaises(ValueError):
            parse_baseline_tsv(["GT/ALL/f1 0.979359\n"])

    def test_formatting_then_parsing_round_trips(self):
        metrics = {"GT/ALL/f1": 0.979359, "GT/chrX/recall": 0.987256, "BASEPAIR/ALL/f1": 0.993238}
        parsed = parse_baseline_tsv(format_baseline_tsv(metrics, "a comment").splitlines())
        self.assertEqual(parsed, metrics)

    def test_the_source_description_is_written_into_the_header(self):
        self.assertIn("# Source: some.vcf.gz vs GIAB", format_baseline_tsv({}, "some.vcf.gz vs GIAB"))

    def test_the_guidance_header_survives_regeneration(self):
        """--update-baseline rewrites the whole file, so the header has to come from the formatter, not the old file."""
        self.assertIn("only when the truth set legitimately changed", format_baseline_tsv({}, "some.vcf.gz"))


class FormatResultsTableTests(unittest.TestCase):

    def test_passing_chromosome_rows_are_summarized_rather_than_listed(self):
        observed = {"GT/ALL/f1": 0.98, "GT/chr1/recall": 0.98, "GT/chr2/recall": 0.98}
        table = format_results_table(compare_to_baseline(observed, observed, 0.002), verbose=False)
        self.assertIn("GT/ALL/f1", table)
        self.assertNotIn("GT/chr1/recall", table)
        self.assertIn("2 per-chromosome metrics within tolerance", table)

    def test_verbose_lists_every_row(self):
        observed = {"GT/ALL/f1": 0.98, "GT/chr1/recall": 0.98}
        table = format_results_table(compare_to_baseline(observed, observed, 0.002), verbose=True)
        self.assertIn("GT/chr1/recall", table)

    def test_a_failing_chromosome_row_is_always_listed(self):
        results = compare_to_baseline({"GT/chrX/recall": 0.0}, {"GT/chrX/recall": 0.987}, tolerance=0.002)
        table = format_results_table(results, verbose=False)
        self.assertIn("GT/chrX/recall", table)
        self.assertIn("regressed", table)

    def test_a_missing_metric_renders_without_an_observed_value(self):
        results = compare_to_baseline({}, {"GT/ALL/f1": 0.98}, tolerance=0.002)
        self.assertIn("n/a", format_results_table(results, verbose=False))

    def test_an_unchanged_metric_does_not_render_as_negative_zero(self):
        results = compare_to_baseline({"GT/ALL/f1": 0.9793594817}, {"GT/ALL/f1": 0.979359}, tolerance=0.002)
        table = format_results_table(results, verbose=False)
        self.assertIn("+0.000000", table)
        self.assertNotIn("-0.000000", table)


if __name__ == "__main__":
    unittest.main()
