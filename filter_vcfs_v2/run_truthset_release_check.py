"""Check the HG002 truth set against the GIAB T2T-Q100 benchmark before publishing a release.

The truth set is scored with Aardvark against an independent benchmark, and the result is compared to a checked-in
baseline. Anything that drops a guarded metric by more than the tolerance fails the check.

GIAB is the truth and this project's own truth set is the query, so a false negative is a variant the truth set is
missing and a false positive is one it invented. Both are truth-set defects, and every tool ranking downstream
inherits them.

What is guarded, and why each metric is here:

    BASEPAIR ALL recall/precision/f1  The headline quality number, representation-invariant, so a caller that writes
                                      the same haplotype with different records is not penalized.
    GT ALL recall/precision/f1        The hap.py-like per-genotype score, which is what most other benchmarking
                                      tools report.
    GT and BASEPAIR per-chromosome    A per-chromosome recall guard. A genome-wide number is too blunt for a whole
    recall                            chromosome breaking: when dipcall's ".|1" chrX/chrY encoding gave both
                                      chromosomes a GT recall of exactly 0, genome-wide BASEPAIR F1 only moved from
                                      0.9932 to 0.9830. Per chromosome the same failure is unmissable.
    GT ALL zygosity_agreement         1 - truth_fn_gt / truth_total, the share of truth variants that are not pure
                                      zygosity errors. The chrX/chrY encoding drove this from 0.9994 to 0.9755, so
                                      it catches the same class of bug on an autosome, where no per-chromosome
                                      recall guard would fire hard enough to notice.

Every guarded metric is a rate in [0, 1] where higher is better, so one tolerance applies to all of them.

Requires bcftools, tabix, bedtools, curl and an aardvark binary. Aardvark v1.0.0 ships only a linux-x86_64 binary,
so on macOS it has to be built from source, and --aardvark then points at the built binary.

Usage:
    python3 run_truthset_release_check.py \
        --truth-vcf results/HG002.high_confidence_regions.vcf.gz \
        --catalog-bed results/combined.321_catalogs.merged.tandem_repeats.bed.gz \
        -R hg38.fa \
        --dipcall-confident-bed gs://str-truth-set-v2/dipcall_pipeline/HG002/HG002.dip.bed.gz
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

GIAB_BASE_URL = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/"
                 "NIST_HG002_DraftBenchmark_defrabbV0.020-20250117")
GIAB_TRUTH_VCF_FILENAME = "GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz"
GIAB_BENCHMARK_BED_FILENAME = "GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed"
GIAB_FILENAMES = (GIAB_TRUTH_VCF_FILENAME, GIAB_TRUTH_VCF_FILENAME + ".tbi", GIAB_BENCHMARK_BED_FILENAME)

# Aardvark drops any variant that is not fully contained in a region, so the catalog intervals are widened before
# being intersected with the two confident region sets. Without this, a variant that straddles a locus edge would be
# dropped from both the truth and the query, which is unbiased but shrinks coverage for no reason.
CATALOG_PADDING_BP = 50

# The clustering window that splits regions into independently solved sub-regions. PacBio's recommended setting for
# tandem repeats, well above the 50bp default, so that neighboring repeat variation is scored together.
MIN_VARIANT_GAP_BP = 1000

CHROMOSOMES = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])

DEFAULT_BASELINE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "truthset_release_check_baseline.tsv")

# A drop this large in a rate that normally sits above 0.98 is a real regression rather than run-to-run noise.
# Aardvark is deterministic, so the only reason a metric moves at all is that an input changed.
DEFAULT_TOLERANCE = 0.002


def run(command, **kwargs):
    """Run a shell command, raising if it fails.

    Args:
        command (str): the command to run.
        **kwargs: passed through to subprocess.run.

    Returns:
        subprocess.CompletedProcess: the finished process.
    """
    # Aardvark logs to stderr, which is unbuffered, so without this flush the commands and their output interleave
    # in the wrong order whenever this script's stdout is redirected to a file.
    print(f"> {command}", flush=True)
    return subprocess.run(command, shell=True, check=True, **kwargs)


def cat_command(path):
    """A shell command that writes the file to stdout, decompressing it if it is gzipped.

    Args:
        path (str): the file's path.

    Returns:
        str: the command.
    """
    if not path.endswith(".gz"):
        return f"cat {path}"

    # macOS ships gzcat rather than zcat, and its zcat only handles .Z files.
    return f"{'gzcat' if sys.platform == 'darwin' else 'zcat'} {path}"


def metric_sort_key(metric):
    """Order metrics so the genome-wide rows come first and chromosomes fall in their natural order.

    Args:
        metric (str): a metric key, like "BASEPAIR/chr10/recall".

    Returns:
        tuple: a sort key.
    """
    comparison, region, name = metric.split("/")
    if region == "ALL":
        return comparison, 0, 0, name
    return comparison, 1, CHROMOSOMES.index(region) if region in CHROMOSOMES else len(CHROMOSOMES), name


def guarded_metrics(summary_rows):
    """Pull the guarded metrics out of Aardvark's summary.tsv rows.

    Args:
        summary_rows (iterable): dicts, one per summary.tsv row, keyed by column name.

    Returns:
        dict: metric key (like "BASEPAIR/chrX/recall") to its value, a rate in [0, 1] where higher is better.
    """
    metrics = {}
    for row in summary_rows:
        if row["variant_type"] != "ALL" or row["comparison"] not in ("GT", "BASEPAIR"):
            continue

        if row["region_label"] == "ALL":
            for name in ("recall", "precision", "f1"):
                metrics[f"{row['comparison']}/ALL/{name}"] = float(row[f"metric_{name}"])

            # truth_fn_gt counts truth variants whose alleles matched but whose zygosity did not. Only the GT
            # comparison reports it; the BASEPAIR rows leave the column empty.
            if row["comparison"] == "GT" and int(row["truth_total"]) > 0:
                metrics["GT/ALL/zygosity_agreement"] = 1 - int(row["truth_fn_gt"]) / int(row["truth_total"])

        elif row["region_label"] in CHROMOSOMES:
            metrics[f"{row['comparison']}/{row['region_label']}/recall"] = float(row["metric_recall"])

    return metrics


def compare_to_baseline(observed, baseline, tolerance):
    """Compare observed metrics to the baseline.

    Args:
        observed (dict): metric key to observed value.
        baseline (dict): metric key to baseline value.
        tolerance (float): how far below its baseline a metric may fall before it counts as a regression.

    Returns:
        list: one dict per metric, with "metric", "baseline", "observed", "delta" and "status". A baseline metric
            that the run did not produce has status "missing" and no observed value; a metric with no baseline has
            status "new". Both "regressed" and "missing" are failures.
    """
    results = []
    for metric in sorted(set(baseline) | set(observed), key=metric_sort_key):
        if metric not in observed:
            results.append({"metric": metric, "baseline": baseline[metric], "observed": None, "delta": None,
                            "status": "missing"})
        elif metric not in baseline:
            results.append({"metric": metric, "baseline": None, "observed": observed[metric], "delta": None,
                            "status": "new"})
        else:
            delta = observed[metric] - baseline[metric]
            results.append({"metric": metric, "baseline": baseline[metric], "observed": observed[metric],
                            "delta": delta, "status": "regressed" if delta < -tolerance else "ok"})

    return results


def failures(results):
    """The results that fail the check.

    Args:
        results (list): dicts from compare_to_baseline.

    Returns:
        list: the subset whose status is a failure.
    """
    return [result for result in results if result["status"] in ("regressed", "missing")]


def format_result_line(result):
    """One aligned line of the results table.

    Args:
        result (dict): a result from compare_to_baseline.

    Returns:
        str: the formatted line, without a trailing newline.
    """
    baseline = "n/a" if result["baseline"] is None else f"{result['baseline']:.6f}"
    observed = "n/a" if result["observed"] is None else f"{result['observed']:.6f}"

    # Both values are rounded to 6 decimals in the baseline, so a delta smaller than that is a rounding artifact.
    # Printing it as +0.000000 rather than -0.000000 keeps an unchanged metric from looking like a small drop.
    delta = "" if result["delta"] is None else f"{result['delta'] if abs(result['delta']) >= 5e-7 else 0.0:+.6f}"
    marker = "  " if result["status"] == "ok" else "!!"
    return f"{marker} {result['metric']:38s} {baseline:>10s} {observed:>10s} {delta:>11s}  {result['status']}"


def format_results_table(results, verbose):
    """Render the results as a table.

    Args:
        results (list): dicts from compare_to_baseline.
        verbose (bool): whether to include the per-chromosome rows that passed. Failing rows are always included.

    Returns:
        str: the table.
    """
    lines = [f"   {'metric':38s} {'baseline':>10s} {'observed':>10s} {'delta':>11s}  status"]
    passing_chromosome_rows = 0
    for result in results:
        if not verbose and result["status"] == "ok" and result["metric"].split("/")[1] != "ALL":
            passing_chromosome_rows += 1
            continue
        lines.append(format_result_line(result))

    if passing_chromosome_rows:
        lines.append(f"   ({passing_chromosome_rows} per-chromosome metrics within tolerance, "
                     f"not shown; pass --verbose to list them)")

    return "\n".join(lines)


def parse_baseline_tsv(lines):
    """Read a baseline file.

    Args:
        lines (iterable): the file's lines. Blank lines and lines starting with "#" are ignored.

    Returns:
        dict: metric key to baseline value.
    """
    baseline = {}
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 2:
            raise ValueError(f"Expected 2 tab-delimited fields in a baseline line, found {len(fields)}: {line!r}")
        baseline[fields[0]] = float(fields[1])

    return baseline


def format_baseline_tsv(metrics, source_description):
    """Render metrics as a baseline file.

    The explanatory header is rewritten on every regeneration rather than preserved from the previous file, so a
    baseline written by --update-baseline always carries the same guidance as the one it replaces.

    Args:
        metrics (dict): metric key to value.
        source_description (str): what was scored to produce these numbers.

    Returns:
        str: the file's contents.
    """
    lines = [
        "# Baseline metrics for run_truthset_release_check.py.",
        "#",
        f"# Source: {source_description}",
        "#",
        "# Regenerate with --update-baseline, and only when the truth set legitimately changed. Say why in the",
        "# commit message, because these numbers are what every later release is judged against.",
        "# metric\tbaseline",
    ]
    lines += [f"{metric}\t{metrics[metric]:.6f}" for metric in sorted(metrics, key=metric_sort_key)]
    return "\n".join(lines) + "\n"


def read_summary_tsv(path):
    """Read Aardvark's summary.tsv.

    Args:
        path (str): the file's path.

    Returns:
        list: one dict per row, keyed by column name.
    """
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def download_giab_files(work_dir):
    """Download the GIAB T2T-Q100 benchmark files, skipping any that are already present.

    Args:
        work_dir (str): the check's working directory.

    Returns:
        str: the directory holding the downloaded files.
    """
    giab_dir = os.path.join(work_dir, "giab")
    os.makedirs(giab_dir, exist_ok=True)
    for filename in GIAB_FILENAMES:
        path = os.path.join(giab_dir, filename)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            run(f"curl -fsS -o {path} {GIAB_BASE_URL}/{filename}")

    return giab_dir


def prepare_query_vcf(truth_vcf_path, work_dir):
    """Uppercase the truth set's REF and ALT alleles so they can be compared to GIAB's.

    dipcall carries the reference's soft-masking through into REF and ALT, and the catalog lives in repeat regions,
    which is exactly where hg38 is soft-masked. GIAB writes uppercase, so the two disagree on case wherever that
    happens.

    Args:
        truth_vcf_path (str): the truth set VCF to check.
        work_dir (str): the check's working directory.

    Returns:
        tuple: (path to the uppercased VCF, number of records that had a lowercase base).
    """
    output_path = os.path.join(work_dir, "query.uppercased.vcf.gz")
    lowercase_count_path = os.path.join(work_dir, "lowercase_record_count.txt")
    run(f"bcftools view {truth_vcf_path} "
        f"| awk 'BEGIN{{OFS=\"\\t\"}} /^#/{{print; next}} "
        f"{{if ($4 ~ /[acgtn]/ || $5 ~ /[acgtn]/) n++; $4=toupper($4); $5=toupper($5); print}} "
        f"END{{print n+0 > \"{lowercase_count_path}\"}}' "
        f"| bgzip -@ 4 > {output_path}")
    run(f"tabix -f -p vcf {output_path}")

    with open(lowercase_count_path) as f:
        return output_path, int(f.read().strip())


def prepare_audit_regions(catalog_bed_path, giab_benchmark_bed_path, dipcall_confident_bed_path, work_dir):
    """Build the regions to score: catalog loci padded, inside both the GIAB and the dipcall confident regions.

    Args:
        catalog_bed_path (str): the merged tandem repeat catalog BED.
        giab_benchmark_bed_path (str): the GIAB benchmark BED.
        dipcall_confident_bed_path (str): the sample's dipcall confident regions BED, local or a gs:// path.
        work_dir (str): the check's working directory.

    Returns:
        str: the audit regions BED path.
    """
    regions_dir = os.path.join(work_dir, "regions")
    os.makedirs(regions_dir, exist_ok=True)

    if dipcall_confident_bed_path.startswith("gs://"):
        local_dipcall_bed_path = os.path.join(regions_dir, os.path.basename(dipcall_confident_bed_path))
        if not os.path.exists(local_dipcall_bed_path):
            run(f"gsutil -q cp {dipcall_confident_bed_path} {local_dipcall_bed_path}")
        dipcall_confident_bed_path = local_dipcall_bed_path

    keep_chromosomes = "grep -E '^chr([0-9]+|X|Y)\\b'"

    catalog_padded_path = os.path.join(regions_dir, "catalog_padded.bed")
    run(f"{cat_command(catalog_bed_path)} "
        f"| awk 'BEGIN{{OFS=\"\\t\"}} {{s=$2-{CATALOG_PADDING_BP}; if(s<0)s=0; print $1, s, $3+{CATALOG_PADDING_BP}}}' "
        f"| {keep_chromosomes} | sort -k1,1 -k2,2n | bedtools merge -i - > {catalog_padded_path}")

    giab_confident_path = os.path.join(regions_dir, "giab_confident.bed")
    run(f"{cat_command(giab_benchmark_bed_path)} | {keep_chromosomes} | sort -k1,1 -k2,2n | bedtools merge -i - "
        f"> {giab_confident_path}")

    dipcall_confident_path = os.path.join(regions_dir, "dipcall_confident.bed")
    run(f"{cat_command(dipcall_confident_bed_path)} | {keep_chromosomes} | sort -k1,1 -k2,2n "
        f"| bedtools merge -i - > {dipcall_confident_path}")

    audit_regions_path = os.path.join(regions_dir, "audit_regions.bed")
    run(f"bedtools intersect -a {catalog_padded_path} -b {giab_confident_path} "
        f"| bedtools intersect -a - -b {dipcall_confident_path} "
        f"| sort -k1,1 -k2,2n | bedtools merge -i - > {audit_regions_path}")

    return audit_regions_path


def write_stratification(audit_regions_path, work_dir):
    """Split the audit regions by chromosome so Aardvark reports a per-chromosome row.

    Aardvark takes hap.py's stratification format: a root TSV whose rows are a label and a BED path relative to it.

    Args:
        audit_regions_path (str): the audit regions BED.
        work_dir (str): the check's working directory.

    Returns:
        str: the root stratification TSV path.
    """
    stratification_dir = os.path.join(work_dir, "stratification")
    beds_dir = os.path.join(stratification_dir, "beds")
    os.makedirs(beds_dir, exist_ok=True)

    regions_by_chromosome = {}
    with open(audit_regions_path) as f:
        for line in f:
            regions_by_chromosome.setdefault(line.split("\t")[0], []).append(line)

    root_path = os.path.join(stratification_dir, "stratification.tsv")
    with open(root_path, "w") as root_file:
        for chromosome in CHROMOSOMES:
            if chromosome not in regions_by_chromosome:
                continue
            with open(os.path.join(beds_dir, f"{chromosome}.bed"), "w") as bed_file:
                bed_file.writelines(regions_by_chromosome[chromosome])
            root_file.write(f"{chromosome}\tbeds/{chromosome}.bed\n")

    return root_path


def run_aardvark(aardvark_path, reference_path, giab_dir, query_vcf_path, query_sample, audit_regions_path,
                 stratification_path, threads, output_dir):
    """Score the truth set against the GIAB benchmark.

    Args:
        aardvark_path (str): the aardvark binary.
        reference_path (str): the reference FASTA.
        giab_dir (str): the directory holding the downloaded GIAB files.
        query_vcf_path (str): the uppercased truth set VCF.
        query_sample (str): the sample name to read from the query VCF.
        audit_regions_path (str): the audit regions BED.
        stratification_path (str): the root stratification TSV.
        threads (int): how many threads to give Aardvark.
        output_dir (str): where Aardvark writes summary.tsv.

    Returns:
        str: the summary.tsv path.
    """
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    run(f"{aardvark_path} compare "
        f"--reference {reference_path} "
        f"--truth-vcf {os.path.join(giab_dir, GIAB_TRUTH_VCF_FILENAME)} --truth-sample HG002 "
        f"--query-vcf {query_vcf_path} --query-sample {query_sample} "
        f"--regions {audit_regions_path} "
        f"--stratifications {stratification_path} "
        f"--min-variant-gap {MIN_VARIANT_GAP_BP} "
        f"--threads {threads} "
        f"--compare-label truthset_release_check "
        f"--output-dir {output_dir}")

    return os.path.join(output_dir, "summary.tsv")


def haploid_normalization_hint(truth_vcf_path, results):
    """A diagnostic for the failure mode this check was written for, or None if it does not apply.

    Args:
        truth_vcf_path (str): the truth set VCF being checked.
        results (list): dicts from compare_to_baseline.

    Returns:
        str: a hint to print alongside the failures, or None.
    """
    regressed_sex_chromosomes = sorted({result["metric"].split("/")[1] for result in failures(results)
                                        if result["metric"].split("/")[1] in ("chrX", "chrY")})
    if not regressed_sex_chromosomes:
        return None

    header = subprocess.run(f"bcftools view -h {truth_vcf_path}", shell=True, check=True, capture_output=True,
                            text=True).stdout
    if "##normalizeHaploidGenotypes=" in header:
        return None

    return (f"{' and '.join(regressed_sex_chromosomes)} regressed and {os.path.basename(truth_vcf_path)} carries no "
            f"##normalizeHaploidGenotypes header, so this VCF was very likely built without the haploid genotype "
            f"normalization step. dipcall writes a male sample's chrX outside the PAR and all of its chrY as '.|1', "
            f"which no benchmarking tool can match. Check that run_filter_vcf_to_tandem_repeats.py ran "
            f"normalize_haploid_genotypes.py and that this sample's sex is set in the metadata table.")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-vcf", required=True, help="The truth set VCF to check.")
    p.add_argument("--catalog-bed", required=True, help="The merged tandem repeat catalog BED, used to pick the "
                                                        "regions to score.")
    p.add_argument("-R", "--reference", required=True, help="The hg38 reference FASTA.")
    p.add_argument("--dipcall-confident-bed", required=True,
                   help="The sample's dipcall confident regions BED. A gs:// path is copied locally first.")
    p.add_argument("--query-sample", default="syndip", help="The sample name to read from the truth set VCF.")
    p.add_argument("--aardvark", default="aardvark",
                   help="The aardvark binary. Upstream ships only a linux-x86_64 build, so on macOS point this at a "
                        "binary built from source.")
    p.add_argument("--work-dir", default="truthset_release_check",
                   help="Where to cache the GIAB download and write intermediate files.")
    p.add_argument("--baseline", default=DEFAULT_BASELINE_PATH, help="The baseline TSV to compare against.")
    p.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE,
                   help="How far below its baseline a metric may fall before it counts as a regression.")
    p.add_argument("--threads", type=int, default=os.cpu_count(), help="How many threads to give Aardvark.")
    p.add_argument("--update-baseline", action="store_true",
                   help="Write this run's metrics to the baseline file instead of checking them. Only do this when "
                        "the truth set legitimately changed, and say why in the commit message.")
    p.add_argument("--verbose", action="store_true", help="List every metric, not just the genome-wide ones and the "
                                                          "failures.")
    args = p.parse_args()

    if not shutil.which(args.aardvark) and not os.path.isfile(args.aardvark):
        p.error(f"aardvark not found at {args.aardvark!r}. Pass --aardvark, or see docs/install.md in "
                f"https://github.com/PacificBiosciences/aardvark. Building from source is required on macOS.")

    os.makedirs(args.work_dir, exist_ok=True)

    giab_dir = download_giab_files(args.work_dir)
    query_vcf_path, lowercase_record_count = prepare_query_vcf(args.truth_vcf, args.work_dir)
    audit_regions_path = prepare_audit_regions(args.catalog_bed, os.path.join(giab_dir, GIAB_BENCHMARK_BED_FILENAME),
                                               args.dipcall_confident_bed, args.work_dir)
    summary_tsv_path = run_aardvark(
        args.aardvark, args.reference, giab_dir, query_vcf_path, args.query_sample, audit_regions_path,
        write_stratification(audit_regions_path, args.work_dir), args.threads,
        os.path.join(args.work_dir, "aardvark_out"))

    observed = guarded_metrics(read_summary_tsv(summary_tsv_path))
    if not observed:
        raise ValueError(f"{summary_tsv_path} has no rows this check knows how to read")

    print()
    if lowercase_record_count:
        # Worth knowing but not worth failing on: it costs about 0.0004 of BASEPAIR F1, and this check uppercases
        # before scoring so it never reaches the metrics. It does reach anyone who benchmarks the published VCF.
        print(f"Note: {lowercase_record_count:,d} records carry soft-masked lowercase REF/ALT bases, which this "
              f"check uppercases before scoring. Any tool comparing the published VCF against an uppercase "
              f"benchmark sees them as mismatches.")
        print()

    if args.update_baseline:
        with open(args.baseline, "w") as f:
            # The basename rather than the path, so a baseline generated from a scratch copy of the VCF still reads
            # as a description of the release rather than of one machine's directory layout.
            f.write(format_baseline_tsv(observed, f"{os.path.basename(args.truth_vcf)} scored against GIAB "
                                                  f"{GIAB_TRUTH_VCF_FILENAME} over catalog loci padded "
                                                  f"{CATALOG_PADDING_BP}bp and intersected with both the GIAB and "
                                                  f"the dipcall confident regions."))
        print(f"Wrote {len(observed):,d} metrics to {args.baseline}")
        return

    with open(args.baseline) as f:
        results = compare_to_baseline(observed, parse_baseline_tsv(f), args.tolerance)

    print(format_results_table(results, args.verbose))
    print()

    if not failures(results):
        print(f"PASS: all {len(observed):,d} metrics are within {args.tolerance} of {args.baseline}")
        return

    print(f"FAIL: {len(failures(results)):,d} of {len(results):,d} metrics regressed by more than {args.tolerance}")
    hint = haploid_normalization_hint(args.truth_vcf, results)
    if hint:
        print()
        print(hint)

    sys.exit(1)


if __name__ == "__main__":
    main()
