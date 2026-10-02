"""Work out which DipCall batch produced each sample's high-confidence VCF in gs://str-truth-set-v2/filter_vcf_v2/.

gs://str-truth-set-v2/dipcall_pipeline/ keeps each DipCall run in a batch directory: the top level (the original
HPRC release 1, HGSVC and other assemblies), HPRC_release2/ and human579_assemblies/. A sample assembled in more
than one batch has a DipCall run in each. filter_vcf_v2/ was written flat, one directory per sample, so its paths do
not say which run a sample's truth came from.

The evidence used here: filter_vcf_v2's <sample>.high_confidence_regions.vcf.gz holds only variants inside the
high-confidence regions of the DipCall run it was made from (run_filter_vcf_to_tandem_repeats.py intersects the
run's VCF with that run's .dip.bed.gz). So every variant in it lies inside the source run's BED, while a different
run's BED, made from a different assembly, normally leaves some of them outside. For each sample, every candidate
batch's BED is tested against an evenly spaced subset of the VCF's variants, and the batch that contains all of
them is the source. all_assemblies_v2.321_samples.tsv, the metadata table the filter step was run with, is reported
alongside as a cross-check rather than used as the answer.

Writes sample_to_dipcall_batch.tsv with one row per sample directory in filter_vcf_v2/ (directories holding no
high-confidence VCF are not samples and are left out):
    sample_id                     sample directory name
    dipcall_batch                 "top_level", "HPRC_release2" or "human579_assemblies"; empty when no candidate
                                  batch contains every tested variant, or more than one does and the metadata
                                  table cannot break the tie
    evidence                      how dipcall_batch was decided
    metadata_table_batch          the batch all_assemblies_v2.321_samples.tsv lists, or empty if not listed
    candidate_batches             every batch with a DipCall BED for this sample
    n_variants_tested             variants tested against each candidate BED
    n_outside_by_batch            per candidate batch, how many of them fall outside its BED

Usage:
    python3 build_sample_to_dipcall_batch_table.py [--threads 8] [--every-nth-variant 50]
"""
import argparse
import bisect
import gzip
import io
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

FILTER_VCF_V2_DIR = "gs://str-truth-set-v2/filter_vcf_v2"
DIPCALL_DIR = "gs://str-truth-set-v2/dipcall_pipeline"
# Batch name -> its subdirectory under DIPCALL_DIR ("" for the top level).
DIPCALL_BATCH_SUBDIRECTORIES = {"top_level": "", "HPRC_release2": "HPRC_release2", "human579_assemblies": "human579_assemblies"}
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
METADATA_TSV_PATH = os.path.join(SCRIPT_DIR, "all_assemblies_v2.321_samples.tsv")
OUTPUT_PATH = os.path.join(SCRIPT_DIR, "sample_to_dipcall_batch.tsv")


def gsutil_ls(pattern):
    """Return the gs:// paths matching pattern, or [] when nothing matches."""
    result = subprocess.run(["gsutil", "ls", pattern], capture_output=True, text=True)
    if result.returncode != 0 and "matched no objects" not in result.stderr:
        raise RuntimeError(f"gsutil ls {pattern} failed: {result.stderr[-1000:]}")
    return result.stdout.split()


def read_gcs_bytes(path):
    return subprocess.run(["gsutil", "cat", path], capture_output=True, check=True).stdout


def load_bed_intervals(bed_bytes):
    """Return {chrom: (sorted starts, ends)} for a BED, merging nothing (DipCall's regions are disjoint)."""
    intervals = {}
    with gzip.open(io.BytesIO(bed_bytes), "rt") as f:
        for line in f:
            chrom, start, end = line.split("\t")[:3]
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    return {chrom: ([s for s, _ in sorted(v)], [e for _, e in sorted(v)]) for chrom, v in intervals.items()}


def count_outside(intervals, positions):
    """Count (chrom, 0-based position) pairs not inside any interval."""
    n_outside = 0
    for chrom, position in positions:
        if chrom not in intervals:
            n_outside += 1
            continue
        starts, ends = intervals[chrom]
        index = bisect.bisect_right(starts, position) - 1
        if index < 0 or position >= ends[index]:
            n_outside += 1
    return n_outside


def sampled_variant_positions(vcf_path, every_nth_variant):
    """Return every nth variant's (chrom, 0-based position) from a gzipped VCF on GCS, streamed."""
    command = (f"gsutil cat {vcf_path} | gzip -dc | "
               f"awk -v n={every_nth_variant} '!/^#/ {{ i++; if (i % n == 0) print $1\"\\t\"$2 }}'")
    output = subprocess.run(["bash", "-c", "set -o pipefail; " + command], capture_output=True, text=True, check=True).stdout
    return [(chrom, int(position) - 1) for chrom, position in (line.split("\t") for line in output.splitlines())]


def candidate_bed_paths(sample_id):
    """Return {batch: gs:// path of the sample's .dip.bed.gz} for every batch that has one."""
    paths = {}
    for batch, subdirectory in DIPCALL_BATCH_SUBDIRECTORIES.items():
        path = "/".join(p for p in (DIPCALL_DIR, subdirectory, sample_id, f"{sample_id}.dip.bed.gz") if p)
        if gsutil_ls(path):
            paths[batch] = path
    return paths


def decide_batch(sample_id, every_nth_variant, metadata_batch_by_sample):
    vcf_path = f"{FILTER_VCF_V2_DIR}/{sample_id}/{sample_id}.high_confidence_regions.vcf.gz"
    row = {"sample_id": sample_id, "dipcall_batch": "", "evidence": "",
           "metadata_table_batch": metadata_batch_by_sample.get(sample_id, ""),
           "candidate_batches": "", "n_variants_tested": 0, "n_outside_by_batch": ""}
    if not gsutil_ls(vcf_path):
        row["evidence"] = "no high_confidence_regions.vcf.gz in this directory"
        return row
    bed_path_by_batch = candidate_bed_paths(sample_id)
    row["candidate_batches"] = ",".join(bed_path_by_batch)
    if not bed_path_by_batch:
        row["evidence"] = "no DipCall BED for this sample in any batch"
        return row
    positions = sampled_variant_positions(vcf_path, every_nth_variant)
    row["n_variants_tested"] = len(positions)
    n_outside_by_batch = {batch: count_outside(load_bed_intervals(read_gcs_bytes(path)), positions)
                          for batch, path in bed_path_by_batch.items()}
    row["n_outside_by_batch"] = ",".join(f"{batch}:{n}" for batch, n in n_outside_by_batch.items())
    containing_batches = [batch for batch, n in n_outside_by_batch.items() if n == 0]
    if len(containing_batches) == 1:
        row["dipcall_batch"] = containing_batches[0]
        row["evidence"] = ("only candidate batch; all tested variants inside its BED" if len(bed_path_by_batch) == 1
                           else "the only candidate batch whose BED contains all tested variants")
    elif len(containing_batches) > 1 and row["metadata_table_batch"] in containing_batches:
        row["dipcall_batch"] = row["metadata_table_batch"]
        row["evidence"] = "several candidate BEDs contain all tested variants; metadata table breaks the tie"
    elif containing_batches:
        row["evidence"] = "several candidate BEDs contain all tested variants and the metadata table does not decide"
    else:
        row["evidence"] = "no candidate BED contains all tested variants"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--every-nth-variant", type=int, default=50)
    parser.add_argument("-s", "--sample-id", action="append", help="Only these samples (for testing).")
    args = parser.parse_args()

    metadata = pd.read_table(METADATA_TSV_PATH)
    metadata_batch_by_sample = dict(zip(metadata.sample_id, metadata.subdirectory.fillna("top_level")))

    sample_ids = args.sample_id or sorted(
        path.rstrip("/").split("/")[-1] for path in gsutil_ls(f"{FILTER_VCF_V2_DIR}/") if path.endswith("/"))
    print(f"{len(sample_ids)} sample directories in {FILTER_VCF_V2_DIR}")
    with ThreadPoolExecutor(args.threads) as pool:
        rows = list(pool.map(lambda s: decide_batch(s, args.every_nth_variant, metadata_batch_by_sample), sample_ids))

    df = pd.DataFrame(rows).sort_values("sample_id")
    # A directory with no high-confidence VCF is not a sample (e.g. backup__2026_09_03__combined_1_catalogs/).
    not_samples = df[df.n_variants_tested == 0]
    if len(not_samples):
        print(f"Leaving out {len(not_samples)} directories that hold no high-confidence VCF: "
              f"{', '.join(not_samples.sample_id)}")
    df = df[df.n_variants_tested > 0]
    df.to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {OUTPUT_PATH}")
    print(df.dipcall_batch.replace("", "(undecided)").value_counts().to_string())
    print(df.evidence.value_counts().to_string())
    disagree = df[(df.dipcall_batch != "") & (df.metadata_table_batch != "") & (df.dipcall_batch != df.metadata_table_batch)]
    print(f"{len(disagree)} samples where the evidence disagrees with the metadata table")


if __name__ == "__main__":
    main()
