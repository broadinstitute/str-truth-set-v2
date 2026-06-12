"""One-off migration: move HG002 long-read TRGT tool_results into the TRGTv5 path.

The single pre-split "TRGT" tool used the TRGT v5.0.0 docker image (commit b22c413
renamed run_genotyping_tools.py's create_trgt_step default DOCKER_IMAGE to
TRGT_V5_DOCKER_IMAGE and split the tool into TRGTv3/TRGTv5). The VCF headers in the
old TRGT dirs confirm ##trgtVersion=5.0.0-3cf8401, identical to the new TRGTv5 run.
So the orphaned TRGT plots/data belong under TRGTv5.

For every object under tool_results/HG002/{pacbio,ONT,pacbio_isoseq}/TRGT/ this
renames the "TRGT" token to "TRGTv5" in both the directory segment and the filename
(e.g. .../TRGT/37x_coverage/HG002.TRGT.vcf.gz -> .../TRGTv5/37x_coverage/HG002.TRGTv5.vcf.gz
and tool_accuracy_....TRGT.svg -> tool_accuracy_....TRGTv5.svg), via a server-side
copy followed by deletion of the source (i.e. a move). Existing partial TRGTv5/37x
objects are overwritten with the equivalent same-version data.

Usage:
    python3 migrate_TRGT_to_TRGTv5.py            # dry run: print sample mapping + counts
    python3 migrate_TRGT_to_TRGTv5.py --execute  # perform the move
"""

import argparse
import concurrent.futures
import time
from google.api_core import exceptions as gax
from google.cloud import storage

BUCKET = "str-truth-set-v2"
PREFIXES = [
    "tool_results/HG002/pacbio/TRGT/",
    "tool_results/HG002/ONT/TRGT/",
    "tool_results/HG002/pacbio_isoseq/TRGT/",
]

# Transient GCS/network errors worth retrying (the first run died on a ReadTimeout).
TRANSIENT = (gax.ServiceUnavailable, gax.GatewayTimeout, gax.TooManyRequests,
             gax.InternalServerError, ConnectionError, TimeoutError)


def dst_name(src_name):
    """Map a source blob name to its TRGTv5 destination (rename the TRGT token)."""
    return src_name.replace("TRGT", "TRGTv5")


def move_blob(bucket, src_blob, attempts=6):
    """Server-side copy src -> TRGTv5 path, then delete src. Idempotent and resumable:
    copy overwrites an existing dest, and an already-deleted source is treated as done."""
    for attempt in range(attempts):
        try:
            bucket.copy_blob(src_blob, bucket, new_name=dst_name(src_blob.name), timeout=300)
            try:
                src_blob.delete(timeout=120)
            except gax.NotFound:
                pass  # already deleted on a prior (timed-out) attempt
            return
        except TRANSIENT as e:
            if attempt == attempts - 1:
                raise
            time.sleep(min(2 ** attempt, 30))  # exponential backoff, capped
        except Exception as e:
            # for ReadTimeout (requests) and similar non-gax timeouts, retry too
            if "timeout" in str(e).lower() or "timed out" in str(e).lower():
                if attempt == attempts - 1:
                    raise
                time.sleep(min(2 ** attempt, 30))
            else:
                raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Perform the move (default is a dry run).")
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--passes", type=int, default=5, help="Re-list and retry remaining sources up to this many times.")
    args = parser.parse_args()

    bucket = storage.Client().bucket(BUCKET)

    if not args.execute:
        blobs = [b for prefix in PREFIXES for b in bucket.list_blobs(prefix=prefix)]
        print(f"Found {len(blobs)} objects under {len(PREFIXES)} TRGT prefixes")
        for src_blob in blobs[:5]:
            print(f"  {src_blob.name}\n    -> {dst_name(src_blob.name)}")
        if any(dst_name(b.name) == b.name for b in blobs):
            raise SystemExit("ERROR: at least one object has no TRGT token to rename")
        print(f"\nDRY RUN — re-run with --execute to move all {len(blobs)} objects.")
        return

    for p in range(1, args.passes + 1):
        blobs = [b for prefix in PREFIXES for b in bucket.list_blobs(prefix=prefix)]
        if not blobs:
            print("No remaining TRGT objects — migration complete.")
            return
        print(f"Pass {p}: {len(blobs)} source objects remaining")
        done = failed = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(move_blob, bucket, b): b for b in blobs}
            for fut in concurrent.futures.as_completed(futs):
                try:
                    fut.result()
                    done += 1
                    if done % 1000 == 0:
                        print(f"  moved {done}/{len(blobs)}")
                except Exception as e:
                    failed += 1
                    if failed <= 10:
                        print(f"  FAILED {futs[fut].name}: {type(e).__name__}: {e}")
        print(f"Pass {p} done: moved {done}, failed {failed}")
        if failed == 0:
            break

    remaining = sum(1 for prefix in PREFIXES for _ in bucket.list_blobs(prefix=prefix))
    if remaining:
        raise SystemExit(f"ERROR: {remaining} TRGT objects still remain after {args.passes} passes")
    print("Done. All TRGT objects moved to TRGTv5.")


if __name__ == "__main__":
    main()
