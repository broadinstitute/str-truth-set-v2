"""One-off Hail Batch job: copy the HG002 short-read RNA-seq (SRR29437757) from ENA into GCS.

Downloads the paired FASTQs on a Hail Batch worker (fast cloud egress from ENA), verifies md5, and
delocalizes them to gs://str-truth-set-v2/raw_data/HG002/rnaseq_short_read/. SRR29437757 is the HG002
lymphoblastoid (LCL) bulk RNA-seq, Illumina NovaSeq X (BioProject PRJNA1124992).
"""

import os

import hailtop.batch as hb
from step_pipeline import pipeline, Backend, files_exist

# step_pipeline calls Batch.run(wait=True, disable_progress_bar=False), which deadlocks under
# hailtop's nest_asyncio event loop (the submit hangs in select() forever). A plain non-blocking
# hailtop submit works, so force wait=False here: the batch is still submitted and runs in the
# cloud; we just don't block the local process waiting for it to finish.
_orig_batch_run = hb.Batch.run
def _nonblocking_batch_run(self, *args, **kwargs):
    kwargs["wait"] = False
    kwargs["disable_progress_bar"] = True
    kwargs["open"] = False
    return _orig_batch_run(self, *args, **kwargs)
hb.Batch.run = _nonblocking_batch_run

# long-reads image (has wget + coreutils md5sum); same digest as the Iso-Seq alignment pipeline.
DOCKER_IMAGE = "weisburd/long-reads@sha256:99cfef9b5ff7562ddfb57986a966aa7c63be42c93d6c49b5ca3444db23668ac2"

OUTPUT_DIR = "gs://str-truth-set-v2/raw_data/HG002/rnaseq_short_read"

ENA_PREFIX = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/057/SRR29437757"
FASTQS = {
    "SRR29437757_1.fastq.gz": "34ec659eb9c229cbb94146c286c5885a",
    "SRR29437757_2.fastq.gz": "d378f79a7c2351c6df1217edce46de18",
}


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    bp.parse_known_args()
    bp.set_name("copy HG002 short-read RNA-seq (SRR29437757)")

    if files_exist([os.path.join(OUTPUT_DIR, f) for f in FASTQS]):
        print("FASTQs already present in", OUTPUT_DIR)
        return

    s = bp.new_step(
        "copy SRR29437757 fastqs",
        arg_suffix="copy",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=2,
        storage="50Gi",
        output_dir=OUTPUT_DIR,
    )
    s.command("set -ex")
    s.command("cd /io")
    for fastq, md5 in FASTQS.items():
        s.command(f"wget -nv {ENA_PREFIX}/{fastq} -O {fastq}")
        s.command(f"echo '{md5}  {fastq}' | md5sum -c -")
        s.output(f"/io/{fastq}")

    bp.run()


if __name__ == "__main__":
    main()
