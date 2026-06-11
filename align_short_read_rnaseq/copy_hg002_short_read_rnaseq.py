"""One-off Hail Batch job: copy the HG002 short-read RNA-seq (SRR29437757) from ENA into GCS.

Downloads the paired FASTQs on a Hail Batch worker (fast cloud egress from ENA), verifies md5, and
delocalizes them to gs://str-truth-set-v2/raw_data/HG002/rnaseq_short_read/. SRR29437757 is the HG002
lymphoblastoid (LCL) bulk RNA-seq, Illumina NovaSeq X (BioProject PRJNA1124992).
"""

import os

from step_pipeline import pipeline, Backend, files_exist

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
