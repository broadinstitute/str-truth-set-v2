"""One-off pipeline that aligns the public GIAB HG002 PacBio Iso-Seq read sets to hg38.

Unlike pbmm2_pipeline.py (which aligns WGS HiFi reads with the SUBREAD/HIFI preset and then downsamples
to a target genome coverage), Iso-Seq reads are spliced mRNA transcripts. They are aligned here with the
splice-aware `--preset ISOSEQ`, and the coverage/downsample steps are intentionally omitted because genome
coverage is meaningless for transcript data.

Input read sets are the unaligned (CCS) Iso-Seq BAMs released by GIAB/Coriell on the human-pangenomics S3
bucket (3 cell sources for HG002). Each is aligned to hg38 and indexed.
"""

import os

from step_pipeline import pipeline, Backend, Localize, files_exist

# weisburd/long-reads image rebuilt with pbmm2 26.1.0 (see align_pacbio/docker_isoseq/Dockerfile).
# Pinned to the digest pushed by the build_isoseq_image GitHub Actions workflow.
DOCKER_IMAGE = "weisburd/long-reads:pbmm2_26.1.0"

REFERENCE_FASTA = "gs://str-truth-set/hg38/ref/hg38.fa"

OUTPUT_DIR = "gs://str-truth-set-v2/raw_data/HG002/pacbio_isoseq"

S3_PREFIX = ("https://human-pangenomics.s3.amazonaws.com/submissions/"
             "35540CB6-EEFE-465B-8A6E-7D3EF7AC8900--HG002-IsoSeq")

# sample_id -> unaligned Iso-Seq CCS BAM on the human-pangenomics S3 bucket
SAMPLE_METADATA = {
    "HG002-NA24385-LCL":           f"{S3_PREFIX}/HG002-NA24385-LCL-collapsed_hifi_reads.bam",
    "HG002-NA26105-iPSC_from_LCL": f"{S3_PREFIX}/HG002-NA26105-iPSC_from_LCL-collapsed_hifi_reads.bam",
    "HG002-NA27730-iPSC_from_PBMC": f"{S3_PREFIX}/HG002-NA27730-iPSC_from_PBMC-collapsed_hifi_reads.bam",
}


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("-s", "--sample-id", action="append", help="only process the given sample id(s)",
                        choices=SAMPLE_METADATA.keys())
    args = bp.parse_known_args()

    if args.sample_id:
        bp.set_name(f"pbmm2 Iso-Seq align pipeline: {len(args.sample_id)} samples")
    else:
        bp.set_name(f"pbmm2 Iso-Seq align pipeline: all {len(SAMPLE_METADATA)} samples")

    for sample_id, reads_url in SAMPLE_METADATA.items():
        if args.sample_id and sample_id not in args.sample_id:
            continue

        output_bam_filename = f"{sample_id}.hg38.bam"
        if files_exist([
            os.path.join(OUTPUT_DIR, output_bam_filename),
            os.path.join(OUTPUT_DIR, f"{output_bam_filename}.bai"),
        ]):
            continue

        s1 = bp.new_step(
            f"pbmm2 isoseq: align {sample_id}",
            arg_suffix="align",
            step_number=1,
            image=DOCKER_IMAGE,
            cpu=16,
            memory="highmem",
            storage="250Gi",
            output_dir=OUTPUT_DIR,
        )

        local_fasta = s1.input(REFERENCE_FASTA, localize_by=Localize.COPY)

        s1.command("set -ex")
        s1.command("cd /io")
        s1.command(f"wget -nv {reads_url} -O {sample_id}.unaligned.bam")
        s1.command(f"pbmm2 align --preset ISOSEQ --sort --num-threads 16 "
                   f"{local_fasta} {sample_id}.unaligned.bam {output_bam_filename}")
        s1.command(f"samtools index {output_bam_filename}")
        s1.command("ls -lh")

        s1.output(f"/io/{output_bam_filename}")
        s1.output(f"/io/{output_bam_filename}.bai")

    bp.run()


if __name__ == "__main__":
    main()
