"""Exapmle step pipeline that converts CRAM/BAM file(s) to FASTQ format."""

import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize, files_exist


DOCKER_IMAGE = "weisburd/bwa@sha256:f0457b9fd9b6ad6ba1337cfad07fcc3bd6af241f5842cd701784c5e297ceb921"

REFERENCE_FASTA_PATH = {
    "37": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta",
    #"38": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
    "38": "gs://str-truth-set/hg38/ref/hg38.fa",
    "t2t": "gs://gcp-public-data--broad-references/t2t/v2/chm13v2.0.maskedY.fasta",
}
REFERENCE_FASTA_FAI_PATH = {
    "37": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
    #"38": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
    "38": "gs://str-truth-set/hg38/ref/hg38.fa.fai",
    "t2t": "gs://gcp-public-data--broad-references/t2t/v2/chm13v2.0.maskedY.fasta.fai",
}

def main():
    batch_pipeline = pipeline(backend=Backend.HAIL_BATCH_SERVICE)
    p = batch_pipeline.get_config_arg_parser()
    p.add_argument("-o", "--output-cram", required=True, help="gs:// path of output CRAM")
    p.add_argument("-g", "--reference", choices={"37", "38", "t2t"}, default="38", help="Reference genome fasta.")
    p.add_argument("fastq_r1", help="Input FASTQ file _R1")
    p.add_argument("fastq_r2", help="Input FASTQ file _R2")
    args, _ = p.parse_known_args()

    batch_pipeline.set_name(f"align to {args.reference}")

    if not args.output_cram.endswith(".cram"):
        p.error("--output-cram path must end with .cram")

    cram_filename = os.path.basename(args.output_cram)
    cram_prefix = re.sub(".cram$", "", cram_filename)
    output_dir = os.path.dirname(args.output_cram)

    cpus = 8
    s1 = batch_pipeline.new_step(f"bwa mem: {cram_filename}", step_number=1, arg_suffix="bwa",
        image=DOCKER_IMAGE, cpu=cpus, memory="standard", storage="300G")
    #s1.switch_gcloud_auth_to_user_account()

    local_r1_fastq, local_r2_fastq = s1.inputs(args.fastq_r1, args.fastq_r2, localize_by=Localize.COPY)

    reference_fasta_path = REFERENCE_FASTA_PATH[args.reference]
    reference_fai_path = REFERENCE_FASTA_FAI_PATH[args.reference]

    local_reference_fasta, _, _, _, _, _, _  = s1.inputs([
        reference_fasta_path,
        f"{reference_fasta_path}.fai",
        f"{reference_fasta_path}.amb",
        f"{reference_fasta_path}.ann",
        f"{reference_fasta_path}.bwt",
        f"{reference_fasta_path}.pac",
        f"{reference_fasta_path}.sa",
    ], localize_by=Localize.COPY)

    s1.command("set -ex")
    s1.command("cd /io")
    #s1.command(f"bwa index {local_reference_fasta}")
    s1.command(f"bwa mem -t {cpus} {local_reference_fasta} {local_r1_fastq} {local_r2_fastq} | "
               f"samtools sort --reference {local_reference_fasta} -o {cram_filename} - ")
    s1.command(f"samtools index {cram_filename}")

    s1.output(cram_filename, output_path=args.output_cram)
    s1.output(f"{cram_filename}.crai", output_path=f"{args.output_cram}.crai")

    print(f"Output: {args.output_cram}")

    # measure coverage
    s2 = batch_pipeline.new_step(
        f"coverage: {cram_filename}", step_number=2, image=DOCKER_IMAGE, cpu=1, memory="standard", storage="20Gi", output_dir=output_dir)
    s2.depends_on(s1)

    local_fasta = s2.input(reference_fasta_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    s2.input(reference_fai_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    local_cram = s2.input(args.output_cram, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    s2.command("set -ex")
    s2.command("wget https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -O /usr/local/bin/mosdepth")
    s2.command("chmod 777 /usr/local/bin/mosdepth")

    s2.command(f"mosdepth -f {local_fasta} -x {cram_prefix}.coverage {local_cram}")
    #s2.output(f"{bam_or_cram_prefix}.coverage.mosdepth.summary.txt")

    s2.command(f"cat {cram_prefix}.coverage.mosdepth.summary.txt | cut -f 4 | tail -n +2 | head -n 23")
    s2.command(f"grep total {cram_prefix}.coverage.mosdepth.summary.txt > {cram_prefix}.total_coverage.txt")
    s2.command(f"cat {cram_prefix}.total_coverage.txt")

    s2.output(f"{cram_prefix}.total_coverage.txt")

    batch_pipeline.run()



if __name__ == "__main__":
    main()

