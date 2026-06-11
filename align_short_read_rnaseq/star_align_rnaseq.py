"""Hail Batch pipeline that aligns the HG002 short-read RNA-seq (SRR29437757) to hg38 with STAR.

Companion to copy_hg002_short_read_rnaseq.py, which downloads the paired FASTQs from ENA into
gs://str-truth-set-v2/raw_data/HG002/rnaseq_short_read/. SRR29437757 is the HG002 lymphoblastoid (LCL)
bulk RNA-seq, Illumina NovaSeq X, 149bp paired-end (BioProject PRJNA1124992).

Two steps:
  1. Build a STAR index for hg38 + GENCODE v49 (sjdbOverhang = read_length - 1). Skipped if it already
     exists on GCS, so re-runs reuse it.
  2. STAR 2-pass alignment (GTEx v8 parameters) -> coordinate-sorted BAM + per-gene read counts
     (--quantMode GeneCounts) + splice junctions (SJ.out.tab), then samtools index.

There is no public hg38 BAM for SRR29437757 (SRA/ENA host raw FASTQ only), so it has to be aligned here.
"""

import os

from step_pipeline import pipeline, Backend, Localize, files_exist

# STAR 2.7.10b + samtools (see tgg-rnaseq-pipelines/star/docker/Dockerfile).
DOCKER_IMAGE = "weisburd/star@sha256:6c0da40e33f50341e4dc3dc428df23b8036f25be488999665f004ac8b39007c9"

REFERENCE_FASTA = "gs://str-truth-set/hg38/ref/hg38.fa"

OUTPUT_DIR = "gs://str-truth-set-v2/raw_data/HG002/rnaseq_short_read/star"
REF_DIR = f"{OUTPUT_DIR}/ref"
# GENCODE v49 GTF (GRCh38, chr-prefixed). Staged here once from the canonical release:
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
GENCODE_V49_GTF = f"{REF_DIR}/gencode.v49.annotation.gtf.gz"

SAMPLE_ID = "SRR29437757"
FASTQ_DIR = "gs://str-truth-set-v2/raw_data/HG002/rnaseq_short_read"
FASTQ_R1 = f"{FASTQ_DIR}/{SAMPLE_ID}_1.fastq.gz"
FASTQ_R2 = f"{FASTQ_DIR}/{SAMPLE_ID}_2.fastq.gz"
READ_LENGTH = 149

# Files written by STAR --runMode genomeGenerate (used to detect an existing index and to localize it).
STAR_INDEX_FILES = [
    "Genome", "Log.out", "SA", "SAindex", "chrLength.txt", "chrName.txt", "chrNameLength.txt",
    "chrStart.txt", "exonGeTrInfo.tab", "exonInfo.tab", "geneInfo.tab", "genomeParameters.txt",
    "sjdbInfo.txt", "sjdbList.fromGTF.out.tab", "sjdbList.out.tab", "transcriptInfo.tab",
]


def add_star_index_step(bp, index_dir):
    """Build a STAR index for hg38 + GENCODE v49. Returns the step (outputs the index files to index_dir)."""
    s = bp.new_step(
        f"STAR index: gencode v49, {READ_LENGTH}bp reads",
        arg_suffix="star-index",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=16,
        memory="highmem",
        storage="75Gi",
        output_dir=index_dir,
    )
    local_fasta = s.input(REFERENCE_FASTA, localize_by=Localize.COPY)
    local_gtf_gz = s.input(GENCODE_V49_GTF, localize_by=Localize.COPY)

    s.command("set -ex")
    s.command("cd /io")
    s.command(f"gunzip -c {local_gtf_gz} > gencode.v49.annotation.gtf")
    s.command("mkdir STAR_index")
    s.command(
        "STAR --runMode genomeGenerate "
        "--runThreadN 16 "
        "--genomeDir STAR_index "
        f"--genomeFastaFiles {local_fasta} "
        "--sjdbGTFfile gencode.v49.annotation.gtf "
        f"--sjdbOverhang {READ_LENGTH - 1}"
    )
    s.command("ls -lh STAR_index")
    for f in STAR_INDEX_FILES:
        s.output(f"/io/STAR_index/{f}")
    return s


def add_star_align_step(bp, index_dir, index_step):
    """Run STAR 2-pass alignment of the SRR29437757 FASTQs. index_step is None if the index already exists."""
    s = bp.new_step(
        f"STAR align: {SAMPLE_ID}",
        arg_suffix="star-align",
        step_number=2,
        image=DOCKER_IMAGE,
        cpu=16,
        memory="highmem",
        storage="200Gi",
        output_dir=OUTPUT_DIR,
    )

    if index_step is not None:
        # Index built in this run: wire the batch-level dependency and localize its outputs directly.
        index_inputs = s.use_previous_step_outputs_as_inputs(index_step, localize_by=Localize.COPY)
        s.depends_on(index_step)
    else:
        # Index already on GCS: copy it onto the worker.
        index_inputs = s.inputs(
            [os.path.join(index_dir, f) for f in STAR_INDEX_FILES], localize_by=Localize.COPY)
    local_index_dir = index_inputs[0].local_dir

    local_r1, local_r2 = s.inputs([FASTQ_R1, FASTQ_R2], localize_by=Localize.COPY)

    s.command("set -ex")
    s.command("cd /io")
    s.command(
        "STAR --runMode alignReads "
        "--twopassMode Basic "
        f"--genomeDir {local_index_dir} "
        "--runThreadN 16 "
        "--readFilesCommand zcat "
        f"--readFilesIn {local_r1} {local_r2} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterType BySJout "
        f"--outFileNamePrefix {SAMPLE_ID}. "
        "--outSAMunmapped Within "
        "--outSAMattributes NH HI AS nM NM ch "
        f"--outSAMattrRGline ID:rg1 SM:{SAMPLE_ID} "
        "--quantMode GeneCounts "
        "--alignSJoverhangMin 8 "
        "--alignSJDBoverhangMin 1 "
        "--alignIntronMin 20 "
        "--alignIntronMax 1000000 "
        "--alignMatesGapMax 1000000 "
        "--limitBAMsortRAM 50000000000"
    )
    s.command(f"samtools index {SAMPLE_ID}.Aligned.sortedByCoord.out.bam")
    s.command(f"gzip {SAMPLE_ID}.SJ.out.tab {SAMPLE_ID}.ReadsPerGene.out.tab")
    s.command("ls -lh")

    s.output(f"/io/{SAMPLE_ID}.Aligned.sortedByCoord.out.bam")
    s.output(f"/io/{SAMPLE_ID}.Aligned.sortedByCoord.out.bam.bai")
    s.output(f"/io/{SAMPLE_ID}.SJ.out.tab.gz")
    s.output(f"/io/{SAMPLE_ID}.ReadsPerGene.out.tab.gz")
    s.output(f"/io/{SAMPLE_ID}.Log.final.out")
    return s


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    bp.parse_known_args()
    bp.set_name(f"STAR align HG002 RNA-seq ({SAMPLE_ID}) to hg38")

    output_bam = f"{OUTPUT_DIR}/{SAMPLE_ID}.Aligned.sortedByCoord.out.bam"
    if files_exist([output_bam, f"{output_bam}.bai"]):
        print("Aligned BAM already present:", output_bam)
        return

    index_dir = f"{REF_DIR}/STAR_index_gencode_v49_{READ_LENGTH}bp"
    index_step = None
    if files_exist([os.path.join(index_dir, f) for f in ("SA", "SAindex", "Genome")]):
        print("Reusing existing STAR index:", index_dir)
    else:
        index_step = add_star_index_step(bp, index_dir)

    add_star_align_step(bp, index_dir, index_step)
    bp.run()
    print("Output BAM:", output_bam)


if __name__ == "__main__":
    main()
