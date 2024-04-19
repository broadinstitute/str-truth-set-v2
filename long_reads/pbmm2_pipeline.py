import os
import pandas as pd
import re
import sys
from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/long-reads@sha256:cff73666379fdf0ab122ee66f614d13dfdff97f99297a563eda74a3f5d08266f"
GATK_DOCKER_IMAGE = "weisburd/gatk:4.3.0.0"

REFERENCE_FASTA = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

TEMP_DIR = "gs://bw2-delete-after-60-days/long-reads"


SAMPLE_METADATA = {
    "CHM13": [
        # from [Vollger 2023]: https://www.ncbi.nlm.nih.gov/sra/SRX5633451[accn]
        "https://sra-pub-src-1.s3.amazonaws.com/SRR9087600/m64011_190228_190319.Q20.fastq.1",
        "https://sra-pub-src-1.s3.amazonaws.com/SRR9087599/m64015_190225_155953.Q20.fastq.1",
        "https://sra-pub-src-1.s3.amazonaws.com/SRR9087598/m64015_190221_025712.Q20.fastq.1",
        "https://sra-pub-src-1.s3.amazonaws.com/SRR9087597/m64015_190224_013150.Q20.fastq.1",
    ],
    "CHM1": [
        # from [Vollger 2023]: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR14407676&display=data-access
        "https://sra-pub-src-2.s3.amazonaws.com/SRR14407676/m64076_200617_023348.ccs.bam.1",
        "https://sra-pub-src-2.s3.amazonaws.com/SRR14407676/m64076_200619_043255.ccs.bam.1",
        "https://sra-pub-src-2.s3.amazonaws.com/SRR14407676/m64076_200621_051102.ccs.bam.1",
        "https://sra-pub-src-2.s3.amazonaws.com/SRR14407676/m64076_200622_112328.ccs.bam.1",
    ],
}

df = pd.read_table("HPRC_AnVIL_samples.tsv")
df = df[df["entity:sample_id"].isin({
    "HG002",
    "HG00438",
    "HG005",
    "HG00621",
    "HG00673",
    "HG00733",
    "HG00735",
    "HG00741",
    "HG01071",
    "HG01106",
    "HG01109",
    "HG01123",
    "HG01175",
    "HG01243",
    "HG01258",
    "HG01358",
    "HG01361",
    "HG01891",
    "HG01928",
    "HG01952",
    "HG01978",
    "HG02055",
    "HG02080",
    "HG02109",
    "HG02145",
    "HG02148",
    "HG02257",
    "HG02486",
    "HG02559",
    "HG02572",
    "HG02622",
    "HG02630",
    "HG02717",
    "HG02723",
    "HG02818",
    "HG02886",
    "HG03098",
    "HG03453",
    "HG03486",
    "HG03492",
    "HG03516",
    "HG03540",
    "HG03579",
    "NA18906",
    "NA19240",
    "NA20129",
    "NA21309",
    "CHM1_CHM13",
    "HG00514",
    "HG03125",
    "NA12878",
})]

for _, row in df.iterrows():
    sample_id = row["entity:sample_id"]
    url_list = eval(row["hifi"])
    if any(".ccs.bam" in u for u in url_list) and any(".fastq" in u for u in url_list):
        # if a sample has both .ccs.bams and .fastq, only keep the .ccs.bams
        url_list = [u for u in url_list if u.endswith(".ccs.bam")]

    SAMPLE_METADATA[sample_id] = url_list
    if len(url_list) >= 3:
        print(f"{sample_id} has {len(SAMPLE_METADATA[sample_id])} PacBio files:")
        for url in url_list:
            print(f"  {url}")

print(f"Processing PacBio data for {len(SAMPLE_METADATA)} HPRC samples")
#sys.exit(0)

def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("-s", "--sample-id", action="append", help="only process the given sample id(s)", choices=SAMPLE_METADATA.keys())
    args = bp.parse_known_args()

    s0 = bp.new_step(f"pbmm2: index hg38",
                     arg_suffix=f"index",
                     step_number=0,
                     image=DOCKER_IMAGE,
                     cpu=4,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=TEMP_DIR)

    local_fasta = s0.input(REFERENCE_FASTA, localize_by=Localize.COPY)
    print(f"Indexing {local_fasta}")
    s0.command(f"pbmm2 index {local_fasta} Homo_sapiens_assembly38.mmi --preset SUBREAD")
    s0.output("Homo_sapiens_assembly38.mmi")
    s0.output(f"{local_fasta}")

    aligned_bam_files_for_CHM1_CHM13 = []
    align_bam_files_for_CHM1_CHM13_steps = []
    for sample_id, unaligned_reads_urls in SAMPLE_METADATA.items():
        if args.sample_id and sample_id not in args.sample_id:
            continue
        
        # if the data is @ SRA / EBI, download it to a gs:// bucket
        download_steps = []
        gs_paths = []
        for url in unaligned_reads_urls:
            if url.startswith("gs://"):
                gs_paths.append(url)
                continue

            s1 = bp.new_step(
                f"pbmm2: download {os.path.basename(url)}",
                arg_suffix=f"download",
                step_number=1,
                image=DOCKER_IMAGE,
                cpu=0.5,
                storage="50Gi",
                output_dir=TEMP_DIR,
                delocalize_by=Delocalize.GSUTIL_COPY,
            )

            local_filename = re.sub("[.]1$", "", os.path.basename(url))
            print(f"{sample_id}: Downloading {local_filename}")
            s1.command("set -ex")
            s1.command("cd /io")
            s1.command(f"wget {url} -O {local_filename}")
            s1.output(f"/io/{local_filename}")

            download_steps.append(s1)
            gs_paths.append(os.path.join(TEMP_DIR, local_filename))

        # align reads
        s2 = bp.new_step(
            f"pbmm2: align {sample_id}",
            arg_suffix=f"align",
            step_number=2,
            image=DOCKER_IMAGE,
            cpu=16,
            memory="highmem",
            storage="500Gi",
            output_dir=TEMP_DIR,
            #delocalize_by=Delocalize.GSUTIL_COPY,
        )
        s2.depends_on(download_steps)

        local_fasta = s2.input(REFERENCE_FASTA, localize_by=Localize.COPY)
        local_read_data_paths = s2.inputs(*gs_paths, localize_by=Localize.COPY)

        s2.command("set -ex")
        s2.command("cd /io")
        for local_path in local_read_data_paths:
            s2.command(f"echo '{local_path}' >> read_data_paths.fofn")

        if all(".ccs.bam" in path for path in gs_paths):
            preset = "HIFI"
        elif all(".ccs.bam" not in path for path in gs_paths):
            preset = "SUBREAD"
        else:
            raise ValueError(f"{sample_id} has a mix of .ccs.bam and .fastq files. Expecting one or the other")

        s2.command(f"pbmm2 align --sort --strip --preset {preset} {local_fasta} read_data_paths.fofn {sample_id}.aligned.bam")
        s2.output(f"/io/{sample_id}.aligned.bam")
        s2.output(f"/io/{sample_id}.aligned.bam.bai")

        # compute depth and other stats
        s3 = bp.new_step(
            f"pbmm2: depth of coverage {sample_id}",
            arg_suffix=f"stats",
            step_number=3,
            image=DOCKER_IMAGE,
            cpu=1,
            memory="standard",
            output_dir=TEMP_DIR,
            localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
        )
        local_bam_path, _ = s3.use_previous_step_outputs_as_inputs(s2)
        s3.command("set -ex")
        s3.command(f"mosdepth -x {sample_id} {local_bam_path}")
        s3.command("ls -lh")
        s3.command(f"cat {sample_id}.mosdepth.summary.txt")

        s3.output(f"{sample_id}.mosdepth.summary.txt")
        s3.output(f"{sample_id}.mosdepth.global.dist.txt")

        # downsample to 30x
        s4 = bp.new_step(
            f"pbmm2: downsample {sample_id}",
            arg_suffix=f"downsample",
            step_number=4,
            image=GATK_DOCKER_IMAGE,
            cpu=2,
            memory="highmem",
            storage="200Gi",
            output_dir=TEMP_DIR,
            localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
        )
        s4.depends_on(s3)

        local_fasta =  s4.input(REFERENCE_FASTA)
        local_bam, _ = s4.inputs(os.path.join(TEMP_DIR, f"{sample_id}.aligned.bam"),
                                 os.path.join(TEMP_DIR, f"{sample_id}.aligned.bam.bai"))

        local_mosdepth_summary = s4.input(os.path.join(TEMP_DIR, f"{sample_id}.mosdepth.summary.txt"))

        if sample_id in ("CHM1", "CHM13"):
            target_coverage = 15
        else:
            target_coverage = 30

        output_bam_filename = f"{sample_id}.downsampled_to_{target_coverage}x.bam"
        if sample_id in ("CHM1", "CHM13"):
            aligned_bam_files_for_CHM1_CHM13.append(output_bam_filename)
            align_bam_files_for_CHM1_CHM13_steps.append(s4)

        s4.command("set -ex")
        s4.command("cd /io")

        s4.command(f"time gatk --java-options '-Xmx11G' DownsampleSam "
                   f"REFERENCE_SEQUENCE={local_fasta} "
                   f"I={local_bam} "
                   f"O={output_bam_filename} "
                   f"""P=$(echo "{target_coverage} / $(grep total {local_mosdepth_summary} | cut -f 4)" | bc -l | awk '{{printf "%.4f", $0}}') """
                   f"CREATE_INDEX=true")

        s4.command(f"mv {output_bam_filename.replace('.bam', '.bai')} {output_bam_filename}.bai")  # rename the .bai file
        s4.command("ls -lh")

        #s4.command(f"gatk CollectWgsMetrics STOP_AFTER={5*10**7} I={output_bam_filename} O=metrics.txt R={local_fasta}")
        #s4.command(f"cat metrics.txt | head -n 8 | tail -n 2")

        s4.output(f"/io/{output_bam_filename}")
        s4.output(f"/io/{output_bam_filename}.bai")

    if args.sample_id:
        bp.run()        
        return
    
    # use samtools to merge bam files in aligned_bam_files_for_CHM1_CHM13
    merged_sample_id = "CHM1_CHM13"
    s5 = bp.new_step(
        f"pbmm2: merge {merged_sample_id}",
        arg_suffix=f"merge",
        step_number=5,
        image=DOCKER_IMAGE,
        cpu=2,
        memory="highmem",
        storage="200Gi",
        output_dir=TEMP_DIR,
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
    )
    for s in align_bam_files_for_CHM1_CHM13_steps:
        s5.depends_on(s)

    local_bam_paths = [
        s5.input(os.path.join(TEMP_DIR, bam_file))
        for bam_file in aligned_bam_files_for_CHM1_CHM13
    ]
    s5.command("set -ex")
    s5.command("cd /io")

    merged_bam_filename = f"{merged_sample_id}.aligned.bam"
    s5.command(f"samtools merge -t 2 -f {merged_bam_filename} {' '.join(map(str, local_bam_paths))}")
    s5.command(f"samtools index {merged_bam_filename}")

    s5.output(merged_bam_filename)
    s5.output(f"{merged_bam_filename}.bai")

    # compute depth and other stats
    s6 = bp.new_step(
        f"pbmm2: depth of coverage {sample_id}",
        arg_suffix=f"stats3",
        step_number=6,
        image=DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        output_dir=TEMP_DIR,
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
    )
    local_bam_path, _ = s6.use_previous_step_outputs_as_inputs(s5)
    s6.command("set -ex")
    s6.command(f"mosdepth -x {merged_sample_id} {local_bam_path}")
    s6.command("ls -lh")
    s6.command(f"cat {merged_sample_id}.mosdepth.summary.txt")

    s6.output(f"{merged_sample_id}.mosdepth.summary.txt")
    s6.output(f"{merged_sample_id}.mosdepth.global.dist.txt")

    bp.run()


if __name__ == "__main__":
    main()

