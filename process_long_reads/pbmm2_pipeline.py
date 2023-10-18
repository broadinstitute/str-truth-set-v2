import collections
import os
import pandas as pd
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

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

#for _, row in df.iterrows():
#    sample_id = row["entity:sample_id"]
#    url_list = eval(row["hifi"])
#    SAMPLE_METADATA[sample_id] = url_list
#    if len(url_list) >= 5:
#        print(f"{sample_id} has {len(SAMPLE_METADATA[sample_id])} PacBio files:")
#        for url in url_list:
#            print(f"  {url}")

print(f"Processing PacBio data for {len(SAMPLE_METADATA)} HPRC samples")


DOCKER_IMAGE = "weisburd/process-long-reads@sha256:9cc49dac607d6aa08147563daccf94726861c5f90809bc7bea44b4d1999aff87"

TEMP_DIR = "gs://bw2-delete-after-30-days/long-reads"
OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/hipstr"

REFERENCE_FASTA = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    #parser = bp.get_config_arg_parser()
    #args = bp.parse_known_args()

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


    for sample_id, unaligned_reads_urls in SAMPLE_METADATA.items():
        download_steps = []
        gs_paths = collections.defaultdict(list)
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
            s1.command("cd /io")
            s1.command(f"wget {url} -O {local_filename}")
            s1.output(f"/io/{local_filename}")

            download_steps.append(s1)
            gs_paths.append(os.path.join(TEMP_DIR, local_filename))

        # align gs_paths
        s2 = bp.new_step(
            f"pbmm2: align {sample_id}",
            arg_suffix=f"align",
            step_number=2,
            image=DOCKER_IMAGE,
            cpu=8,
            memory="highmem",
            storage="250Gi",
            output_dir=TEMP_DIR,
            delocalize_by=Delocalize.GSUTIL_COPY,
        )
        s2.depends_on(download_steps)

        local_fasta = s2.input(REFERENCE_FASTA, localize_by=Localize.COPY)
        local_read_data_paths = s2.inputs(gs_paths, localize_by=Localize.COPY)

        s2.command("cd /io")
        for local_path in local_read_data_paths:
            s2.command(f"echo '{local_path}' >> read_data_paths.txt")

        s2.command(f"pbmm2 align --sort --preset SUBREAD {local_fasta} read_data_paths.txt {sample_id}.subreads.bam")
        s2.output(f"/io/{sample_id}.subreads.bam")
        s2.output(f"/io/{sample_id}.subreads.bam.bai")

    bp.run()


if __name__ == "__main__":
    main()

