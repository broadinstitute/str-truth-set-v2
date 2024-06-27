"""This script creates a table of GRCh38-aligned read data files. The columns are:

sample_id
sequencing_data_type :  "ONT", "pacbio", "illumina", "element", "ultima"
read_data_path
read_data_index_path
"""

import hailtop.fs as hfs
import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize, files_exist

FILTER_VCFS_DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:bce1d8d478808ced1bacebfe20ad3226e581c698f3ae8a8f4f72597f9414c5ca"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_INDEX_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

df = pd.read_table("broad_short_read_cram_paths_for_HPRC_samples.txt")
df["sequencing_data_type"] = "illumina"
df.rename(columns={
	"cram_path": "read_data_path",
	"crai_path": "read_data_index_path",
}, inplace=True)

df["depth_stats_path"] = df["sample_id"].apply(
	lambda sample_id: os.path.join(f"gs://str-truth-set-v2/raw_data/{sample_id}/illumina", f"{sample_id}.total_depth.txt")) 

df = pd.concat([df, pd.DataFrame([{
	"sample_id": "HG002",  # aka. NA24385  (https://www.coriell.org/1/NIGMS/Collections/NIST-Reference-Materials)
	"sequencing_data_type": "illumina",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.total_depth.txt",
}, {
	"sample_id": "HG002",
	"sequencing_data_type": "element",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.total_depth.txt",
}, {
	# Ultima HG002 downloaded from "https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram"
	"sample_id": "HG002",
	"sequencing_data_type": "ultima",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/ultima/HG002.ultima.cram",      #"https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/ultima/HG002.ultima.cram.crai", #"https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG002/ultima/HG002.ultima.total_depth.txt",
}, {
	"sample_id": "HG005",  # aka. NA24631  (https://www.coriell.org/1/NIGMS/Collections/NIST-Reference-Materials)
	"sequencing_data_type": "illumina",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.total_depth.txt",
}, {
	"sample_id": "NA12878",
	"sequencing_data_type": "illumina",
	"read_data_path":       "gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06/Sample_NA12878/analysis/NA12878.final.cram",
	"read_data_index_path": "gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06/Sample_NA12878/analysis/NA12878.final.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/NA12878/illumina/NA12878.pcr_free.total_depth.txt",
	
	# gs://fc-703b4766-6342-4eec-b16f-d9299617a380/gnomad_Misc__wgs_Aug2019/G100862/WGS/NA12878/v3/NA12878.cram
	# gs://fc-703b4766-6342-4eec-b16f-d9299617a380/gnomad_Misc__wgs_Aug2019/G100862/WGS/NA12878/v3/NA12878.cram.crai

}])], ignore_index=True)

bp = pipeline("coverage", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
for _, row in df.iterrows():
	if not files_exist([row.read_data_path]):
		print(f"{row.sample_id} read data file {row.read_data_path} not found. Skipping...")
		continue

	stats = hfs.ls(row.read_data_path)
	read_data_size = max(50, int(stats[0].size/10**9))  # at least 50 Gb
	s1 = bp.new_step(f"depth: {row.sample_id}", image=FILTER_VCFS_DOCKER_IMAGE, arg_suffix="depth", cpu=1, storage=f"{read_data_size + 20}Gi")
	s1.switch_gcloud_auth_to_user_account()
	local_fasta, _ = s1.inputs(REFERENCE_FASTA_PATH, REFERENCE_FASTA_INDEX_PATH, localize_by=Localize.COPY)
	local_bam, _ = s1.inputs(row.read_data_path, row.read_data_index_path, localize_by=Localize.COPY)

	s1.command("set -ex")
	s1.command("curl -L https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -o /usr/local/bin/mosdepth")
	s1.command("chmod 777 /usr/local/bin/mosdepth")

	s1.command("cd /io/")
	s1.command(f"mosdepth -f {local_fasta} -x {row.sample_id}.coverage {local_bam}")
	#s1.output(f"{row.sample_id}.coverage.mosdepth.summary.txt")

	s1.command(f"cat {row.sample_id}.coverage.mosdepth.summary.txt | cut -f 4 | tail -n +2 | head -n 23")
	s1.command(f"grep total {row.sample_id}.coverage.mosdepth.summary.txt > {row.sample_id}.total_coverage.txt")
	s1.command(f"cat {row.sample_id}.total_coverage.txt")

	s1.output(f"/io/{row.sample_id}.total_coverage.txt", row.depth_stats_path)

bp.run()

# download depth stats and add them to the table
df["depth_of_coverage"] = df["depth_stats_path"].apply(
	lambda path: pd.read_table(path, names=["1", "2", "3", "coverage", "5", "6"])["coverage"].iloc[0] if files_exist([path]) else None)


# launch mosdepth
# missing: "HG002", "HG005, HG01123, HG02109,  HG02486, HG02559, NA12878, NA21309
set(df.sample_id)

# gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/long_read/winnowmap_alignments/HG005vGRCh38_wm_ONT.sort.bam
# gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/long_read/winnowmap_alignments/HG002vGRCh38_wm_ONT.sort.bam

# Revio data - from https://github.com/marbl/HG002/blob/main/Sequencing_data.md
# "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/hifi_revio_pbmay24/hg002v1.0.1_hifi_revio_pbmay24.bam"

#%%