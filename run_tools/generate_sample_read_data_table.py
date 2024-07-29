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
from str_analysis.utils.file_utils import open_file


FILTER_VCFS_DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:bce1d8d478808ced1bacebfe20ad3226e581c698f3ae8a8f4f72597f9414c5ca"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_INDEX_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

if not os.getcwd().endswith("run_tools"):
	os.chdir("run_tools/")

df_meta = pd.read_table("../20130606_sample_info_1kGP.tsv")

df_meta.set_index("Sample", inplace=True)
sample_id_sex_lookup = dict(df_meta["Gender"])  # "male" or "female"
sample_id_sex_lookup["HG002"] = "male"
sample_id_sex_lookup["HG005"] = "male"
sample_id_sex_lookup["CHM1_CHM13"] = "female"

#%%

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
	"sample_id": "HG002",  # aka. NA24385
	"sequencing_data_type": "element",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.total_depth.txt",
}, {
	# Ultima HG002 downloaded from "https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram"
	"sample_id": "HG002",  # aka. NA24385
	"sequencing_data_type": "ultima",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/ultima/HG002.ultima.cram",      #"https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/ultima/HG002.ultima.cram.crai", #"https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG002/ultima/HG002.ultima.total_depth.txt",
}, {
	#"sample_id": "HG005",  # aka. NA24631  (https://www.coriell.org/1/NIGMS/Collections/NIST-Reference-Materials)
	#"sequencing_data_type": "illumina",
	#"read_data_path":       "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.cram",
	#"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.cram.crai",
	#"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.total_depth.txt",
#}, {
	"sample_id": "NA12878",
	"sequencing_data_type": "illumina",
	"read_data_path":       "gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06/Sample_NA12878/analysis/NA12878.final.cram",
	"read_data_index_path": "gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06/Sample_NA12878/analysis/NA12878.final.cram.crai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/NA12878/illumina/NA12878.pcr_free.total_depth.txt",
	
	# gs://fc-703b4766-6342-4eec-b16f-d9299617a380/gnomad_Misc__wgs_Aug2019/G100862/WGS/NA12878/v3/NA12878.cram
	# gs://fc-703b4766-6342-4eec-b16f-d9299617a380/gnomad_Misc__wgs_Aug2019/G100862/WGS/NA12878/v3/NA12878.cram.crai

}, {
	"sample_id": "CHM1_CHM13",
	"sequencing_data_type": "pacbio",
	"read_data_path": "gs://str-truth-set-v2/raw_data/CHM1_CHM13/pacbio/CHM1_CHM13.bam",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/CHM1_CHM13/pacbio/CHM1_CHM13.bam.bai",
	"depth_stats_path": "gs://str-truth-set-v2/raw_data/CHM1_CHM13/pacbio/CHM1_CHM13.total_depth.txt",
}, {
	"sample_id": "HG005",
	"sequencing_data_type": "pacbio",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG005/pacbio/HG005.downsampled_to_30x.bam",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG005/pacbio/HG005.downsampled_to_30x.bam.bai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG005/pacbio/HG005.downsampled_to_30x.total_depth.txt",
}, {
	"sample_id": "HG002",   # aka. NA24385
	"sequencing_data_type": "pacbio",
	#"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.bam",
	#"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.bam.bai",
	#"read_data_path":		"gs://fc-7891e5cf-0a7a-4c2f-8a18-0d05b27c53ab/GRCh38/PBCCSWholeGenome/NA24385/alignments/NA24385.bam",
	#"read_data_index_path": "gs://fc-7891e5cf-0a7a-4c2f-8a18-0d05b27c53ab/GRCh38/PBCCSWholeGenome/NA24385/alignments/NA24385.bam.bai",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.bam",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.bam.bai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.total_depth.txt",
}, {
	"sample_id": "HG002",   # aka. NA24385
	"sequencing_data_type": "ONT",
	#"read_data_path":       "gs://fc-7891e5cf-0a7a-4c2f-8a18-0d05b27c53ab/GRCh38/ONTWholeGenome/NA24385/alignments/NA24385.bam",
	#"read_data_index_path": "gs://fc-7891e5cf-0a7a-4c2f-8a18-0d05b27c53ab/GRCh38/ONTWholeGenome/NA24385/alignments/NA24385.bam.bai",
	"read_data_path":       "gs://str-truth-set-v2/raw_data/HG002/ONT/HG002.bam",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/ONT/HG002.bam.bai",
	"depth_stats_path":     "gs://str-truth-set-v2/raw_data/HG002/ONT/HG002.total_depth.txt",
}])], ignore_index=True)

for sequencing_data_type, downsampled_bam in [
	("pacbio", "gs://str-truth-set-v2/raw_data/HG005/pacbio/HG005.downsampled_to_30x.bam"),
	("illumina", "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.downsampled_to_20x.bam"),
	("illumina", "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.downsampled_to_10x.bam"),
	("element", "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.downsampled_to_30x.bam"),
	("pacbio", "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_30x.bam"),
	("pacbio", "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_20x.bam"),
	("pacbio", "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_8x.bam"),
	("ONT", "gs://str-truth-set-v2/raw_data/HG002/ONT/HG002.downsampled_to_20x.bam"),
	("ONT", "gs://str-truth-set-v2/raw_data/HG002/ONT/HG002.downsampled_to_8x.bam"),
]:
	df = pd.concat([df, pd.DataFrame([{
		"sample_id": "HG002",
		"sequencing_data_type": sequencing_data_type,
		"read_data_path":       downsampled_bam,
		"read_data_index_path": f"{downsampled_bam}.bai",
		"depth_stats_path":     downsampled_bam.replace(".bam", ".total_depth.txt"),
	}])], ignore_index=True)

bp = pipeline("coverage", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
for _, row in df[~df["depth_stats_path"].apply(lambda p: files_exist([p]))].iterrows():
	if not files_exist([row.read_data_path]):
		print(f"{row.sample_id} read data file {row.read_data_path} not found. Skipping...")
		continue

	stats = hfs.ls(row.read_data_path)
	read_data_size = max(50, int(stats[0].size/10**9))  # at least 50 Gb
	s1 = bp.new_step(f"depth: {row.sample_id} {row.sequencing_data_type}", image=FILTER_VCFS_DOCKER_IMAGE, arg_suffix="depth", cpu=1, storage=f"{read_data_size + 20}Gi")
	s1.switch_gcloud_auth_to_user_account()
	local_fasta, _ = s1.inputs(REFERENCE_FASTA_PATH, REFERENCE_FASTA_INDEX_PATH, localize_by=Localize.COPY)
	local_bam, _ = s1.inputs(row.read_data_path, row.read_data_index_path, localize_by=Localize.GSUTIL_COPY)

	s1.command("set -ex")
	s1.command("curl -L https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -o /usr/local/bin/mosdepth")
	s1.command("chmod 777 /usr/local/bin/mosdepth")

	s1.command("cd /io/")
	s1.command(f"mosdepth -f {local_fasta} -x {row.sample_id}.coverage {local_bam}")
	#s1.output(f"{row.sample_id}.coverage.mosdepth.summary.txt")

	s1.command(f"cat {row.sample_id}.coverage.mosdepth.summary.txt | cut -f 4 | tail -n +2 | head -n 23")
	s1.command(f"grep total {row.sample_id}.coverage.mosdepth.summary.txt > {row.sample_id}.total_depth.txt")
	s1.command(f"cat {row.sample_id}.total_depth.txt")

	s1.output(f"/io/{row.sample_id}.total_depth.txt", row.depth_stats_path)

bp.run()

# download depth stats and add them to the table
def read_depth_stats(path):
	if not files_exist([path]):
		return None

	with open_file(path) as f:
		df = pd.read_table(f, names=["label", "2", "3", "coverage", "5", "6"])
		assert df["label"].iloc[0] == "total"
		depth_of_coverage = df["coverage"].iloc[0]
		return float(depth_of_coverage)

df["depth_of_coverage"] = df["depth_stats_path"].apply(read_depth_stats)
df["male_or_female"] = df["sample_id"].apply(lambda s: sample_id_sex_lookup.get(s))
if sum(df["male_or_female"].isna()) > 0:
	raise ValueError(f"Missing sex metadata for samples {list(df.sample_id)}")

df.sort_values(["sample_id", "sequencing_data_type"], ascending=[False, True], inplace=True)
for _, row in df.iterrows():
	print(f"{float(row.depth_of_coverage):5.1f}x coverage {row.sample_id:<15s} {row.sequencing_data_type:10s} sample: {row.read_data_path} ")

output_table_path = "gs://str-truth-set-v2/HPRC_all_aligned_short_read_and_long_read_samples.tsv"
df.to_csv(os.path.basename(output_table_path), sep="\t", index=False)
os.system(f"gsutil -m cp {os.path.basename(output_table_path)} {output_table_path}")

print(f"Wrote {len(df):,d} rows to {output_table_path}")
for data_type in set(df["sequencing_data_type"]):
	print(f"   {len(df[df['sequencing_data_type'] == data_type]):,d} {data_type} samples")


# missing short read data: HG01123, HG02109,  HG02486, HG02559, NA12878, NA21309
#set(df.sample_id)

#%%