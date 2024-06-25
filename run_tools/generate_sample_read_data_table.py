"""This script creates a table of GRCh38-aligned read data files. The columns are:

sample_id
sequencing_data_type :  "ONT", "pacbio", "illumina", "element", "ultima"
read_data_path
read_data_index_path
"""

import pandas as pd


df = pd.read_table("run_tools/broad_short_read_cram_paths_for_HPRC_samples.txt")
df["sequencing_data_type"] = "illumina"
df.rename(columns={
	"cram_path": "read_data_path",
	"crai_path": "read_data_index_path",
}, inplace=True)


df = df.append([{
	"sample_id": "HG002",  # aka. NA24385  (https://www.coriell.org/1/NIGMS/Collections/NIST-Reference-Materials)
	"sequencing_data_type": "illumina",
	"read_data_path": "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram.crai",
}, {
	"sample_id": "HG002",
	"sequencing_data_type": "element",
	"read_data_path": "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram.crai",
}, {
	"sample_id": "HG002",
	"sequencing_data_type": "ultima",
	"read_data_path": "https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram",
	"read_data_index_path": "https://ultima-ashg-2023-reference-set.s3.amazonaws.com/crams/030945-NA24385-Z0114-CAACATACATCAGAT.cram.crai",
}, {
	"sample_id": "HG005",  # aka. NA24631  (https://www.coriell.org/1/NIGMS/Collections/NIST-Reference-Materials)
	"sequencing_data_type": "illumina",
	"read_data_path": "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.cram",
	"read_data_index_path": "gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.cram.crai",
}, {
	"sample_id": "NA12878",
	"sequencing_data_type": "illumina",
	"read_data_path": "gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06/Sample_NA12878/analysis/NA12878.final.cram",
	"read_data_index_path": "gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06/Sample_NA12878/analysis/NA12878.final.cram.crai",
	# gs://fc-703b4766-6342-4eec-b16f-d9299617a380/gnomad_Misc__wgs_Aug2019/G100862/WGS/NA12878/v3/NA12878.cram
	# gs://fc-703b4766-6342-4eec-b16f-d9299617a380/gnomad_Misc__wgs_Aug2019/G100862/WGS/NA12878/v3/NA12878.cram.crai

}], ignore_index=True)

# missing: "HG002", "HG005, HG01123, HG02109,  HG02486, HG02559, NA12878, NA21309
set(df.sample_id)

# gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/long_read/winnowmap_alignments/HG005vGRCh38_wm_ONT.sort.bam
# gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/long_read/winnowmap_alignments/HG002vGRCh38_wm_ONT.sort.bam

# Revio data - from https://github.com/marbl/HG002/blob/main/Sequencing_data.md
# "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/hifi_revio_pbmay24/hg002v1.0.1_hifi_revio_pbmay24.bam"

#%%