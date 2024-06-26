
#%%

python3 align_reads_using_bwa.py \
  gs://str-truth-set-v2/raw_data/HG002/illumina/HG002_HiSeq30x_subsampled_R1.fastq.gz \
  gs://str-truth-set-v2/raw_data/HG002/illumina/HG002_HiSeq30x_subsampled_R2.fastq.gz \
  --use-non-preemptibles \
  -o gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram

#python3 align_reads_using_bwa.py \
#  https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz \
#  https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz \
#  --use-non-preemptibles \
#  -o gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram

python3 align_reads_using_bwa.py \
  gs://str-truth-set-v2/raw_data/HG002/element/GAT-LI-C044_R1.fastq.gz \
  gs://str-truth-set-v2/raw_data/HG002/element/GAT-LI-C044_R2.fastq.gz \
  --use-non-preemptibles \  
-o gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram

#python3 align_reads_using_bwa.py \
#  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/element/trio/HG002/ins500_600/GAT-LI-C044_R1.fastq.gz \
#  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/element/trio/HG002/ins500_600/GAT-LI-C044_R2.fastq.gz \
#  --use-non-preemptibles \  
#  -o gs://str-truth-set-v2/raw_data/HG002/element/HG002.element.cram

python3 align_reads_using_bwa.py \
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/028/SRR14724528/SRR14724528_1.fastq.gz \
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/028/SRR14724528/SRR14724528_2.fastq.gz \
  --use-non-preemptibles \  
  -o gs://str-truth-set-v2/raw_data/HG005/illumina/HG005.pcr_free.cram  # novoseq


# HG005 HiSeq-X: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR14724509&display=data-access
#"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/009/SRR14724509/SRR14724509_1.fastq.gz"
#"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/009/SRR14724509/SRR14724509_2.fastq.gz"

# HG005 novoseq: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR14724528&display=data-access
#"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/028/SRR14724528/SRR14724528_1.fastq.gz"
#"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/028/SRR14724528/SRR14724528_2.fastq.gz"



# missing: "HG002", "HG005, HG01123, HG02109,  HG02486, HG02559, NA12878, NA21309

# gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/long_read/winnowmap_alignments/HG005vGRCh38_wm_ONT.sort.bam
# gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/long_read/winnowmap_alignments/HG002vGRCh38_wm_ONT.sort.bam

