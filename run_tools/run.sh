set -ex

# Targeted validation run (HG002, homopolymers included):
#   pacbio 30x: TRGTv3 + TRGTv5 + inquiSTR + vamos
#   ONT 26x:    TRGTv3 + TRGTv5 + inquiSTR + vamos
# The --filename-keyword values select the 30x pacbio and the 26x ONT bams.
python3 run_genotyping_tools.py \
  --data-type pacbio --data-type ONT \
  --tool TRGTv3 --tool TRGTv5 --tool inquiSTR --tool vamos \
  --sample-id HG002 \
  --filename-keyword downsampled_to_30x.bam --filename-keyword ONT/HG002.bam

# Illumina short-read validation run (HG002): both ExpansionHunter v5 variants
#   EHv5               = bw2 fork, --analysis-mode low-mem-streaming
#   EHv5-bw2-optimized = bw2 fork, --analysis-mode optimized-streaming --improved-genotyping
python3 run_genotyping_tools.py \
  --data-type illumina \
  --tool EHv5 --tool EHv5-bw2-optimized \
  --sample-id HG002

exit 0

# Full run (all tools, all data types), homopolymers included:
python3 run_genotyping_tools.py \
  --data-type pacbio --data-type ONT --data-type illumina --data-type illumina_exome --data-type element --data-type ultima --data-type illumina_rnaseq --data-type pacbio_isoseq \
  --tool EHv5 --tool EHv5-bw2-optimized --tool IlluminaEHv5 --tool GangSTR --tool HipSTR --tool TRGTv3 --tool TRGTv5 --tool LongTR --tool inquiSTR --tool vamos \
  --sample-id HG002 --sample-id CHM1_CHM13
