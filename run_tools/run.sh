set -ex

# Targeted validation run (HG002, homopolymers included):
#   pacbio 30x: TRGT + inquiSTR + vamos
#   ONT 26x:    inquiSTR + vamos  (TRGT is skipped on ONT since it only supports pacbio)
# The --filename-keyword values select the 30x pacbio and the 26x ONT bams.
python3 run_genotyping_tools.py \
  --data-type pacbio --data-type ONT \
  --tool TRGT --tool inquiSTR --tool vamos \
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
  --data-type pacbio --data-type ONT --data-type illumina --data-type illumina_exome --data-type element \
  --tool EHv5 --tool EHv5-bw2-optimized --tool GangSTR --tool HipSTR --tool TRGT --tool LongTR --tool inquiSTR --tool vamos \
  --sample-id HG002 --sample-id CHM1_CHM13
