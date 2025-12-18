set -ex

python3 -u run_dipcall_pipeline.py \
    --sample-table hprc_assemblies_release2.tsv \
    --urls-table hprc_assemblies_release2_urls.tsv \
    --output-dir gs://str-truth-set-v2/dipcall_pipeline/HPRC_release2  --more-memory
