set -ex

python3 -u run_dipcall_pipeline.py \
    --sample-table hprc_assemblies.tsv \
    --urls-table hprc_assembly_urls.tsv \
    --output-dir gs://str-truth-set-v2/dipcall_pipeline  # --more-memory
