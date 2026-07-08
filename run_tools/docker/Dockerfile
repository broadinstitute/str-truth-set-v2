# Image for standalone scripts that the run_tools pipeline(s) (run_tools/run_genotyping_tools.py) run as
# Hail Batch steps. Based on the str-analysis toolchain image, which already provides python, pandas, gcloud
# (for GSUTIL_COPY localize/delocalize), and the installed str_analysis package the scripts import. The repo's
# run_tools scripts are copied in so steps invoke them directly instead of localizing a loose file at runtime.
#
# Built from the repo root so the COPY below can see run_tools/ (see Makefile / build_run_tools_image.yml).
FROM weisburd/str-analysis@sha256:63260a1a684a929abfc0da525279fe43211bc49cd36bf5c1c90354af79b83ad4

# builds per-tool repeat catalogs (EHv5/GangSTR/HipSTR/TRGT/LongTR + inquiSTR bed) from a
# filter_vcf_to_tandem_repeats {sample}.tandem_repeat_genotypes.tsv.gz truth set
COPY run_tools/convert_truth_set_to_variant_catalogs.py /run_tools/convert_truth_set_to_variant_catalogs.py

ENV PYTHONPATH=/:.:$PYTHONPATH
ENTRYPOINT ["/bin/bash"]
