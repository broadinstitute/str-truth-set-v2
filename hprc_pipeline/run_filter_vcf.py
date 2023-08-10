"""Hail Batch pipeline for running dipcall on HPRC assemblies.

Relevant links:

HPRC assemblies
  https://projects.ensembl.org/hprc/

The design and construction of reference pangenome graphs with minigraph by Li et al. 2020
  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z

Increased mutation and gene conversion within human segmental duplications by Vollger et al.
  https://www.nature.com/articles/s41586-023-05895-y

"""
import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize


DOCKER_IMAGE = "weisburd/str-analysis@sha256:aac8013e830e48e2033d6933ccf4bc8556cb141a918bcf79f806c91c6e99a0d2"

bp = pipeline("filter_vcf_to_STRs", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

parser = bp.get_config_arg_parser()
parser.add_argument("--exclude-homopolymers", action="store_true")
parser.add_argument("-s", "--sample-id", action="append",
                    help="Process only this sample. Can be specified more than once.")
parser.add_argument("--output-dir", default="gs://str-truth-set-v2/filter_vcf")
args = bp.parse_known_args()

df = pd.read_table("hprc_assemblies.tsv")
if args.sample_id:
    df = df[df.sample_id.isin(args.sample_id)]

suffix = ".STRs"
if args.exclude_homopolymers:
    suffix += "_excluding_homopolymers"

for i, (_, row) in enumerate(df.iterrows()):

    s1 = bp.new_step(
        f"filter_vcf: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="step1",
        cpu=1,
        memory="highmem",
        #storage="50Gi",
        output_dir=os.path.join(args.output_dir, row.sample_id))

    hg38_fasta_input, _ = s1.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.COPY)

    dipcall_vcf_input, dipcall_high_confidence_regions_bed_input = s1.inputs(
        f"gs://str-truth-set-v2/hprc_dipcall/{row.sample_id}/{row.sample_id}.dip.vcf.gz",
        f"gs://str-truth-set-v2/hprc_dipcall/{row.sample_id}/{row.sample_id}.dip.bed",
    )

    s1.command("set -exuo pipefail")

    s1.command(f"[ -s {dipcall_high_confidence_regions_bed_input} ] || exit 1")  # check that the bed file isn't emtpy

    s1.command(f"bedtools intersect -header -f 1 -wa -u \
        -a {dipcall_vcf_input}  \
        -b {dipcall_high_confidence_regions_bed_input} \
        | bgzip > {row.sample_id}.high_confidence_regions.vcf.gz")
    s1.command(f"tabix -f {row.sample_id}.high_confidence_regions.vcf.gz")

    min_repeat_unit_length = 2 if args.exclude_homopolymers else 1
    s1.command(f"python3 -u -m str_analysis.filter_vcf_to_STR_variants \
        -R {hg38_fasta_input} \
        --allow-interruptions only-if-pure-repeats-not-found \
        --write-bed-file \
        --min-str-length 9 \
        --min-str-repeats 3 \
        --min-repeat-unit-length {min_repeat_unit_length} \
        --max-repeat-unit-length 50 \
        --output-prefix {row.sample_id}{suffix} \
        --verbose \
        {row.sample_id}.high_confidence_regions.vcf.gz |& tee {row.sample_id}{suffix}.filter_vcf.log")

    s1.command("ls -lhtr")

    s1.output(f"{row.sample_id}.high_confidence_regions.vcf.gz")
    s1.output(f"{row.sample_id}.high_confidence_regions.vcf.gz.tbi")
    s1.output(f"{row.sample_id}{suffix}.vcf.gz")
    s1.output(f"{row.sample_id}{suffix}.vcf.gz.tbi")
    s1.output(f"{row.sample_id}{suffix}.variants.tsv.gz")
    s1.output(f"{row.sample_id}{suffix}.alleles.tsv.gz")
    s1.output(f"{row.sample_id}{suffix}.variants.bed.gz")
    s1.output(f"{row.sample_id}{suffix}.variants.bed.gz.tbi")
    s1.output(f"{row.sample_id}{suffix}.filter_vcf.log")

bp.run()


#%%
