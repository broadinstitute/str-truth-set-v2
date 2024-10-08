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


DOCKER_IMAGE = "weisburd/dipcall-pipeline@sha256:6e13af8c4008fbea3ca1498e22c00df7265988f2a3673532f0b97f24c85c0dcc"

bp = pipeline("dipcall pipeline", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

parser = bp.get_config_arg_parser()
parser.add_argument("-s", "--sample-id", action="append",
                    help="Process only this sample. Can be specified more than once.")
parser.add_argument("--more-memory", action="store_true", help="Run with 2x more memory")
parser.add_argument("--sample-table", default="hprc_assemblies.tsv")
parser.add_argument("--urls-table", default="hprc_assembly_urls.tsv")
parser.add_argument("--output-dir", default="gs://str-truth-set-v2/dipcall_pipeline")
args = bp.parse_known_args()

df = pd.read_table(args.sample_table)
df_urls = pd.read_table(args.urls_table)
accession_to_url_map = dict(zip(df_urls.accession, df_urls.url))
df["url_pat"] = df["accession_pat"].map(accession_to_url_map)
df["url_mat"] = df["accession_mat"].map(accession_to_url_map)

if args.sample_id:
    df = df[df.sample_id.isin(args.sample_id)]

s1_steps = []
for i, (_, row) in enumerate(df.iterrows()):

    s1 = bp.new_step(
        f"dipcall: {row.sample_id}",
        image=DOCKER_IMAGE,
        arg_suffix="step1",
        cpu=16 if args.more_memory else 8,
        memory="highmem",
        storage="50Gi",
        output_dir=os.path.join(args.output_dir, row.sample_id))

    s1_steps.append(s1)

    hg38_fasta_input, _ = s1.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.COPY)

    s1.command("set -exuo pipefail")
    s1.command(f"cd /io")
    s1.command(f"wget --quiet {row.url_pat} & wget --quiet {row.url_mat} & wait")

    pseudoautosomal_region_arg = "-x /dipcall.kit/hs38.PAR.bed" if row["sex"] == "male" else ""
    s1.command(f"/dipcall.kit/run-dipcall "
               f"{pseudoautosomal_region_arg} "
               f"{row.sample_id} "
               f"{hg38_fasta_input} "
               f"{os.path.basename(row.url_pat)} "
               f"{os.path.basename(row.url_mat)} > out.mak")

    s1.command(f"make -j2 -f out.mak")

    s1.command("ls -lhtr")
    s1.command(f"[ -s {row.sample_id}.dip.bed ] || exit 1")  # check that the bed file isn't emtpy

    s1.command(f"bgzip {row.sample_id}.dip.bed")
    s1.command(f"tabix {row.sample_id}.dip.bed.gz")

    s1.output(f"{row.sample_id}.dip.vcf.gz")
    s1.output(f"{row.sample_id}.dip.bed.gz")
    s1.output(f"{row.sample_id}.dip.bed.gz.tbi")

    #s1.output(f"{row.sample_id}.hap1.bed")
    #s1.output(f"{row.sample_id}.hap2.bed")
    #s1.output(f"{row.sample_id}.pair.vcf.gz")

# combine high-confidence regions  (intersection, union)
s2 = bp.new_step(
    f"Combine high-confidence regions",
    image=DOCKER_IMAGE,
    arg_suffix="step2",
    cpu=1,
    output_dir=args.output_dir)

s2.depends_on(s1_steps)

s2.command("set -exuo pipefail")
for _, row in df.iterrows():
    current_input = s2.input(os.path.join(args.output_dir, row.sample_id, f"{row.sample_id}.dip.bed.gz"))
    s2.command(f"zcat {current_input} >> all.dip.bed")
output_bed_filename = f"combined.high_confidence_regions.{len(df)}_samples.union.bed.gz"
s2.command(f"bedtools sort -i all.dip.bed | bedtools merge -i - | bgzip > {output_bed_filename}")
s2.command(f"tabix {output_bed_filename}")
s2.output(output_bed_filename)
s2.output(f"{output_bed_filename}.tbi")

bp.run()


#%%
