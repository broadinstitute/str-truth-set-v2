"""Hail Batch pipeline for running dipcall on HPRC assemblies.

Relevant links:

HPRC assemblies
  https://projects.ensembl.org/hprc/

The design and construction of reference pangenome graphs with minigraph by Li et al. 2020
  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z

Increased mutation and gene conversion within human segmental duplications by Vollger et al.
  https://www.nature.com/articles/s41586-023-05895-y

A collection of 579 high-quality human assemblies (including HPRC release 1 and 2) is a available at:
https://github.com/lh3/OpenHGL
https://zenodo.org/records/17772289
"""

import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize


DOCKER_IMAGE = "weisburd/dipcall-pipeline@sha256:bf74b72c4b0a44932f5d659350153526b3a56015463ea63d0121371b769dc76e"

# dipcall requires the chrY pseudoautosomal regions (PARs) to be hard-masked in the reference, and takes a BED of the
# chrX PARs for male samples (see the dipcall README). hg38.fa's chrY PARs are masked, and the dipcall kit in the
# docker image includes hs38.PAR.bed. For CHM13, Broad's public copy of chm13v2.0.maskedY.fasta has its chrY PARs
# masked, unlike gs://str-truth-set/chm13/ref/chm13v2.0.fa.
REFERENCE_GENOMES = {
    "hg38": {
        "fasta": "gs://str-truth-set/hg38/ref/hg38.fa",
        "pseudoautosomal_regions_bed": "/dipcall.kit/hs38.PAR.bed",
        "output_dir": "gs://str-truth-set-v2/dipcall_pipeline",
    },
    "chm13": {
        "fasta": "gs://gcp-public-data--broad-references/t2t/v2/chm13v2.0.maskedY.fasta",
        "pseudoautosomal_regions_bed": "gs://gcp-public-data--broad-references/t2t/v2/chm13v2.0_PAR.bed",
        "output_dir": "gs://str-truth-set-v2/chm13/dipcall_pipeline",
    },
}

bp = pipeline("dipcall pipeline", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

parser = bp.get_config_arg_parser()
parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample. Can be specified more than once.")
parser.add_argument("-n", "--num-samples", type=int, help="Process only this number of samples.")
parser.add_argument("--offset", type=int, help="Skip the first N samples.")
parser.add_argument("--more-memory", action="store_true", help="Run with 2x more memory")
parser.add_argument("--sample-table", default="hprc_assemblies.tsv", help="Sample table path")
g = parser.add_mutually_exclusive_group(required=True)
g.add_argument("--urls-table", help="URLs table path")   # default="hprc_assembly_urls.tsv"
g.add_argument("--agc-archive-url", help="Output table path")   # default="https://zenodo.org/records/17772289/files/human579.agc" or "gs://str-truth-set/hg38/ref/human579.agc"
parser.add_argument("--reference", choices=REFERENCE_GENOMES.keys(), default="hg38",
                    help="Reference genome to align the assemblies to")
parser.add_argument("--output-dir", help="Output bucket path. Defaults to " + ", ".join(
    f"{reference} => {files['output_dir']}" for reference, files in REFERENCE_GENOMES.items()))
args = bp.parse_known_args()

reference_genome = REFERENCE_GENOMES[args.reference]
if not args.output_dir:
    args.output_dir = reference_genome["output_dir"]
elif (args.reference != "hg38" and args.reference not in args.output_dir) or any(
        reference in args.output_dir for reference in REFERENCE_GENOMES if reference != args.reference):
    # step_pipeline skips a job whose outputs already exist, so pointing a run at a folder of another reference's
    # results would silently report those results as this reference's.
    parser.error(f"--output-dir {args.output_dir} doesn't match --reference {args.reference}. To keep each "
                 f"reference's results apart, a non-hg38 folder must include the reference's name, and no folder may "
                 f"include another reference's name.")

df = pd.read_table(args.sample_table)
if args.urls_table:
    df_urls = pd.read_table(args.urls_table)
    accession_to_url_map = dict(zip(df_urls.accession, df_urls.url))
    df["url_pat"] = df["accession_pat"].map(accession_to_url_map)
    df["url_mat"] = df["accession_mat"].map(accession_to_url_map)

if args.sample_id:
    df = df[df["sample_id"].isin(args.sample_id)]

if args.offset:
    df = df.iloc[args.offset:]

if args.num_samples:
    df = df.iloc[:args.num_samples]

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

    reference_fasta_input, _ = s1.inputs(
        reference_genome["fasta"],
        f"{reference_genome['fasta']}.fai",
        localize_by=Localize.COPY)

    s1.command("set -exuo pipefail")
    s1.command(f"cd /io")
    if args.urls_table:
        s1.command(f"wget --quiet {row.url_pat} & wget --quiet {row.url_mat} & wait")
        father_fasta_path = os.path.basename(row.url_pat)
        mother_fasta_path = os.path.basename(row.url_mat)
    elif args.agc_archive_url:
        agc_archive_path = s1.input(args.agc_archive_url)
        father_fasta_path = f"{row.assembly_pat}.fasta"
        mother_fasta_path = f"{row.assembly_mat}.fasta"
        s1.command(f"agc getset {agc_archive_path} {row.assembly_pat} > {father_fasta_path} &")
        s1.command(f"agc getset {agc_archive_path} {row.assembly_mat} > {mother_fasta_path} &")
        s1.command(f"wait")

    else:
        parser.error("Must specify either --urls-table or --agc-archive-url")

    pseudoautosomal_region_arg = ""
    if row.sex == "male":
        pseudoautosomal_regions_bed = reference_genome["pseudoautosomal_regions_bed"]
        if pseudoautosomal_regions_bed.startswith("gs://"):
            pseudoautosomal_regions_bed = s1.input(pseudoautosomal_regions_bed, localize_by=Localize.COPY)
        pseudoautosomal_region_arg = f"-x {pseudoautosomal_regions_bed}"
    s1.command(f"/dipcall.kit/run-dipcall "
               f"{pseudoautosomal_region_arg} "
               f"{row.sample_id} "
               f"{reference_fasta_input} "
               f"{father_fasta_path} "
               f"{mother_fasta_path} > out.mak")

    s1.command(f"make -j2 -f out.mak")

    s1.command("ls -lhtr")
    s1.command(f"[ -s {row.sample_id}.dip.bed ] || exit 1")  # check that the bed file isn't emtpy

    s1.command(f"bgzip {row.sample_id}.dip.bed")
    s1.command(f"tabix {row.sample_id}.dip.bed.gz")

    s1.command(f"tabix {row.sample_id}.dip.vcf.gz")

    s1.output(f"{row.sample_id}.dip.vcf.gz")
    s1.output(f"{row.sample_id}.dip.vcf.gz.tbi")

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
