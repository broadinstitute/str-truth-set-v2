import os
import pandas as pd

# Links: https://projects.ensembl.org/hprc/
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z

from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "weisburd/hprc-pipeline@sha256:70360d3cb05a49afcf303ca36aa794462159e2e758142f81ec0f1e4bcaf60039"

bp = pipeline("HPRC dipcall pipeline", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline_gnomad")

parser = bp.get_config_arg_parser()
parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample id. Can be specified more than once.")
parser.add_argument("--output-dir", default="gs://str-truth-set-v2/hprc_dipcall")
args = bp.parse_known_args()

df = pd.read_table("hprc_assemblies.tsv")
df_urls = pd.read_table("hprc_assembly_urls.tsv")
accession_to_url_map = dict(zip(df_urls.accession, df_urls.url))
df["url_pat"] = df["accession_pat"].map(accession_to_url_map)
df["url_mat"] = df["accession_mat"].map(accession_to_url_map)

if args.sample_id:
    df = df[df.sample_id.isin(args.sample_id)]

for i, (_, row) in enumerate(df.iterrows()):

    s1 = bp.new_step(
        f"HPRC dipcall: {row.sample_id}",
        image=DOCKER_IMAGE,
        cpu=8,
        memory="highmem",
        storage="50Gi",
        output_dir=os.path.join(args.output_dir, row.sample_id))

    hg38_fasta_input, _ = s1.inputs(
        "gs://str-truth-set/hg38/ref/hg38.fa",
        "gs://str-truth-set/hg38/ref/hg38.fa.fai",
        localize_by=Localize.COPY)

    s1.command("set -exuo pipefail")
    s1.command(f"cd /io")
    s1.command(f"wget --quiet {row.url_pat} & wget --quiet {row.url_mat} & wait")
    s1.command(f"/dipcall.kit/run-dipcall "
               f"{row.sample_id} "
               f"{hg38_fasta_input} "
               f"{os.path.basename(row.url_pat)} "
               f"{os.path.basename(row.url_mat)} > out.mak")

    s1.command(f"make -j2 -f out.mak")
    s1.command("ls -lhtr")

    s1.output(f"{row.sample_id}.dip.vcf.gz")
    s1.output(f"{row.sample_id}.dip.bed")
    #s1.output(f"{row.sample_id}.hap1.bed")
    #s1.output(f"{row.sample_id}.hap2.bed")
    #s1.output(f"{row.sample_id}.pair.vcf.gz")

bp.run()


#%%
