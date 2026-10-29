import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Combine the dipcall assembly tables into the metadata table for the "
                                 "filter step. An individual in more than one table gets the first table's assemblies.")
parser.add_argument("--prefer-hprc-release2-and-hgsvc3-assemblies", action="store_true",
                    help="List HPRC release 2 first and HGSVC3 second, so that the individuals in HPRC release 1 or "
                         "HGSVC2 (all_assemblies.tsv) get their newer assemblies when they have them. The CHM13 truth "
                         "set uses this. The hg38 truth set doesn't, so those individuals keep their HPRC release 1 "
                         "and HGSVC2 assemblies there.")
args = parser.parse_args()

if args.prefer_hprc_release2_and_hgsvc3_assemblies:
    tables = [
        ("../dipcall_pipeline/hprc_assemblies_release2.tsv", "HPRC_release2"),
        ("../dipcall_pipeline/hgsvc3_assemblies.tsv", "HGSVC3"),
        ("../dipcall_pipeline/all_assemblies.tsv", ""),
        ("../dipcall_pipeline/human579_assemblies_excluding_hprc.tsv", "human579_assemblies"),
    ]
else:
    tables = [
        ("../dipcall_pipeline/all_assemblies.tsv", ""),
        ("../dipcall_pipeline/hprc_assemblies_release2.tsv", "HPRC_release2"),
        ("../dipcall_pipeline/human579_assemblies_excluding_hprc.tsv", "human579_assemblies"),
    ]

df_combined = None
for path, subdirectory in tables:
    df = pd.read_table(path)
    df["subdirectory"] = subdirectory
    print(f"Parsed {len(df):,d} samples from {path}")
    if df_combined is None:
        df_combined = df
        continue

    before = len(df)
    df = df[~df.sample_id.isin(df_combined.sample_id)]
    if before != len(df):
        print(f"Kept {len(df)} out of {before} rows")
    df_combined = pd.concat([df_combined, df])

output_path = f"all_assemblies_v2.{len(df_combined)}_samples" + (
    ".preferring_hprc_release2_and_hgsvc3.tsv" if args.prefer_hprc_release2_and_hgsvc3_assemblies else ".tsv")
df_combined.to_csv(output_path, header=True, sep="\t", index=False)
print(f"Wrote {len(df_combined):,d} rows to {output_path}")
