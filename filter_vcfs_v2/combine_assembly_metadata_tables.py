import pandas as pd

df_combined = None
for path, subdirectory in [
    ("../dipcall_pipeline/all_assemblies.tsv", ""),
    ("../dipcall_pipeline/hprc_assemblies_release2.tsv", "HPRC_release2"),
    ("../dipcall_pipeline/human579_assemblies_excluding_hprc.tsv", "human579_assemblies"),
]:
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

output_path = f"all_assemblies_v2.{len(df_combined)}_samples.tsv"
df_combined.to_csv(output_path, header=True, sep="\t", index=False)
print(f"Wrote {len(df_combined):,d} rows to {output_path}")
