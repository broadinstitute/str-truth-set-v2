"""Concatenates input TSVs into a single combined TSV"""

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description='Concatenate per-sample tables')
parser.add_argument('-o', '--output-tsv', help='Combined tsv file path')
parser.add_argument("input_tsvs", nargs='+', help="Input tsv files")
args = parser.parse_args()

if not args.output_tsv:
    args.output_tsv = f"concatenated.{len(args.input_tsvs)}_tables.tsv.gz"

if not args.output_tsv.endswith('.gz'):
    args.output_tsv += ".gz"

combined_df = None
for i, input_tsv in enumerate(args.input_tsvs):
    sample_id = os.path.basename(input_tsv).split('.')[0]

    df = pd.read_table(input_tsv)
    print(f"Read {sample_id:10s} table with {len(set(df.LocusId)):10,d} unique loci")

    df["SampleId"] = sample_id
    if combined_df is None:
        combined_df = df
    else:
        combined_df = pd.concat([combined_df, df])

print(f"Writing combined table with {len(combined_df):,d} rows and {len(set(combined_df.LocusId)):,d} unique loci to",
      args.output_tsv)

combined_df.to_csv(args.output_tsv, sep='\t', index=False)
