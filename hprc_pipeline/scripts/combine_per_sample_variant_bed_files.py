"""Performs an outer-join on per-sample tables and outputs a single combined table with one record per locus"""

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description="Perform an outer-join on multiple per-sample bed files")
parser.add_argument("-o", "--output-bed", help="Combined output bed file path")
parser.add_argument("--output-stats-tsv", action="store", const="+", default="-", nargs="?", help="If specified, will "
                    "output a table with stats. The optional value can be the path of this output file")
parser.add_argument("input_beds", nargs="+", help="Input bed files")
args = parser.parse_args()

if not args.output_bed:
    args.output_bed = f"combined.{len(args.input_beds)}_bed_files.bed"

if args.output_bed.endswith(".gz"):
    args.output_bed = args.output_bed.replace(".gz", "")

if args.output_stats_tsv == "-":
    args.output_stats_tsv = None
elif args.output_stats_tsv == "+":
    args.output_stats_tsv = f"combined.{len(args.input_beds)}_bed_files.stats.tsv"


for input_bed in args.input_beds:
    if not os.path.exists(input_bed):
        parser.error(f"Input file {input_bed} does not exist")

args.input_beds.sort(key=os.path.getsize, reverse=True)  # largest to smallest

def get_motif(name_field):
    tokens = name_field.split(":")
    if len(tokens) < 2:
        print("WARNING: Unexpected name field format:", name_field)
        return tokens[0]
    return tokens[1]

combined_df = None
output_stats = []
for input_i, input_bed in enumerate(args.input_beds):
    sample_id = os.path.basename(input_bed).split(".")[0]

    df = pd.read_table(input_bed, header=None, names=["chrom", "start_0based", "end", "motif", "score"])

    df["motif"] = df["motif"].apply(get_motif)
    df.set_index(["chrom", "start_0based", "end", "motif", "score"], inplace=True)

    if combined_df is None:
        locus_ids_before_join = 0
        combined_df = df
    else:
        locus_ids_before_join = len(combined_df)
        combined_df = combined_df.join(df, how="outer", rsuffix=f":{sample_id}")
        combined_df = combined_df.copy()  # intended to avoid the "DataFrame is highly fragmented" warning.

    duplicate_locus_id_count = len(combined_df) - len(set(combined_df.index))
    if duplicate_locus_id_count > 0:
        raise ValueError(f"{duplicate_locus_id_count:,d} duplicate locus ids found after adding table #{input_i+1}: {input_bed}")

    new_locus_id_count = len(combined_df) - locus_ids_before_join
    print(f"#{input_i+1}: Added {sample_id:10s} with {len(df):8,d} loci which yielded {new_locus_id_count:8,d} new loci"
          f" ({new_locus_id_count/len(combined_df):6.1%}) for an overall total of {len(combined_df):10,d} loci in the "
          f"combined table.")
    if args.output_stats_tsv:
        output_stats.append({
            "Id": input_i + 1,
            "SampleId": sample_id,
            "Loci": len(df),
            "NewLoci": new_locus_id_count,
            "FractionNewLoci": new_locus_id_count/len(combined_df),
            "CumulativeTotalLoci": len(combined_df),
        })

if args.output_stats_tsv:
    pd.DataFrame(output_stats).to_csv(args.output_stats_tsv, sep="\t", index=False)
    print(f"Wrote {len(output_stats):,d} rows to {args.output_stats_tsv}")

combined_df = combined_df.reset_index()
combined_df.sort_values(["chrom", "start_0based", "end"], inplace=True, ascending=True)
combined_df.to_csv(args.output_bed, sep="\t", index=False, header=False)
os.system(f"bgzip -f {args.output_bed}")
os.system(f"tabix -f {args.output_bed}.gz")
print(f"Wrote {len(combined_df):,d} loci to {args.output_bed}.gz")
