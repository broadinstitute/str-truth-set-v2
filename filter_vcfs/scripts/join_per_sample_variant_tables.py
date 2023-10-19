"""Performs an outer-join on per-sample tables and outputs a single combined table with one record per locus"""

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description="Perform an outer-join on multiple per-sample tables")
parser.add_argument("-o", "--output-tsv", help="Combined tsv file path")
parser.add_argument("--output-stats-tsv", action="store", const="+", default="-", nargs="?", help="If specified, will "
                    "output a table with stats. The optional value can be the path of this output file")
parser.add_argument("input_tsvs", nargs="+", help="Input tsv files")
args = parser.parse_args()

if not args.output_tsv:
    args.output_tsv = f"joined.{len(args.input_tsvs)}_tables.tsv.gz"

if not args.output_tsv.endswith(".gz"):
    args.output_tsv += ".gz"

if args.output_stats_tsv == "-":
    args.output_stats_tsv = None
elif args.output_stats_tsv == "+":
    args.output_stats_tsv = f"combined.{len(args.input_beds)}_bed_files.stats.tsv"

"""
$1                    Chrom : chr1
$2              Start1Based : 674824
$3                End1Based : 674832
$4                    Locus : 1:674824-674832
$5                  LocusId : 1-674823-674832-AGG
$6               INS_or_DEL : INS
$7                    Motif : AGG
$8   MotifInterruptionIndex : 2.0
$9           CanonicalMotif : AGG
$10               MotifSize : 3
$11   NumRepeatsInReference : 3.0
$12                  VcfPos : 674829
$13                  VcfRef : C
$14                  VcfAlt : CAGG
$15             VcfGenotype : 1|0
$16           SummaryString : 3bp:AGG:INS:3=>4:HET:not-pure
$17      IsFoundInReference : True
$18            IsPureRepeat : False
$19          IsMultiallelic : False
$20              NumRepeats : 4
$21         RepeatSize (bp) : 12
$22          NumPureRepeats : 3
$23     PureRepeatSize (bp) : 9
$24     FractionPureRepeats : 0.75
"""
PER_LOCUS_COLUMNS = [
    "Chrom",
    "Start1Based",
    "End1Based",
    "Locus",
    "LocusId",
    #"INS_or_DEL",
    "Motif",
    "CanonicalMotif",
    "MotifSize",
    "NumRepeatsInReference",
    "IsFoundInReference",
]
SAMPLE_SPECIFIC_COLUMNS = [
    #"VcfPos",
    #"VcfRef",
    #"VcfAlt",
    #"VcfGenotype",
    #"SummaryString",
    #"IsMultiallelic",
    #"HET_or_HOM_or_HEMI_or_MULTI",
    "NumRepeatsShortAllele",
    "NumRepeatsLongAllele",
    #"RepeatSizeShortAllele (bp)",
    #"RepeatSizeLongAllele (bp)",
    "IsPureRepeat",
    #"MotifInterruptionIndex",

    #"NumRepeats",
    #"RepeatSize (bp)",
    #"NumPureRepeats",
    #"PureRepeatSize (bp)",
    #"FractionPureRepeats",
]

for input_tsv in args.input_tsvs:
    if not os.path.exists(input_tsv):
        parser.error(f"Input file {input_tsv} does not exist")

args.input_tsvs.sort(key=os.path.getsize, reverse=True)  # largest to smallest

combined_df = None
output_stats = []
for table_i, input_tsv in enumerate(args.input_tsvs):
    sample_id = os.path.basename(input_tsv).split(".")[0]

    df = pd.read_table(input_tsv)
    missing_columns = set(PER_LOCUS_COLUMNS + SAMPLE_SPECIFIC_COLUMNS) - set(df.columns)
    if len(missing_columns) > 0:
        raise ValueError(f"{input_tsv} is missing these columns: {missing_columns}. Its columns are: {df.columns}")

    df = df[PER_LOCUS_COLUMNS + SAMPLE_SPECIFIC_COLUMNS]
    df.set_index(PER_LOCUS_COLUMNS, inplace=True)

    for column in SAMPLE_SPECIFIC_COLUMNS:
        df.rename(columns={
            column: f"{column}:{sample_id}",
        }, inplace=True)
    if combined_df is None:
        locus_ids_before_join = 0
        combined_df = df
    else:
        locus_ids_before_join = len(combined_df)
        combined_df = combined_df.join(df, how="outer", rsuffix=f":{sample_id}")
        combined_df = combined_df.copy()  # intended to avoid the "DataFrame is highly fragmented" warning.

    duplicate_locus_id_count = len(combined_df) - len(set(combined_df.index))
    if duplicate_locus_id_count > 0:
        raise ValueError(f"{duplicate_locus_id_count:,d} duplicate locus ids found after adding table #{table_i+1}: {input_tsv}")

    new_locus_id_count = len(combined_df) - locus_ids_before_join
    print(f"#{table_i+1}: Added {sample_id:10s} with {len(df):8,d} loci which yielded {new_locus_id_count:8,d} new loci"
          f" ({new_locus_id_count/len(combined_df):6.1%}) for an overall total of {len(combined_df):10,d} loci in the "
          f"combined table.")
    if args.output_stats_tsv:
        output_stats.append({
            "Id": table_i + 1,
            "SampleId": sample_id,
            "Loci": len(df),
            "NewLoci": new_locus_id_count,
            "FractionNewLoci": new_locus_id_count/len(combined_df),
            "CumulativeTotalLoci": len(combined_df),
        })

if args.output_stats_tsv:
    pd.DataFrame(output_stats).to_csv(args.output_stats_tsv, sep="\t", index=False)
    print(f"Wrote {len(output_stats):,d} rows to {args.output_stats_tsv}")

is_pure_repeat_columns = [c for c in combined_df.columns if c.startswith("IsPureRepeat:")]
combined_df["IsPureRepeat"] = combined_df[is_pure_repeat_columns].all(axis=1)

combined_df.drop(columns=is_pure_repeat_columns, inplace=True)

combined_df = combined_df.reset_index()
combined_df.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote combined table with {len(combined_df):,d} loci to {args.output_tsv}")
