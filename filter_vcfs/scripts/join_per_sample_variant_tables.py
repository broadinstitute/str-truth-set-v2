"""Performs an outer-join on per-sample tables and outputs a single combined table with one record per locus"""

import argparse
import os
import pandas as pd
import tqdm

parser = argparse.ArgumentParser(description="Perform an outer-join on multiple per-sample tables")
parser.add_argument("-o", "--output-tsv", help="Combined tsv file path")
parser.add_argument("-n", type=int, help="Number of tables to process")
parser.add_argument("--discard-impure-genotypes", action="store_true", help="Discard genotypes that are not pure repeats")
parser.add_argument("--output-stats-tsv", help="If specified, output a table with per-sample join stats to this path.")
parser.add_argument("input_tsvs", nargs="+", help="Input tsv files")
args = parser.parse_args()

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

# Define dtypes for each column to avoid mixed types warning
DTYPES = {
    "Chrom": "string",
    "Start1Based": "Int32",
    "End1Based": "Int32", 
    "Locus": "string",
    "LocusId": "string",
    "Motif": "string",
    "CanonicalMotif": "string",
    "MotifSize": "Int32",
    "NumRepeatsInReference": "Float32",
    "IsFoundInReference": "boolean",
    "NumRepeatsShortAllele": "Float32",
    "NumRepeatsLongAllele": "Float32",
    "IsPureRepeat": "boolean",
}

for input_tsv in args.input_tsvs:
    if not os.path.exists(input_tsv):
        parser.error(f"Input file {input_tsv} does not exist")

# Sort by file size (largest first)
args.input_tsvs.sort(key=os.path.getsize, reverse=True)

if args.n:
    args.input_tsvs = args.input_tsvs[:args.n]


if not args.output_tsv:
    args.output_tsv = f"joined.{len(args.input_tsvs)}_tables.tsv.gz"

if not args.output_tsv.endswith(".gz"):
    args.output_tsv += ".gz"


combined_df = None
output_stats = []
all_allele_size_columns = []


for table_i, input_tsv in tqdm.tqdm(enumerate(args.input_tsvs), total=len(args.input_tsvs), unit=" tables"):

    sample_id = os.path.basename(input_tsv).split(".")[0]

    df = pd.read_table(input_tsv, usecols=PER_LOCUS_COLUMNS + SAMPLE_SPECIFIC_COLUMNS, low_memory=False, dtype=DTYPES)
    missing_columns = set(PER_LOCUS_COLUMNS + SAMPLE_SPECIFIC_COLUMNS) - set(df.columns)
    if len(missing_columns) > 0:
        raise ValueError(f"{input_tsv} is missing these columns: {missing_columns}. Its columns are: {df.columns}")

    if args.discard_impure_genotypes:
        df = df[df["IsPureRepeat"]]

    df.set_index(PER_LOCUS_COLUMNS, inplace=True)

    # Rename columns efficiently - build rename dict once
    rename_dict = {}
    allele_columns = []
    for column in SAMPLE_SPECIFIC_COLUMNS:
        renamed_column = f"{column}:{sample_id}"
        rename_dict[column] = renamed_column
        if column.startswith("NumRepeats") and column.endswith("Allele"):
            allele_columns.append(renamed_column)
    
    df.rename(columns=rename_dict, inplace=True)
    all_allele_size_columns.extend(allele_columns)

    if combined_df is None:
        locus_ids_before_join = 0
        combined_df = df
    else:
        locus_ids_before_join = len(combined_df)
        # Use concat instead of join for better performance with large datasets
        combined_df = pd.concat([combined_df, df], axis=1, join='outer')

    #if table_i % 50 == 0:
    #    combined_df = combined_df.copy()  # intended to avoid the "DataFrame is highly fragmented" warning.

    # Check for duplicates more efficiently
    if combined_df.index.duplicated().any():
        duplicate_count = combined_df.index.duplicated().sum()
        raise ValueError(f"{duplicate_count:,d} duplicate locus ids found after adding table #{table_i+1}: {input_tsv}")

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

for c in all_allele_size_columns:
    num_empty_values = combined_df[c].isna().sum()
    combined_df[c] = combined_df[c].fillna(combined_df["NumRepeatsInReference"])
    print(f"Filled {num_empty_values:,d} empty values in column {c} out of {len(combined_df):,d} total rows")

combined_df.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote combined table with {len(combined_df):,d} loci to {args.output_tsv}")
