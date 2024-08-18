import os
import pandas as pd

if not os.getcwd().endswith("dipcall_pipeline"):
	os.chdir("dipcall_pipeline")

expected_columns = {'sample_id', 'assembly_mat', 'assembly_pat', 'population', 'accession_mat', 'accession_pat', 'sex'}
dfs = []
added_sample_ids = set()
for filename in ["hprc_assemblies.tsv", "hgvsc_assemblies.tsv"]:
	df = pd.read_table(filename)
	if set(df.columns) != expected_columns:
		raise ValueError(f"ERROR in {filename}:\n{expected_columns}\ncolumns were expected\n{set(df.columns)}\nfound instead")
	previously_added_sample_ids = added_sample_ids & set(df.sample_id)
	if previously_added_sample_ids:
		print(f"NOTE: Skipping {len(previously_added_sample_ids)} samples from {filename} that are already in the combined table:",
			  ", ".join(previously_added_sample_ids))
		df = df[~df.sample_id.isin(previously_added_sample_ids)]
	added_sample_ids.update(df.sample_id)

	print(f"Adding {len(df)} samples from {filename}")
	dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

for column in expected_columns:
	missing_values = sum(df[column].isna())
	if missing_values > 0 and column != "population":
		raise ValueError(f"ERROR: '{column}' column is missing for {missing_values} samples")

output_filename = "all_assemblies.tsv"
df.to_csv(output_filename, sep="\t", index=False)
print(f"Wrote {len(df)} samples to {output_filename}")


#%%
