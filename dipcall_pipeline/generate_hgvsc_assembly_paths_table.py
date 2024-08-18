import collections
import pandas as pd
import re
import subprocess

metadata_df = pd.read_table("20130606_sample_info_1kGP.tsv")
sample_id_sex = dict(zip(metadata_df["Sample"], metadata_df["Gender"]))
sample_id_sex["NA24385"] = "male"

#%%
# print ftp:// paths of assebmly fastas from the HGVS Consortium
urls = subprocess.check_output("""for x in \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200828_JAX_assembly-results_CLR_v12/assemblies/phased/ \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200810_HHU_assembly-results_HG00514_v12/assemblies/phased/ \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200717_HHU_assembly-results_CCS_v12/assemblies/phased/  \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200716_UW_assembly-results_CLR_v12/assemblies/phased/ \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200628_HHU_assembly-results_CCS_v12/assemblies/phased/ \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200612_HHU_assembly-results_CLR_v12/assemblies/phased/
do
    for n in $(curl -sl $x )
    do
      echo ${x}/${n}
    done
done
""", shell=True, text=True).strip().split("\n")

#%%

url_table_rows = []
urls_by_sample_id = collections.defaultdict(dict)
for url in urls:
    filename = os.path.basename(url)
    if filename.endswith("fai"):
        continue
    match = re.match("v12_([^_]+)_(hgsvc|giab|hpg)_([^_]+)_1000-([^.]+).(h[12])", filename)
    if not match:
        raise ValueError(f"Unable to parse filename: {filename}")

    print(f"sample_id: {sample_id}, source_project: {source_project} sequencing_type: {sequencing_type}, assembler: {assembler}, h1_or_h2: {h1_or_h2}")

    row = {
        "sample_id": match.group(1),
        "source_project": match.group(2),
        "sequencing_type": match.group(3),
        "assembler": match.group(4),
        "h1_or_h2": match.group(5),
        "url": url,
    }
    url_table_rows.append({
        "accession": row["sample_id"] + "_" + row["h1_or_h2"],
        "url": row["url"],
    })

    urls_by_sample_id[row["sample_id"]][row["h1_or_h2"]] = row["url"]

pd.DataFrame(url_table_rows).to_csv("hgvsc_assembly_urls.tsv", sep="\t", index=False)

metadata_table_rows = []
for sample_id, sample_urls in urls_by_sample_id.items():
    if "h1" not in sample_urls or "h2" not in sample_urls:
        print(f"WARNING: missing h1 or h2 for sample_id: {sample_id}")
        continue
    h1_url = sample_urls["h1"]
    h2_url = sample_urls["h2"]
    row = {
        "sample_id": sample_id,
        "assembly_mat": row["sample_id"] + "_h1",
        "assembly_pat": row["sample_id"] + "_h2",
        "population": None,
        "accession_mat": row["sample_id"] + "_h1",
        "accession_pat": row["sample_id"] + "_h2",
        "sex": sample_id_sex.get(sample_id, None),
    }
    print(row)
    metadata_table_rows.append(row)

    if row["sex"] is None:
        print(f"WARNING: missing sex for sample_id: {sample_id}")

pd.DataFrame(metadata_table_rows).to_csv("hgvsc_assemblies.tsv", sep="\t", index=False)


#

#%%