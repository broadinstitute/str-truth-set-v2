"""Write hgsvc3_assemblies.tsv and hgsvc3_assembly_urls.tsv for the HGSVC3 Verkko assemblies of the HGSVC2 individuals.

The CHM13 truth set uses these instead of the HGSVC2 assemblies in hgvsc_assemblies.tsv. Each HGSVC2 haplotype
assembly (h1-un, h2-un) also includes the reads that couldn't be assigned to a haplotype, which in a male are all of
the chrX and chrY reads. So both haplotypes contain chrX and chrY, and dipcall, which expects chrX only in the maternal
and chrY only in the paternal assembly, calls almost nothing on them. HGSVC3's Strand-seq-phased Verkko assemblies put
a male's chrX and chrY in different haplotypes, but which one is hap1 varies by sample. So this script reads HGSVC3's
contig-to-CHM13 chromosome reports and, for males, lists the haplotype that carries chrY as the paternal assembly.
For females, hap1 is listed as paternal and hap2 as maternal, which doesn't matter to dipcall.

HGSVC3 data reuse policy: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/README_HGSVC3_datareuse_statement.md

Usage:
    python3 generate_hgsvc3_assembly_tables.py
"""

import io
import re
import subprocess
import pandas as pd

HGSVC3_WORKING_DIR = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working"

# HG00514 is in both 20240201_verkko_batch3 and 20241001_verkko_HG00514_fix. HGSVC3's Assembly_Info/v2.0 README says
# the batch3 version has corrupted data and only the fix should be used, so a later folder overrides an earlier one.
# 20230711_verkko_pilot is left out because all 5 of its samples were reassembled in 20230818_verkko_batch1.
HGSVC3_VERKKO_ASSEMBLY_FOLDERS = [
    "20230818_verkko_batch1",
    "20230927_verkko_batch2",
    "20240201_verkko_batch3",
    "20241001_verkko_HG00514_fix",
]

CHROMOSOME_LENGTHS = {"chrX": 154_259_566, "chrY": 62_460_029}  # in T2T-CHM13v2.0, which the reports align to

# A male's chrX haplotype must cover at least this much of chrX and its chrY haplotype this much of chrY, and neither
# may cover more than MAX_PERCENT_OF_OTHER_SEX_CHROMOSOME_IN_MALE_HAPLOTYPE of the other one. A female's haplotypes must
# each cover chrX, and neither may cover more than MAX_PERCENT_OF_CHRY_IN_FEMALE of chrY.
MIN_PERCENT_OF_CHRX_IN_CHRX_HAPLOTYPE = 90
MIN_PERCENT_OF_CHRY_IN_MALE_CHRY_HAPLOTYPE = 50
MAX_PERCENT_OF_OTHER_SEX_CHROMOSOME_IN_MALE_HAPLOTYPE = 10
MAX_PERCENT_OF_CHRY_IN_FEMALE = 1

# HGSVC3 couldn't assign one of HG00732's two chrX copies to a haplotype: hap1 has none of chrX, and 83.5% of chrX is
# in its unassigned contigs. So dipcall leaves HG00732's chrX out of its confident regions. It's used anyway, since its
# autosomes are much more contiguous than in its HGSVC2 assembly (contig N50 134 vs 22 Mb).
FEMALES_WITH_CHRX_IN_ONLY_ONE_HAPLOTYPE = {"HG00732"}


def run_curl(url):
    """Returns the contents of the given URL (for an ftp:// folder, the names of the files in it)."""
    return subprocess.run(["curl", "--silent", "--fail", "--list-only" if url.endswith("/") else "--globoff",
                           "--max-time", "300", url], check=True, capture_output=True, text=True).stdout


def get_hgsvc3_verkko_assembly_folder_of_each_sample():
    """Returns {sample_id: the HGSVC3 working-directory folder that holds its Verkko assembly}."""
    folder_of_each_sample = {}
    for folder in HGSVC3_VERKKO_ASSEMBLY_FOLDERS:
        for sample_id in run_curl(f"{HGSVC3_WORKING_DIR}/{folder}/assemblies/").split():
            folder_of_each_sample[sample_id] = folder
    return folder_of_each_sample


def find_haplotype_fasta_url(folder, sample_id, haplotype):
    """Returns the URL of the haplotype's FASTA, which the HG00514 fix names <sample>.vrk-ps-sseq.asm-<hap>.v2.fasta.gz.

    Args:
        folder (str): HGSVC3 working-directory folder that holds the sample's Verkko assembly
        sample_id (str): sample id, eg. HG00096
        haplotype (str): "hap1" or "hap2"

    Returns:
        str: ftp:// URL
    """
    assembly_dir = f"{HGSVC3_WORKING_DIR}/{folder}/assemblies/{sample_id}/"
    filenames = [filename for filename in run_curl(assembly_dir).split() if re.fullmatch(
        rf"{sample_id}\.vrk-ps-sseq\.asm-{haplotype}(\.v2)?\.fasta\.gz", filename)]
    if len(filenames) != 1:
        raise ValueError(f"Expected 1 {haplotype} FASTA in {assembly_dir}, found {len(filenames)}: {filenames}")
    return assembly_dir + filenames[0]


def compute_percent_of_chrx_and_chry_in_haplotype(folder, sample_id, haplotype):
    """Returns the percent of CHM13's chrX and chrY covered by the haplotype's contigs that mostly align there.

    Args:
        folder (str): HGSVC3 working-directory folder that holds the sample's Verkko assembly
        sample_id (str): sample id, eg. HG00096
        haplotype (str): "hap1" or "hap2"

    Returns:
        dict: {"chrX": percent, "chrY": percent}, or None if HGSVC3 has no chromosome report for the haplotype (as for
            the HG00514 fix, whose README says its reports "will be updated in a later upload").
    """
    report_url = (f"{HGSVC3_WORKING_DIR}/{folder}/reports/ref_chrom_assign/{sample_id}/"
                  f"{sample_id}.vrk-ps-sseq.asm-{haplotype}.t2tv2.chrom-assign-by-target.tsv")
    try:
        report_text = run_curl(report_url)
    except subprocess.CalledProcessError as e:
        if e.returncode in (9, 78):  # the FTP server's "access denied" and "file not found", which mean no such file
            return None
        raise
    report_df = pd.read_table(io.StringIO(report_text), comment="#")
    percents = {}
    for chrom, length in CHROMOSOME_LENGTHS.items():
        # skip short alignments of contigs that belong to another chromosome
        rows = report_df[(report_df.target_name == chrom) & (report_df.match_cov_query > 0.5)]
        percents[chrom] = 100 * rows.align_matching.sum() / length
    return percents


def choose_paternal_and_maternal_haplotypes(sample_id, sex, percents_by_haplotype,
                                            allow_chrx_in_only_one_haplotype=False):
    """Returns (paternal haplotype, maternal haplotype), checking that chrX and chrY are where the sex says they are.

    Args:
        sample_id (str): sample id, used in error messages
        sex (str): "male" or "female"
        percents_by_haplotype (dict): {"hap1": {"chrX": percent, "chrY": percent}, "hap2": {...}}
        allow_chrx_in_only_one_haplotype (bool): for a female, require chrX in just one of the haplotypes

    Returns:
        tuple: ("hap1", "hap2") or ("hap2", "hap1"). For a male, the first one carries chrY.
    """
    description = ", ".join(f"{haplotype} has {percents['chrX']:.1f}% of chrX and {percents['chrY']:.1f}% of chrY"
                            for haplotype, percents in percents_by_haplotype.items())
    if sex == "male":
        paternal = max(percents_by_haplotype, key=lambda haplotype: percents_by_haplotype[haplotype]["chrY"])
        maternal = "hap1" if paternal == "hap2" else "hap2"
        if (percents_by_haplotype[paternal]["chrY"] < MIN_PERCENT_OF_CHRY_IN_MALE_CHRY_HAPLOTYPE
                or percents_by_haplotype[maternal]["chrX"] < MIN_PERCENT_OF_CHRX_IN_CHRX_HAPLOTYPE
                or percents_by_haplotype[paternal]["chrX"] > MAX_PERCENT_OF_OTHER_SEX_CHROMOSOME_IN_MALE_HAPLOTYPE
                or percents_by_haplotype[maternal]["chrY"] > MAX_PERCENT_OF_OTHER_SEX_CHROMOSOME_IN_MALE_HAPLOTYPE):
            raise ValueError(f"{sample_id} is male, but its chrX and chrY aren't in separate haplotypes: {description}")
        return paternal, maternal

    if sex == "female":
        number_of_haplotypes_with_chrx = sum(percents["chrX"] >= MIN_PERCENT_OF_CHRX_IN_CHRX_HAPLOTYPE
                                             for percents in percents_by_haplotype.values())
        if (number_of_haplotypes_with_chrx < (1 if allow_chrx_in_only_one_haplotype else 2)
                or any(percents["chrY"] > MAX_PERCENT_OF_CHRY_IN_FEMALE for percents in percents_by_haplotype.values())):
            raise ValueError(f"{sample_id} is female, but its haplotypes don't each have one chrX and no chrY: "
                             f"{description}")
        return "hap1", "hap2"

    raise ValueError(f"{sample_id} has unexpected sex: {sex}")


def main():
    hgsvc2_df = pd.read_table("hgvsc_assemblies.tsv")
    folder_of_each_sample = get_hgsvc3_verkko_assembly_folder_of_each_sample()
    print(f"Found HGSVC3 Verkko assemblies of {len(folder_of_each_sample)} samples")
    missing_sample_ids = sorted(set(hgsvc2_df.sample_id) - set(folder_of_each_sample))
    print(f"{len(missing_sample_ids)} of the {len(hgsvc2_df)} HGSVC2 samples have no HGSVC3 Verkko assembly: "
          + ", ".join(missing_sample_ids))

    assembly_rows = []
    url_rows = []
    for _, row in hgsvc2_df[hgsvc2_df.sample_id.isin(folder_of_each_sample)].iterrows():
        folder = folder_of_each_sample[row.sample_id]
        percents_by_haplotype = {
            haplotype: compute_percent_of_chrx_and_chry_in_haplotype(folder, row.sample_id, haplotype)
            for haplotype in ("hap1", "hap2")
        }
        if None in percents_by_haplotype.values():
            if row.sex != "male":
                print(f"WARNING: {row.sample_id} ({row.sex}) has no HGSVC3 chromosome report, so its sex can't be "
                      f"checked. Using hap1 as paternal and hap2 as maternal.")
                paternal, maternal = "hap1", "hap2"
            else:
                raise ValueError(f"{row.sample_id} is male but has no HGSVC3 chromosome report, so its chrY "
                                 f"haplotype is unknown")
        else:
            paternal, maternal = choose_paternal_and_maternal_haplotypes(
                row.sample_id, row.sex, percents_by_haplotype,
                allow_chrx_in_only_one_haplotype=row.sample_id in FEMALES_WITH_CHRX_IN_ONLY_ONE_HAPLOTYPE)
            print(f"{row.sample_id} ({row.sex}): paternal = {paternal}, maternal = {maternal}. " + ", ".join(
                f"{haplotype}: chrX {percents['chrX']:.1f}%, chrY {percents['chrY']:.1f}%"
                for haplotype, percents in percents_by_haplotype.items()))

        accessions = {haplotype: f"{row.sample_id}_hgsvc3_{haplotype}" for haplotype in ("hap1", "hap2")}
        assembly_rows.append({
            "sample_id": row.sample_id,
            "assembly_mat": accessions[maternal],
            "assembly_pat": accessions[paternal],
            "population": row.population,
            "accession_mat": accessions[maternal],
            "accession_pat": accessions[paternal],
            "sex": row.sex,
        })
        for haplotype, accession in accessions.items():
            url_rows.append({
                "accession": accession,
                "url": find_haplotype_fasta_url(folder, row.sample_id, haplotype),
            })

    pd.DataFrame(assembly_rows).to_csv("hgsvc3_assemblies.tsv", sep="\t", index=False)
    pd.DataFrame(url_rows).to_csv("hgsvc3_assembly_urls.tsv", sep="\t", index=False)
    print(f"Wrote {len(assembly_rows)} samples to hgsvc3_assemblies.tsv and {len(url_rows)} URLs to "
          f"hgsvc3_assembly_urls.tsv")


if __name__ == "__main__":
    main()
