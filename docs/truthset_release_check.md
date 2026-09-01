# HG002 truth set release check

Run this before publishing a truth set release. It scores the HG002 truth set against an independent benchmark and
compares the result to a checked-in baseline, so a change to the catalog or the filter steps cannot quietly degrade
the truth set that every tool ranking downstream is measured against.

It takes about 30 seconds.

## What it measures

[Aardvark](https://github.com/PacificBiosciences/aardvark) compares two VCFs by edit distance over haplotype
sequences, so two callers that describe the same haplotype with different records score the same. The check uses
GIAB's T2T-Q100 v1.1 benchmark as the truth and this project's own truth set as the query, which means a false
negative is a variant the truth set is missing and a false positive is one it invented. Both are truth-set defects.

The scored regions are the catalog loci padded 50bp, intersected with both the GIAB benchmark BED and dipcall's
confident regions BED.

55 metrics are guarded, all of them rates in `[0, 1]` where higher is better:

| Metric | Why it is guarded |
| --- | --- |
| `BASEPAIR/ALL/{recall,precision,f1}` | The headline quality number, representation-invariant. |
| `GT/ALL/{recall,precision,f1}` | The hap.py-like per-genotype score, which is what most other benchmarking tools report. |
| `BASEPAIR/<chrom>/recall`, `GT/<chrom>/recall` | Per-chromosome recall for chr1 through chr22, chrX and chrY. |
| `GT/ALL/zygosity_agreement` | The share of truth variants that are not pure zygosity errors, `1 - truth_fn_gt / truth_total`. |

The per-chromosome rows are the point. A genome-wide number is too blunt for a whole chromosome breaking: when
dipcall's `.|1` encoding gave chrX and chrY a GT recall of exactly 0, genome-wide BASEPAIR F1 only moved from 0.9932
to 0.9830, which a loose threshold would miss. Per chromosome the same failure shows up as a drop from 0.987 to
0.000.

## Prerequisites

`bcftools`, `tabix`, `bedtools` and `curl` on the `PATH`, plus an `aardvark` binary.

Aardvark v1.0.0 ships only a linux-x86_64 binary, so on macOS it has to be built from source. Its `build.rs` pulls in
`vergen-gitcl` 9.1.0, which needs rustc 1.88, and Homebrew currently ships 1.87. Since `vergen` only supplies a
version string, stubbing it out is enough:

```bash
git clone --depth 1 --branch v1.0.0 https://github.com/PacificBiosciences/aardvark.git
cd aardvark
cp build.rs build.rs.orig && cp Cargo.toml Cargo.toml.orig
cat > build.rs <<'EOF'
fn main() {
    println!("cargo:rustc-env=VERGEN_GIT_DESCRIBE=v1.0.0-localbuild");
    println!("cargo:rerun-if-changed=Cargo.toml");
    println!("cargo:rerun-if-changed=src");
}
EOF
# drop the vergen-gitcl line from [build-dependencies] in Cargo.toml, then:
cargo build --release
```

No comparison logic is touched by this, only the version string baked into `--version`.

## Running it

```bash
cd filter_vcfs_v2

python3 run_truthset_release_check.py \
    --truth-vcf results/HG002.high_confidence_regions.vcf.gz \
    --catalog-bed results/combined.321_catalogs.merged.tandem_repeats.bed.gz \
    -R /path/to/hg38.fa \
    --dipcall-confident-bed gs://str-truth-set-v2/dipcall_pipeline/HG002/HG002.dip.bed.gz \
    --aardvark /path/to/aardvark
```

It exits 0 on pass and 1 on regression. The GIAB benchmark files are downloaded into `--work-dir` on the first run
and reused after that.

Pass `--verbose` to list the per-chromosome metrics that passed, which are summarized as a count by default.

## When it fails

The output names every metric that dropped by more than `--tolerance` (0.002 by default). Two failure shapes are
worth recognizing:

**chrX and chrY regress together, and everything else passes.** The truth set was almost certainly built without the
haploid genotype normalization step. dipcall fills in both haplotype slots on every record, so a male sample's chrX
outside the PAR and all of its chrY come out as `.|1`, which no benchmarking tool can match. The check prints this
diagnostic itself when the VCF carries no `##normalizeHaploidGenotypes` header. Confirm that
`run_filter_vcf_to_tandem_repeats.py` ran `normalize_haploid_genotypes.py` and that the sample's sex is set in the
metadata table.

**One or two autosomes regress.** Look at the catalog first, since the scored regions are derived from it. A catalog
change that moves or drops loci changes which variants are in scope.

## Updating the baseline

Only when the truth set legitimately changed, and never to make a failure go away:

```bash
python3 run_truthset_release_check.py ... --update-baseline
```

Commit `truthset_release_check_baseline.tsv` with a message saying what changed and why the new numbers are correct.
These numbers are what every later release is judged against, so an unexplained baseline bump removes the check's
value entirely.

## Known caveats

- **HG002 only.** It is the only sample with an independent benchmark of this quality. The other 320 assemblies are
  untested by this method.
- **Restricted to the GIAB benchmark BED**, which deliberately excludes the hardest regions. The true error rate over
  the full catalog is higher than these numbers suggest, by an unknown amount.
- **Aardvark drops variants not fully contained in the scored regions.** That is unbiased between truth and query,
  which is why the catalog intervals are padded 50bp before intersecting, but it does reduce coverage.
- **GIAB Q100 is itself a draft benchmark**, not ground truth. Some disagreements will be GIAB's.
- **The check uppercases REF and ALT before scoring.** dipcall carries hg38's soft-masking through into the alleles,
  and the catalog lives in repeat regions, which is exactly where hg38 is soft-masked. The filter step now uppercases
  them (`uppercase_ref_and_alt.py`), so a VCF built by the current pipeline has none left. The check still uppercases,
  both so it can score VCFs produced before that change and so its numbers stay comparable across the change, and it
  reports how many records it had to rewrite. A non-zero count on a fresh release means the filter step did not run
  that stage.
