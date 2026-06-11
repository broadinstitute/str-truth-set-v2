# Plan — updated HG002 tool-comparison viewer (str-truth-set-v2)

Goal: rerun the HG002 genotyping-tool evaluations against the latest HG002 truth
set and publish a **new** tool-comparison viewer (GitHub Pages) in this repo,
showing only the **Accuracy-by-allele-size** plot type, stratified by **purity**.

### Latest HG002 truth set (use THIS)
`gs://str-truth-set-v2/filter_vcf_v2__2025_12_29/HG002/HG002.tandem_repeat_genotypes.tsv.gz`
(+ `.json.gz`), dated **2026-06-09**, from Hail Batch **8417721**
("genotype (cpu=4): HG002") — latest `str-analysis filter_vcf_to_tandem_repeats
genotype`. **2,214,976 loci** (HG002 genotyped vs the combined 321-sample
catalog; the variant/"positive" subset feeds tool comparison).

New-format columns (NOT the old comparison-pipeline schema): `Chrom, Start0Based,
End, Locus, LocusId, Motif, CanonicalMotif, MotifSize, NumRepeatsInReference,
NumRepeats{Short,Long}Allele, RepeatSize{Short,Long}AlleleBp, Zygosity,
IsPureRepeat, RepeatPurity, Allele{1,2}Sequence, Allele{1,2}MotifSequence, …`.
Note **`RepeatPurity` is per-locus**, whereas the plot script reads per-allele
`FractionPureRepeats: Allele: Truth` — a mapping/rename is required.

Status: **PLAN ONLY — nothing launched.** A parallel Claude owns the pipeline
edits (`run_tools/run_genotyping_tools.py`, `str-truth-set` `inquistr-integration`
branch). Coordinate before launching any Hail Batch run.

Budget cap for the whole run: **$20**.

---

## 1. Tools × data types (HG002 rows 44–56)

7 tools. Each runs only on compatible data types; the driver enforces this.

| Tool | Mode | Data types (coverages) | Runs |
|---|---|---|---|
| EHv5 (ExpansionHunter) | short read, optimized-streaming | illumina 31/20/10, element 30, ultima 36, exome 3 | 6 |
| **EHv5-bw2-optimized** | short read, optimized-streaming + `--improved-genotyping` | illumina 31/20/10, element 30, ultima 36, exome 3 | 6 |
| GangSTR | short read | illumina 31/20/10, element 30, exome 3 (skips ultima) | 5 |
| HipSTR | short read | illumina 31/20/10, element 30, ultima 36, exome 3 | 6 |
| TRGT | long read (pacbio only), cpu=16 | pacbio 37/30/20/8 | 4 |
| LongTR | long read, cpu=1 | pacbio 37/30/20/8, ONT 26/20/8 | 7 |
| inquiSTR | long read, cpu=16 | pacbio 37/30/20/8, ONT 26/20/8 | 7 |

**Total ≈ 41 genotyping runs.** Each run = genotype step(s) + combine +
add-columns + plot + image-headers. HG002 truth set = **139,832 loci**
(excluding homopolymers).

---

## 2. Gaps to close BEFORE running (pipeline owner)

1. **Bridge the new truth set → comparison format.** The add-columns step expects
   `…/HG002/HG002.STRs.excluding_homopolymers.annotated.variants.for_comparison.tsv.gz`
   (+ `.alleles`) with per-allele `FractionPureRepeats` and `DiffFromRefRepeats`.
   The latest truth set is the **new `tandem_repeat_genotypes.tsv.gz` format**
   (per-locus `RepeatPurity`, short/long-allele fields). Need a conversion that
   produces the `annotated.variants/alleles` + `for_comparison` tables the
   comparison scripts consume, mapping `RepeatPurity` (per-locus) →
   `FractionPureRepeats: Allele: Truth` (per-allele). Without the purity column
   the plot script collapses purity to a single "all" bin (no stratification).
   The stale 2024-07-29 `for_comparison` table on GCS must NOT be reused.
2. **EHv5 tool token.** Driver passes `--tool EHv5` to
   `plot_tool_accuracy_by_allele_size.py`, whose `--tool` choices are
   `{ExpansionHunter, GangSTR, HipSTR, TRGT, LongTR, inquiSTR, …}` — **`EHv5` is
   not a choice** → argparse error. Map `EHv5 → ExpansionHunter` (and pick a token
   for `EHv5-bw2-optimized`) in the plot + add-columns steps. The viewer config
   (`fileToken`) must match whatever string is finally emitted.
3. **Genotype dimension.** The plot step hardcodes `--genotype all`, so only
   `.all_genotypes` SVGs are produced. The viewer's Genotype facet
   (het/hom/multi) needs the plot step to **drop `--genotype all`** (the script
   then loops all four). Otherwise remove the Genotype facet from the viewer.
4. **Add EHv5-bw2-optimized.** Register in `SHORT_READ_TOOLS`; add a code path
   running ExpansionHunter `--analysis-mode optimized-streaming --improved-genotyping`;
   distinct output dir + plot/add-columns token. Confirm the
   `FILTER_VCFS_DOCKER_IMAGE` / EH docker has an ExpansionHunter build that
   supports `--improved-genotyping`.
5. **Regenerate per-tool catalogs from the new truth set.** The driver reads
   positive-loci catalogs from `--filter-vcf-dir gs://str-truth-set-v2/filter_vcf`
   (`…/HG002/HG002.STRs.excluding_homopolymers.positive_loci.{EHv5,GangSTR,HipSTR,
   TRGT,LongTR,inquiSTR}…`), currently dated **2026-01-26** and derived from the
   OLD truth set. Re-export per-tool catalogs from the **2026-06-09** truth set's
   variant/positive-loci subset, then either overwrite that path or point
   `--filter-vcf-dir` at the new location. The positive-loci count (subset of the
   2.2M genotyped loci) sets tool-genotyping runtime → reconfirm cost once known.

---

## 3. Cost estimate

Hail Batch spot ≈ **$0.04/cpu-hr** (measured from prior HG002 batches: HipSTR
31x = $0.065, EHv5 ≈ $0.06/run). Long-read cpu=16 tools dominate.

| Group | Runs | ~Cost |
|---|---|---|
| Short read (EHv5, EHv5-opt, GangSTR, HipSTR) | 23 | ~$1.5 |
| TRGT (cpu=16, pacbio) | 4 | ~$0.8 |
| LongTR (cpu=1) | 7 | ~$1.4 |
| inquiSTR (cpu=16) | 7 | ~$1.8 |
| Downstream (combine/add-cols/plot/headers) | 41 | ~$1.2 |
| Truth-set table regen (one-time) | — | ~$0.1–0.3 |
| **Total (typical)** | | **~$6–9** |

With retries/failures buffer, still **well under the $20 cap**. Mitigation:
monitor batch cost live; cancel if trending > ~$15.

---

## 4. Execution sequence

1. Pipeline owner closes gaps §2.1–2.4 on `inquistr-integration`; push.
2. **Dry run**: driver with `--only-print-total-number-of-plots`, or one
   `(tool,row)`, to confirm SVG filenames + that purity/genotype variants land on
   GCS under the expected paths.
3. **Full run** (one Hail Batch; samples run in parallel):
   ```
   cd run_tools
   python3 run_genotyping_tools.py \
     --sample-id HG002 --exclude-homopolymers \
     --data-type illumina --data-type illumina_exome --data-type element \
     --data-type ultima --data-type pacbio --data-type ONT \
     --tool EHv5 --tool EHv5-bw2-optimized --tool GangSTR --tool HipSTR \
     --tool TRGT --tool LongTR --tool inquiSTR
   ```
   Watch cost in the Hail Batch UI; cancel if > ~$15 trending.
4. Verify a sample of expected SVG URLs return HTTP 200 from
   `https://storage.googleapis.com/str-truth-set-v2/tool_results/…`.

---

## 5. Viewer build + deploy

- Scaffold: `docs/tool_comparison_viewer.html` (this repo). Facets: **Tool (7),
  Sequencing data, Coverage, Motif size, Genotype, Purity (7 bins + all)** +
  "Hide No-Call loci". **Dropped** vs v1: Q-threshold, Pure/Interrupted.
- URL schema (single `buildImageUrl()` function):
  ```
  BASE/{data_type}/{tool_dir}/{cov}x_coverage/
    tool_accuracy_by_true_allele_size.{motif}_motifs.{genotype}.{cov}x[.purity_{p}][.exclude_no_call_loci].{fileToken}.svg
  ```
- **Lock the template to one verified real SVG** before announcing (the EHv5
  token, motif-token form, and genotype availability are the moving parts).
- GitHub Pages: Settings → Pages → Deploy from branch → `main` `/docs`.
  Page URL: `https://<owner>.github.io/str-truth-set-v2/tool_comparison_viewer.html`.
- Bucket `gs://str-truth-set-v2` is already public-read (verified HTTP 200).
  Ensure new SVG objects are public too (uniform bucket-level access inherits).

---

## 6. Verification checklist

- [ ] for_comparison table regenerated, has `FractionPureRepeats`.
- [ ] Dry run shows purity-suffixed + genotype-suffixed SVGs on GCS.
- [ ] Full run finished under budget; cost recorded.
- [ ] Spot-check SVG URLs return 200.
- [ ] `buildImageUrl()` tokens match real filenames (esp. EHv5 / EHv5-opt).
- [ ] GitHub Pages enabled; page loads; facet switching shows correct plots.
- [ ] No-plot combinations degrade gracefully (the `#plot-status` message).

## 7. Open coordination items

- Who launches the batch (this instance vs the parallel one) — pending your call.
- Who owns the truth-set→catalog/for_comparison conversion (gaps §2.1, §2.5) from
  the 2026-06-09 `tandem_repeat_genotypes.tsv.gz` — likely a str-analysis export
  step; confirm it exists or needs writing.
- Confirm Genotype facet wanted (drives gap §2.3) — keep all four, or "all" only.
