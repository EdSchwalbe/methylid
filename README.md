# Matching methylation arrays to SNP arrays / sequencing by genotype

Identifies which methylation-array samples and which SNP-array or sequencing
samples come from the same patient, using genotypes read from both. No
clinicopathological matching required.

Validated on 763 Cavalli 450K samples (GSE85212) against 1,097 Northcott SNP6
samples (GSE37384): **488 matches called, 0 ambiguous, sex agreement 98.4%**,
with a gap of 92 orders of magnitude between the weakest true match and the
strongest non-match.

## Prior art — read this first

**This is not the first tool to do cross-omics genotype fingerprinting, and if
you are starting fresh you should try the established packages first.**

- **[omicsPrint](https://bioconductor.org/packages/omicsPrint)** (van Iterson
  et al., *Bioinformatics* 2018, [doi:10.1093/bioinformatics/bty062](https://doi.org/10.1093/bioinformatics/bty062))
  — identity-by-state matching across omics types, built for exactly this
  purpose. The closest equivalent to this repository.
- **[meffil](https://github.com/perishky/meffil)** — `meffil.snp.concordance()`
  compares rs-probe betas against PLINK genotypes.
- **[ewastools](https://github.com/hhhh5/ewastools)** (Heiss & Just, *Clin
  Epigenetics* 2018) — mixture-model genotype calling robust to noise and
  contamination, plus duplicate and mislabel detection.
- **[sesame](https://bioconductor.org/packages/sesame)** — `formatVCF()` writes
  methylation-derived genotypes to VCF for comparison with standard tools.
- **Zhou, Laird & Shen, *NAR* 2017** ([doi:10.1093/nar/gkw967](https://doi.org/10.1093/nar/gkw967))
  — the reference for using SNP-affected Infinium probes as genotype markers,
  which is the basis of the extended marker set here.

What this repository adds is narrower: Affymetrix SNP6 Birdseed CHP input,
scoring that tolerates LOH and allelic imbalance in tumours, a platform-portable
per-SNP threshold with a measured minimum marker count, and position-based VCF
joining for files with no rsIDs. It was written for medulloblastoma tumour
cohorts where only 19 of the 65 450K rs probes are genotyped on SNP6. It has
**not** been benchmarked against the packages above — if you do that, the
result would be worth knowing.

## Validation status

| Combination | Status |
|---|---|
| 450K × SNP6 | **Validated end to end** on 763 × 1,097 real samples: 488 matches, 0 ambiguous |
| Few-marker case (303 SNPs, 3 rs probes, seeds supplied) | **Validated in simulation** on real genotypes and betas: 486/488 recovered, 0 false positives |
| EPICv2 extraction | **Works** on real EPICv2 IDATs: 65 rs + 7,383 probes, betas trimodal |
| 450K / EPIC / EPICv2 marker tables | **Verified**: every feature probe present on its array |
| PLINK, VCF, SNP6 loaders | **Round-trip tested** against real data |
| Folder of per-sample VCFs | **Validated** on 1,097 per-sample files: 488/488 the same partner, 0 disagreements |
| EPIC v1 IDAT reading | *Untested* — low risk (routine in minfi; all 11,661 feature probes confirmed on the manifest) |
| **Any GSA data** | *Untested* — no manifest available. Run step 1 below against your `.bim` to get the real overlap |
| EPIC/EPICv2 × genotypes end to end | *Untested* — needs one pilot with paired samples |

Two things decide the GSA case, both answered by one run of `build_reference.R`
with your `.bim`: how many candidate SNPs are directly typed (needs ≥50), and
how many of the 59–65 rs probes are present (needs ≥8 to bootstrap seeds, else
pass `--seeds` with known pairs — which works fine, as the simulation above
shows).

## Why more than the 65 rs probes

The 59–65 dedicated `rs` probes are the obvious markers, but only **19 of them
are directly genotyped on Affymetrix SNP6**, which is not enough on its own
(see the marker-count table below). The usable set is much larger: any probe
with a common SNP at the interrogated C/G or at the single-base-extension
position has a beta that tracks that genotype. Candidate markers per platform:

| Methylation platform | rs probes | candidate SNPs |
|---|---|---|
| 450K | 65 | 7,675 |
| EPIC v1 | 59 | 11,713 |
| EPIC v2 | 62 | 7,416 |

Against SNP6 this yielded 485 shared SNPs, of which 275 survived calibration
(median R² 0.81) — about 250–280 compared per pair, against 12–19 for rs-only.

Note on EPIC v2: it keeps 83% of EPIC v1 probes in general but only **0.4%** of
those with a common SNP at the CpG, because Illumina removed SNP-affected
probes. Its markers therefore come from its own dbSNP 151 / hg38 annotation at
different probes, and the 450K's 65 rs probes are restored.

## How it scores a pair

Each marker's beta ~ dosage relationship is fitted on training pairs, markers
with R² below 0.5 are dropped, and the best marker per SNP is kept. A pair is
then scored

    log10 LR = sum over markers of log10 [ p(beta | genotype) / p(beta) ]

with p a normal mixed with a 2% uniform outlier component. Three consequences
matter for tumour work:

- **Allele orientation is learned, not looked up.** The sign of the fitted
  slope settles which allele a probe's high beta reports, so strand and A/B
  allele conventions never have to be reconciled. This is where hand-rolled
  implementations usually go wrong.
- **LOH and allelic imbalance cannot veto a true match.** The outlier component
  bounds how much any single discordant marker can cost, so losing
  heterozygosity at some loci degrades the score without destroying it. This
  matters for SHH (9q, 10q, 17p) and Group 3/4 (i17q) medulloblastoma.
- **Purity and FFPE degradation are absorbed** by the fitted residual SD rather
  than by fixed beta cut-offs.

Calibration is cross-fitted: seed pairs are split in two, and each seed sample
is scored with markers fitted on the other half, so a seed cannot be scored
with its own fit.

## How many SNPs you need

Sensitivity and false calls against the 488-pair truth set, markers randomly
downsampled, 5 replicates (call rule: per-SNP ≥ 0.08, margin ≥ 5):

| shared SNPs | sensitivity | false calls | weakest true (per SNP) | strongest non-match |
|---|---|---|---|---|
| 10 | 0.00 | 0.0 | 0.140 | 0.436 |
| 20 | 0.37 | 0.6 | 0.044 | 0.237 |
| 30 | 0.94 | 1.0 | 0.032 | 0.142 |
| **50** | **1.00** | **0.0** | 0.054 | 0.030 |
| 100 | 1.00 | 0.0 | 0.067 | −0.047 |
| 275 | 1.00 | 0.0 | 0.111 | −0.038 |

**Use at least 50 shared SNPs.** Below ~30, the best of many candidates can
reach a true pair's score and false matches appear.

The call rule thresholds **log10 LR per SNP**, not total LR: a true pair
averages about +0.27 per SNP regardless of platform or marker count, while an
unrelated pair goes negative. That makes one threshold (0.08) portable across
platform combinations. Total LR is reported too, and is what gives the large
margins, but it scales with marker count so it is not portable.

## Usage

Output paths below are written as `work/...`. That is only a convention for
scratch output — nothing in the repository depends on it, the scripts create
whatever directory you name, and `work/` is gitignored so it is absent from a
fresh clone. Point `--out` anywhere you like.

### A folder of IDATs against a folder of VCFs

Three commands. Set `PLATFORM` to `450k`, `epic` or `epicv2`, and `BUILD` to
the genome build of the VCFs (`grch38` for most modern callsets).

```bash
PLATFORM=epic; BUILD=grch38

# 1. marker positions, so markers are found even when the VCF ID column is "."
python3 fetch_positions.py ref/feat_$PLATFORM.tsv ref/pos_$PLATFORM.tsv $BUILD

# 2. betas and inferred sex from the IDAT folder (last arg = marker probe list)
cut -f1 ref/feat_$PLATFORM.tsv | tail -n +2 | sort -u > /tmp/probes.txt
Rscript extract_meth_snps.R /path/to/idat_folder work/mine 16 /tmp/probes.txt

# 3. match: --geno vcf: takes a directory of VCFs as well as a single file
Rscript match_samples.R \
  --beta     work/mine_beta_all.tsv \
  --features ref/feat_$PLATFORM.tsv \
  --sites    ref/pos_$PLATFORM.tsv \
  --geno     vcf:/path/to/vcf_folder \
  --sex      work/mine_sex.tsv \
  --out      work/mine
```

Results land in `work/mine_matches.tsv`. Validated on 1,097 per-sample VCFs:
identical to the PLINK route, 488/488 the same partner, 0 disagreements.

**If the VCFs are variants-only this will fail**, and deliberately so: without
homozygous-reference calls there is too little left per marker to calibrate.
The reader warns first ("<10% homozygous-reference calls"), then the run stops
with an explanation rather than returning weak matches. See
"Asking a colleague for sequencing data" below.

### Against a SNP array

```bash
# markers your genotyping platform actually carries
cut -f2 your_gsa.bim > gsa_rsids.txt
Rscript build_reference.R epic ref/feat_epic_gsa.tsv 0.05 gsa_rsids.txt

Rscript extract_meth_snps.R /path/to/idats work/mine 16 /tmp/probes.txt

Rscript match_samples.R \
  --beta work/mine_beta_all.tsv \
  --features ref/feat_epic_gsa.tsv \
  --geno plink:/path/to/your_gsa \
  --out work/mine \
  --sex work/mine_sex.tsv --geno-sex work/gsa_sexmap.tsv
```

`--geno` accepts `plink:PREFIX` (reads only the needed variants by seeking, so
a full GSA or imputed set is fine), `vcf:FILE` for sequencing, or
`snp6:DIR` for Affymetrix CHP calls.

For a VCF, also pass `--sites`. Most real VCFs have `.` in the ID column —
every one of the 1,103,547 chr22 variants in the 1000 Genomes phase 3 release
does — so markers have to be found by position:

```bash
python3 fetch_positions.py ref/feat_epic.tsv ref/pos_epic.tsv grch38  # or grch37
Rscript match_samples.R ... --geno vcf:cohort.vcf.gz --sites ref/pos_epic.tsv
```

Match the build to the VCF. Running the full pipeline through a VCF gives
byte-identical results to the PLINK route (488/488 matches, 0 disagreements),
including with phased genotypes, `chr`-prefixed contigs and multi-allelic
sites. Supply `--seeds` with known same-patient
pairs if you have them; otherwise seeds are bootstrapped from the rs probes,
which needs at least 8 of them shared with your genotyping platform.

Output `*_matches.tsv` gives each methylation sample's best genotype sample,
`per_snp` and total `best_llr`, the margin over the runner-up, `called`, a
`dup_geno` flag where one genotype sample is claimed by two methylation
samples, and the sex check if supplied.

### Asking a colleague for sequencing data

Run `Rscript make_site_list.R epic request_epic` and forward the generated
`REQUEST.md`. It asks for genotypes at our site list rather than a whole WGS
callset, and spells out the one requirement people get wrong: **emit
homozygous-reference calls**. A variants-only VCF omits sites matching the
reference, which loses ~33% of usable genotypes here, since hom-reference is as
informative as a variant. rsIDs must be in the ID column, because
`read_vcf_dosage()` joins on it. PLINK output from the sequencing data is
equally good and avoids the ambiguity entirely.

### For GSA specifically

Only a few hundred of the ~7,400–11,700 candidate SNPs are likely to be
directly typed on a 650K GSA, which is comfortably above the 50 needed. If the
direct overlap turns out to be thin, imputing the GSA data (TOPMed or Michigan)
recovers essentially all of them. Check the real number with step 1 above —
that prints the overlap before you run anything else.

## Files

| File | Purpose |
|---|---|
| `match_samples.R` | CLI driver |
| `build_reference.R` | probe → rsID marker table per platform |
| `extract_meth_snps.R` | betas + inferred sex from IDATs |
| `match_lr.R` | calibration, likelihood-ratio scoring, cross-fitting |
| `match_meth_snp.R` | rs-only concordance matcher, used to bootstrap seeds |
| `genotypes.R` | PLINK / VCF / SNP6 dosage loaders |
| `read_calvin_chp.py` | Affymetrix Calvin CHP reader (see below) |
| `make_site_list.R` | builds the site list + `REQUEST.md` to send to a sequencing data provider |
| `fetch_positions.py` | rsID to chrom:pos table from Ensembl, for VCFs with no rsIDs |
| `ref/feat_*.tsv` | prebuilt marker tables for 450K, EPIC and EPICv2 |
| `ref/pos_450k_grch37.tsv` | prebuilt GRCh37 positions for the 450K markers |
| `exploration/` | the original one-off GEO analysis scripts and a PLINK writer used to round-trip the reader; kept for provenance, superseded by `match_samples.R` |

`affxparser::readCcg` fails on the GEO Birdseed CHPs (it mis-sizes their string
columns), hence `read_calvin_chp.py`. It also extracts the chip-summary header
fields, including APT's own `computed_gender` — a free, independent sex call on
the SNP-array side.

## Licence and citation

GPL-3 (`LICENSE`), matching the Bioconductor packages this depends on.

If it is useful in published work, please cite the tools whose annotation and
methods it rests on — `minfi`, `sesame`, and Zhou, Laird & Shen *NAR* 2017 for
the SNP-affected probe content — and consider citing `omicsPrint` if you used
it for comparison. There is no paper for this repository; cite the URL.

## Caveats

- Sex is reported as a check only and never used for matching.
- A tumour/normal pair from one patient matches, so a match means same patient,
  not same sample. Replicates and relapse pairs show up as `dup_geno`.
- First-degree relatives score well above unrelated samples but well below a
  self-match; with ~250 SNPs the margin separates them, with ~50 treat a
  small margin cautiously.
- Seed bootstrapping assumes a reasonable number of true pairs exist in the two
  sets. For a handful of samples, supply `--seeds`.
- GSE37384/GSE85212 sex annotation in GEO disagrees with both arrays in ~27
  cases; the arrays agree with each other. Trust the arrays.
- **A sex mismatch is not evidence against a match.** All 8 sex-discordant
  pairs here are explained by sex-chromosome loss in the tumour, not by
  mislabelling. SNP6's APT sex call uses a Y/X intensity *ratio*, which is
  degenerate: losing an X raises it and losing Y lowers it, so both push a
  sample across the threshold. Every one of the 8 falls in the ambiguous band
  (ratio 0.40–0.60); none of the 464 samples outside that band is discordant.
  The platforms agree on Y content (Spearman 0.785), so the measurements are
  consistent and only the binary calls disagree. Sex-chromosome loss is common
  in these tumours — ~39% of female and ~14% of male samples here — so treat
  the sex check as a weak prior that is uninformative near the boundary, and
  let the genotype score decide. `work/sex_resolved.tsv` has the per-sample
  X/Y intensities and chrX heterozygosity.
