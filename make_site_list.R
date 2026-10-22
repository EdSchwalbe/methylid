# Write the site lists to send to whoever is supplying sequencing data, so they
# can return a small targeted file instead of a whole WGS callset.
#
# Usage: Rscript make_site_list.R PLATFORM OUTDIR
#   PLATFORM  450k | epic | epicv2
#
# Produces in OUTDIR:
#   sites_<platform>_rsids.txt    every candidate marker rsID, one per line.
#                                 Build-agnostic, so this is the primary ask.
#   sites_<platform>_hg19.bed     positions for the subset whose coordinates we
#                                 can resolve locally (those also on SNP6).
#                                 A fallback for callsets with no rsIDs.
#   REQUEST.md                    what to ask for, in a form that can be
#                                 forwarded as-is.
suppressMessages({
  library(data.table)
  library(DBI)
  library(pd.genomewidesnp.6)
})

args <- commandArgs(trailingOnly = TRUE)
platform <- tolower(args[1])
outdir <- if (length(args) > 1) args[2] else "request"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

feat <- fread(file.path("ref", sprintf("feat_%s.tsv", platform)))
rsids <- sort(unique(feat$rsid))
writeLines(rsids, file.path(outdir, sprintf("sites_%s_rsids.txt", platform)))

# hg19 coordinates, for the markers that are also on SNP6 (pd.genomewidesnp.6 is
# HG19). The rest would need a dbSNP lookup, which is why rsIDs are the main ask.
con <- db(pd.genomewidesnp.6)
fs <- as.data.table(dbGetQuery(con,
  "select dbsnp_rs_id, chrom, physical_pos from featureSet where dbsnp_rs_id is not null"))
fs <- fs[!is.na(physical_pos) & dbsnp_rs_id %in% rsids]
fs <- fs[!duplicated(dbsnp_rs_id)][order(chrom, physical_pos)]
bed <- fs[, .(chrom = paste0("chr", chrom), start = physical_pos - 1L,
              end = physical_pos, name = dbsnp_rs_id)]
fwrite(bed, file.path(outdir, sprintf("sites_%s_hg19.bed", platform)),
       sep = "\t", col.names = FALSE)

message(sprintf("%s: %d candidate SNPs; %d with local hg19 coordinates",
                platform, length(rsids), nrow(bed)))

request <- sprintf('# Data request: targeted genotypes for sample matching

We are matching methylation-array samples to sequencing samples by genotype, to
work out which came from the same patient. We need genotypes at a fixed list of
common SNPs -- **not** a whole-genome callset.

## What we need

Genotypes at the %d sites in `sites_%s_rsids.txt` (rsIDs, one per line). An
hg19 BED of the subset we could resolve locally is included as
`sites_%s_hg19.bed` if positions are easier for you than rsIDs.

The resulting file should be a few MB at most.

## Please make sure of these four things

1. **Emit a call at every requested site, including homozygous-reference.**
   This is the one that matters most. A variants-only VCF omits sites where the
   sample matches the reference, and we measured that this loses about a third
   of the usable genotypes, because homozygous-reference is just as informative
   to us as a variant. With `bcftools mpileup -R sites.bed` or GATK
   `--output-mode EMIT_ALL_SITES`, every requested site appears.

2. **Distinguish "homozygous reference" from "no coverage".** These must not
   both come back as absent or as `./.`. We treat a missing genotype as unknown
   and drop it, so an uncovered site that looks like `0/0` would actively
   mislead us, while a genuine `0/0` that is absent merely costs us power.

3. **Put rsIDs in the ID column**, e.g.
   `bcftools annotate -c ID -a dbsnp.vcf.gz`. We join on rsID, which keeps us
   independent of genome build. If the ID column is `.` we cannot use the file
   without a position-based fallback.

4. **One multi-sample VCF is ideal**, with the sample names we can map to
   patient or specimen identifiers. Per-sample files are fine too, as long as
   point 1 holds for each.

## Formats

- **Good:** VCF (plain or bgzipped) with all requested sites called.
- **Also good:** PLINK `.bed/.bim/.fam` — we read it directly and it is compact.
- **Please avoid:** gVCF (reference blocks need converting first), and
  variants-only VCFs (see point 1).
- **Not needed:** BAM/CRAM, full WGS VCFs, or anything beyond these sites.

## A note on scope

These are common germline SNPs chosen only to fingerprint sample identity. The
list carries no disease-associated or clinically actionable content, and the
targeted file is far less disclosive than a whole-genome callset, which usually
makes the governance side easier. We are happy to work from whatever your
information-governance arrangements allow.
', length(rsids), platform, platform)

writeLines(request, file.path(outdir, "REQUEST.md"))
message("wrote ", file.path(outdir, "REQUEST.md"))
