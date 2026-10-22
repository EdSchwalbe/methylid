#!/usr/bin/env Rscript
# Match methylation-array samples to SNP-array or sequencing samples by genotype.
#
# Usage:
#   match_samples.R --beta BETA.tsv --features FEAT.tsv --geno SRC:PATH \
#                   --out PREFIX [options]
#
#   --beta      beta matrix, probes x samples, first column = probe ID. Needs the
#               rs/nv probes and, for the extended marker set, the CpG/SBE-SNP
#               probes; extract_meth_snps.R writes both from IDATs.
#   --features  probe -> rsID table from build_reference.R
#   --geno      plink:PREFIX | vcf:FILE | snp6:DIR (CHP call files from
#               read_calvin_chp.py; needs --snp6-ann)
#   --sites     for vcf: a rsid/chrom/pos table from fetch_positions.py, so
#               markers are found by position when the VCF ID column is "."
#               (the norm for release and pipeline VCFs). Build must match the
#               VCF: pass grch37 or grch38 to fetch_positions.py accordingly.
#   --seeds     TSV with columns meth,snp giving known same-patient pairs. If
#               omitted, seeds are bootstrapped from the rs probes alone.
#   --min-per-snp  call threshold on mean log10 LR per SNP compared (default
#               0.08). This is the portable statistic: a true pair averages
#               about +0.27 per SNP on any platform, an unrelated pair is
#               negative once a few dozen SNPs are compared.
#   --min-gap   minimum total log10 LR lead over the runner-up (default 5),
#               which separates a genuine match from a close relative.
#   --min-snps  warn below this many shared SNPs (default 50): with fewer, the
#               best of many candidates can reach a true pair's per-SNP score.
#   --sex       optional TSV (sample,sex) for methylation samples, and
#   --geno-sex  optional TSV (sample,sex) for genotype samples: reported as an
#               independent check, never used for matching.
suppressMessages(library(data.table))
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "match_lr.R"))
source(file.path(here, "match_meth_snp.R"))
source(file.path(here, "genotypes.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

opt <- local({
  a <- commandArgs(TRUE)
  keys <- grep("^--", a)
  setNames(as.list(a[keys + 1]), sub("^--", "", a[keys]))
})
req <- c("beta", "features", "geno", "out")
if (!all(req %in% names(opt))) stop("missing: ", paste(setdiff(req, names(opt)), collapse = ", "))
min_per_snp <- as.numeric(opt$`min-per-snp` %||% 0.08)
min_gap <- as.numeric(opt$`min-gap` %||% 5)
min_snps <- as.numeric(opt$`min-snps` %||% 50)

feat <- fread(opt$features)
beta <- as.matrix(fread(opt$beta), rownames = 1)
message(sprintf("methylation: %d probes x %d samples", nrow(beta), ncol(beta)))

src <- sub(":.*", "", opt$geno)
path <- sub("^[^:]*:", "", opt$geno)
want <- unique(feat$rsid)
dos <- switch(src,
  plink = read_plink_dosage(path, want),
  vcf   = read_vcf_dosage(path, want,
                          sites = if (!is.null(opt$sites)) fread(opt$sites) else NULL),
  snp6  = read_snp6_dosage(list.files(path, "\\.tsv\\.gz$", full.names = TRUE),
                           fread(opt$`snp6-ann`)),
  stop("--geno must start with plink:, vcf: or snp6:")
)

feat <- feat[feat$feature %in% rownames(beta) & feat$rsid %in% rownames(dos), ]
message(sprintf("usable features: %d probes over %d SNPs (%s)", nrow(feat),
                uniqueN(feat$rsid), paste(names(table(feat$kind)), table(feat$kind),
                                          collapse = ", ")))
if (uniqueN(feat$rsid) < 20) stop("too few shared SNPs (", uniqueN(feat$rsid), ")")
if (uniqueN(feat$rsid) < min_snps)
  warning("only ", uniqueN(feat$rsid), " shared SNPs before calibration; ",
          "matches may not separate cleanly from unrelated samples", call. = FALSE)

## Seeds: supplied truth set, or bootstrapped from the rs probes
if (!is.null(opt$seeds)) {
  seeds <- as.data.frame(fread(opt$seeds))
  seeds <- seeds[seeds$meth %in% colnames(beta) & seeds$snp %in% colnames(dos), ]
  message("seed pairs supplied: ", nrow(seeds))
} else {
  rs_feat <- feat[feat$kind == "rs", ]
  if (nrow(rs_feat) < 8)
    stop("only ", nrow(rs_feat), " shared rs probes: supply --seeds with known pairs")
  b_rs <- beta[rs_feat$feature, , drop = FALSE]
  rownames(b_rs) <- rs_feat$rsid
  sm <- match_meth_snp(b_rs, dos[rs_feat$rsid, , drop = FALSE])$matches
  seeds <- data.frame(meth = sm$meth_sample, snp = sm$best_snp_sample)[
    sm$best_score <= 0.09 & sm$gap >= 0.03, ]
  message(sprintf("bootstrapped %d seed pairs from %d rs probes", nrow(seeds), nrow(rs_feat)))
  if (nrow(seeds) < 30)
    stop("too few confident seeds (", nrow(seeds), "): supply --seeds with known pairs")
}

m <- match_lr(beta, dos, feat, seeds)
res <- as.data.table(m$matches)
res[, margin := best_llr - second_llr]
res[, called := per_snp >= min_per_snp & margin >= min_gap]

## Independent checks, reported only
add_sex <- function(res, file, key, col) {
  if (is.null(file)) return(res)
  s <- fread(file)
  setnames(s, 1:2, c(key, col))
  merge(res, unique(s), by = key, all.x = TRUE)
}
res <- add_sex(res, opt$sex, "meth_sample", "meth_sex")
res <- add_sex(res, opt$`geno-sex`, "best_snp_sample", "geno_sex")
if (all(c("meth_sex", "geno_sex") %in% names(res)))
  res[, sex_agree := toupper(substr(meth_sex, 1, 1)) == toupper(substr(geno_sex, 1, 1))]

# a genotype sample claimed by two methylation samples is either a replicate
# pair or a labelling error: flag rather than silently keep the better one
res[called == TRUE, dup_geno := .N > 1, by = best_snp_sample]
res <- res[order(-best_llr)]

fwrite(res, paste0(opt$out, "_matches.tsv"), sep = "\t")
fwrite(as.data.table(m$features), paste0(opt$out, "_features.tsv"), sep = "\t")
saveRDS(m, paste0(opt$out, "_lr.rds"))

cat(sprintf("\nfeatures calibrated: %d (median R2 %.2f)\n", nrow(m$features),
            median(m$features$r2)))
cat(sprintf("SNPs compared per pair: %d-%d\n", min(res$n_snps), max(res$n_snps)))
null <- m$llr
null[cbind(match(res$meth_sample, rownames(null)),
           match(res$best_snp_sample, colnames(null)))] <- NA
cat(sprintf("best non-match total log10 LR anywhere: %.1f\n", max(null, na.rm = TRUE)))
cat(sprintf("called matches: %d of %d methylation samples\n", sum(res$called), nrow(res)))
if (sum(res$called)) {
  cat(sprintf("  per-SNP log10 LR %.3f to %.3f; total LR %.1f to %.1f\n",
              min(res[called == TRUE, per_snp]), max(res[called == TRUE, per_snp]),
              min(res[called == TRUE, best_llr]), max(res[called == TRUE, best_llr])))
  cat(sprintf("  highest per-SNP score among uncalled samples: %.3f\n",
              suppressWarnings(max(res[called == FALSE, per_snp], na.rm = TRUE))))
  if ("sex_agree" %in% names(res))
    cat(sprintf("  sex agreement: %.1f%% (%d discordant, listed in output)\n",
                100 * mean(res[called == TRUE, sex_agree], na.rm = TRUE),
                res[called == TRUE & sex_agree == FALSE, .N]))
  if (res[called == TRUE & dup_geno == TRUE, .N])
    cat(sprintf("  WARNING: %d genotype samples matched by >1 methylation sample\n",
                uniqueN(res[called == TRUE & dup_geno == TRUE, best_snp_sample])))
}
cat("wrote ", opt$out, "_matches.tsv\n", sep = "")
