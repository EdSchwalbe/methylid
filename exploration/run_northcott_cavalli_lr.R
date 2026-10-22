# Extended-marker matching: Cavalli 2017 450K (GSE85212) to Northcott 2012 SNP6
# (GSE37384), using the rs probes plus CpG probes with a common SNP at the CpG or
# single-base-extension position. Seeds for feature calibration come from the
# rs-only matcher, so no external truth set is needed to get started.
# Run from the repository root: Rscript exploration/run_northcott_cavalli_lr.R
# Needs the intermediate files in work/ (see HANDOVER.md); superseded by
# match_samples.R, kept to document how the original GEO analysis was done.
root <- dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
suppressMessages(library(data.table))
source(file.path(root, "match_lr.R"))
source(file.path(root, "match_meth_snp.R"))

geo_meta <- function(path) {
  l <- strsplit(readLines(gzfile(path)), "\t")
  l <- l[sapply(l, `[`, 1) %in% c("!Sample_geo_accession", "!Sample_title", "!Sample_characteristics_ch1")]
  out <- data.frame(gsm = gsub('"', "", l[[which(sapply(l, `[`, 1) == "!Sample_geo_accession")]][-1]),
                    title = gsub('"', "", l[[which(sapply(l, `[`, 1) == "!Sample_title")]][-1]))
  # characteristics rows are not aligned by key across samples: parse each cell
  for (x in l[sapply(l, `[`, 1) == "!Sample_characteristics_ch1"]) {
    v <- gsub('"', "", x[-1])
    for (j in which(grepl(":", v))) {
      key <- make.names(sub(":.*", "", v[j]))
      if (is.null(out[[key]])) out[[key]] <- NA_character_
      out[[key]][j] <- trimws(sub("^[^:]*:", "", v[j]))
    }
  }
  out
}
geo <- "ref/geo"
meta_s <- geo_meta(file.path(geo, "GSE37384_sm.txt.gz"))
meta_m <- geo_meta(file.path(geo, "GSE85212_sm.txt.gz"))

## SNP-array side: dosage matrix over all candidate probesets, autosomes only
ann_rs <- fread("ref/snp6_rs_overlap.tsv")[, .(probeset = man_fsetid, rsid = dbsnp_rs_id, chrom)]
ann_cs <- unique(fread("ref/450k_cpgsnp_snp6.tsv")[, .(probeset = man_fsetid, rsid, chrom)])
ann <- unique(rbind(ann_rs, ann_cs))[chrom %in% as.character(1:22)]

calls <- rbindlist(lapply(list.files("work", "^chpchunk_[a-z]+\\.tsv\\.gz$", full.names = TRUE), fread))
calls[, gsm := sub("_.*", "", sample)]
calls[, dosage := c(`6` = 0, `8` = 1, `7` = 2)[as.character(call)]]
calls <- merge(calls, ann, by = "probeset")
# a few SNPs are on more than one probeset: keep the highest-confidence call
calls <- calls[order(confidence)][!duplicated(calls[order(confidence), .(gsm, rsid)])]
dos <- as.matrix(dcast(calls, rsid ~ gsm, value.var = "dosage"), rownames = "rsid")
message(sprintf("SNP array: %d SNPs x %d samples, call rate %.1f%%",
                nrow(dos), ncol(dos), 100 * mean(!is.na(dos))))

## Methylation side: rs-probe betas plus the CpG/SBE-SNP probe betas
beta_rs <- as.matrix(fread("work/cavalli_rs_beta.tsv"), rownames = 1)
beta_ex <- as.matrix(fread("work/cavalli_extra_beta.tsv"), rownames = 1)
colnames(beta_rs) <- sub("_.*", "", colnames(beta_rs))
colnames(beta_ex) <- sub("_.*", "", colnames(beta_ex))
beta <- rbind(beta_rs, beta_ex[, colnames(beta_rs)])

# feature table: which methylation probe reports which SNP-array SNP
feat <- rbind(
  data.frame(feature = intersect(rownames(beta_rs), rownames(dos)),
             rsid = intersect(rownames(beta_rs), rownames(dos)), kind = "rs"),
  unique(as.data.frame(fread("ref/450k_cpgsnp_snp6.tsv")[, .(feature = probe, rsid, kind = where)]))
)
feat <- feat[feat$feature %in% rownames(beta) & feat$rsid %in% rownames(dos), ]
message(sprintf("candidate features: %d probes over %d SNPs (%s)",
                nrow(feat), length(unique(feat$rsid)),
                paste(names(table(feat$kind)), table(feat$kind), collapse = ", ")))

## Seeds from the rs-only matcher
seed_fit <- match_meth_snp(beta_rs[intersect(rownames(beta_rs), rownames(dos)), ], dos)
sm <- seed_fit$matches
seeds <- data.frame(meth = sm$meth_sample, snp = sm$best_snp_sample)[sm$best_score <= 0.09 & sm$gap >= 0.03, ]
message(sprintf("seed pairs: %d", nrow(seeds)))

m <- match_lr(beta, dos, feat, seeds)
saveRDS(m, "work/northcott_cavalli_lr.rds")

res <- as.data.table(m$matches)
res[, margin := best_llr - second_llr]

## Independent checks the matcher never saw: subgroup, sex, SNP-array sex
norm_sg <- function(x) sub("^Group ?", "Group", x)
sx <- fread("work/cavalli_sex.tsv")[, .(meth_sample = sub("_.*", "", array), meth_sex = predictedSex)]
snp_sx <- rbindlist(lapply(list.files("work", "summary.tsv.gz$", full.names = TRUE), fread),
                    fill = TRUE)[, .(best_snp_sample = sub("_.*", "", sample),
                                     snp_array_sex = ifelse(computed_gender == "female", "F", "M"))]
res <- merge(res, sx, by = "meth_sample")
res <- merge(res, data.table(meth_sample = meta_m$gsm, meth_title = meta_m$title,
                             meth_subgroup = norm_sg(meta_m$subgroup)), by = "meth_sample")
res <- merge(res, data.table(best_snp_sample = meta_s$gsm, snp_title = meta_s$title,
                             snp_sex = meta_s$gender, snp_age = meta_s$age,
                             snp_subgroup = norm_sg(meta_s$subgroup)), by = "best_snp_sample")
res <- merge(res, unique(snp_sx), by = "best_snp_sample", all.x = TRUE)
res[, sex_agree := meth_sex == snp_array_sex]
res[, subgroup_agree := meth_subgroup == snp_subgroup]
res <- res[order(-best_llr)]
fwrite(res, "work/northcott_cavalli_lr_matches.tsv", sep = "\t")

cat(sprintf("\nfeatures used: %d (median R2 %.2f)\n", nrow(m$features), median(m$features$r2)))
cat("SNPs compared per pair:", paste(range(res$n_snps), collapse = "-"), "\n")

null <- m$llr
null[cbind(match(res$meth_sample, rownames(null)), match(res$best_snp_sample, colnames(null)))] <- NA
cat("\nlog10 LR, every pair except each sample's own best:\n")
print(round(quantile(null, c(0.5, 0.99, 0.9999, 1), na.rm = TRUE), 1))
cat("\nlog10 LR of best match, deciles:\n")
print(round(quantile(res$best_llr, seq(0, 1, 0.1)), 1))

cat("\nagreement with held-out annotation, by log10 LR of best match:\n")
res[, bin := cut(best_llr, c(-Inf, 0, 5, 10, 20, 40, Inf))]
print(res[, .(n = .N, sex_agree = round(mean(sex_agree, na.rm = TRUE), 3),
              subgroup_agree = round(mean(subgroup_agree, na.rm = TRUE), 3),
              median_margin = round(median(margin), 1)), by = bin][order(bin)])

cat("\ncalled matches (LR > 10 and margin > 5):", res[best_llr > 10 & margin > 5, .N], "\n")
cat("SNP-array samples claimed by more than one methylation sample:",
    res[best_llr > 10 & margin > 5, .N, by = best_snp_sample][N > 1, .N], "\n")
