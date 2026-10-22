# Match Cavalli 2017 450K (GSE85212) to Northcott 2012 SNP6 (GSE37384).
# Run from the repository root: Rscript exploration/run_northcott_cavalli.R
# Needs the intermediate files in work/ (see HANDOVER.md); superseded by
# match_samples.R, kept to document how the original GEO analysis was done.
root <- dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
suppressMessages(library(data.table))
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

# SNP6 dosage matrix (allele B count), autosomes only
ann <- fread("ref/snp6_rs_overlap.tsv")
calls <- rbindlist(lapply(list.files("work", "^chpchunk_.*tsv.gz$", full.names = TRUE), fread))
calls[, gsm := sub("_.*", "", sample)]
calls[, dosage := c(`6` = 0, `8` = 1, `7` = 2)[as.character(call)]]
calls <- merge(calls, ann[, .(probeset = man_fsetid, rsid = dbsnp_rs_id, chrom)], by = "probeset")
calls <- calls[chrom %in% as.character(1:22)]
dos <- as.matrix(dcast(calls, rsid ~ gsm, value.var = "dosage"), rownames = "rsid")

beta <- as.matrix(fread("work/cavalli_rs_beta.tsv"), rownames = 1)
colnames(beta) <- sub("_.*", "", colnames(beta))

m <- match_meth_snp(beta, dos)
res <- m$matches
saveRDS(m, "work/northcott_cavalli_match.rds")

# Independent checks: subgroup and sex
sex_m <- fread("work/cavalli_sex.tsv")[, .(meth_sample = sub("_.*", "", array), meth_sex = predictedSex)]
norm_sg <- function(x) sub("^Group ?", "Group", x)
res <- merge(res, sex_m, by = "meth_sample")
res <- merge(res, data.frame(meth_sample = meta_m$gsm, meth_title = meta_m$title,
                             meth_subgroup = norm_sg(meta_m$subgroup)), by = "meth_sample")
res <- merge(res, data.frame(best_snp_sample = meta_s$gsm, snp_title = meta_s$title,
                             snp_sex = meta_s$gender, snp_age = meta_s$age,
                             snp_subgroup = norm_sg(meta_s$subgroup)), by = "best_snp_sample")
res$sex_agree <- res$meth_sex == res$snp_sex
res$subgroup_agree <- res$meth_subgroup == res$snp_subgroup
res <- res[order(res$best_score), ]
fwrite(res, "work/northcott_cavalli_matches.tsv", sep = "\t")

# Null: all non-best pairs
null <- m$score
null[cbind(seq_len(nrow(null)), match(res$best_snp_sample[match(rownames(null), res$meth_sample)], colnames(null)))] <- NA
cat("\nScore distribution, all pairs except each sample's best:\n"); print(quantile(null, c(0, 0.0001, 0.001, 0.01, 0.5), na.rm = TRUE))
cat("\nBest-match score distribution:\n"); print(quantile(res$best_score, seq(0, 1, 0.1)))
cat("\nOrientation flips:", sum(m$orientation$flipped), "of", nrow(m$orientation), "\n")
cat("SNPs per pair:", range(res$n_snps), "\n")
br <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 1)
tab <- with(res, data.frame(bin = cut(best_score, br), sex_agree, subgroup_agree, gap, opp_hom))
print(aggregate(cbind(n = 1, sex_agree, subgroup_agree, gap, opp_hom) ~ bin, tab,
                function(x) round(c(mean = mean(x), n = length(x))[if (identical(x, rep(1, length(x)))) 2 else 1], 3)))
