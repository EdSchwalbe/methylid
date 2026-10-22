# Extract rs-probe betas and inferred sex from 450K/EPIC IDATs.
# Usage: Rscript extract_meth_snps.R IDAT_DIR OUT_PREFIX [NCORES] [EXTRA_PROBES]
# EXTRA_PROBES: optional file of CpG probe IDs (e.g. probes with a SNP at the CpG
# or extension base) whose raw betas are also written, to OUT_PREFIX_extra_beta.tsv
#
# Writes OUT_PREFIX_beta_all.tsv -- the rs and extra probes in one matrix, which
# is what match_samples.R takes as --beta.
suppressMessages({
  library(minfi)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
idat_dir <- args[1]
out_prefix <- args[2]
ncores <- if (length(args) > 2) as.integer(args[3]) else 16
extra <- if (length(args) > 3) readLines(args[4]) else character(0)

bases <- unique(sub("_(Grn|Red)\\.idat(\\.gz)?$", "",
                    list.files(idat_dir, pattern = "_Grn\\.idat(\\.gz)?$", full.names = TRUE)))
message(length(bases), " arrays")

# Read in chunks in parallel, keeping only what we need from each chunk.
# Never make more chunks than arrays, or a chunk ends up empty.
nchunk <- max(1, min(ncores, length(bases)))
chunks <- split(bases, cut(seq_along(bases), nchunk, labels = FALSE))
res <- bplapply(chunks, function(b) {
  rg <- read.metharray(b, force = TRUE)
  snp <- getSnpBeta(rg)
  raw <- preprocessRaw(rg)
  # getSex clusters the X/Y intensities, so it fails on a handful of arrays or
  # when they are all one sex. Sex is only ever an independent check here, so a
  # failure must not cost us the genotypes: fall back to the raw medians.
  sex <- tryCatch(getSex(mapToGenome(raw)), error = function(e) {
    message("getSex failed on this chunk (", conditionMessage(e),
            "); reporting X/Y medians without a call")
    gr <- mapToGenome(raw)
    cn <- getCN(gr)
    chr <- as.character(GenomicRanges::seqnames(gr))
    data.frame(predictedSex = NA_character_,
               xMed = apply(cn[chr == "chrX", , drop = FALSE], 2, median, na.rm = TRUE),
               yMed = apply(cn[chr == "chrY", , drop = FALSE], 2, median, na.rm = TRUE))
  })
  list(snp = snp, extra = getBeta(raw)[intersect(extra, rownames(raw)), , drop = FALSE],
       sex = data.frame(array = colnames(rg), predictedSex = sex$predictedSex,
                        xMed = sex$xMed, yMed = sex$yMed))
}, BPPARAM = MulticoreParam(nchunk))

snp <- do.call(cbind, lapply(res, `[[`, "snp"))
sex <- do.call(rbind, lapply(res, `[[`, "sex"))
colnames(snp) <- basename(colnames(snp))
ex <- NULL
if (length(extra)) {
  ex <- do.call(cbind, lapply(res, `[[`, "extra"))
  colnames(ex) <- basename(colnames(ex))
  write.table(ex, paste0(out_prefix, "_extra_beta.tsv"), sep = "\t", quote = FALSE, col.names = NA)
}

# the combined matrix match_samples.R expects
all_beta <- if (is.null(ex)) snp else rbind(snp, ex[, colnames(snp), drop = FALSE])
write.table(all_beta, paste0(out_prefix, "_beta_all.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)
sex$array <- basename(sex$array)

write.table(snp, paste0(out_prefix, "_rs_beta.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(sex, paste0(out_prefix, "_sex.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("done: ", ncol(all_beta), " samples x ", nrow(all_beta),
        " probes (", nrow(snp), " rs + ", nrow(all_beta) - nrow(snp), " SNP-affected)")
message("  -> ", out_prefix, "_beta_all.tsv  (pass this as --beta)")
message("  -> ", out_prefix, "_sex.tsv       (pass this as --sex)")
