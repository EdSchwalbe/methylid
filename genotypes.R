# Genotype loaders producing a dosage matrix (SNPs x samples), values 0/1/2/NA.
# Rows are named by rsID, to match the feature table from build_reference.R.
suppressMessages(library(data.table))

# PLINK 1 binary. Only the requested variants are read, by seeking, so this is
# usable on a full GSA or imputed dataset without loading it all.
read_plink_dosage <- function(prefix, rsids = NULL) {
  bim <- fread(paste0(prefix, ".bim"), header = FALSE,
               col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  fam <- fread(paste0(prefix, ".fam"), header = FALSE)
  # prefer the IID alone; fall back to FID_IID only if IIDs are not unique
  ids <- as.character(if (ncol(fam) > 1) fam[[2]] else fam[[1]])
  if (anyDuplicated(ids) && ncol(fam) > 1) ids <- paste(fam[[1]], fam[[2]], sep = "_")
  n <- length(ids)
  bim[, idx := .I]
  want <- if (is.null(rsids)) bim else bim[rsid %in% rsids]
  want <- want[!duplicated(rsid)]
  if (!nrow(want)) stop("none of the requested rsIDs are in ", prefix, ".bim")

  bytes_per <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 3)
  if (!identical(as.integer(magic), c(0x6cL, 0x1bL, 0x01L)))
    stop("not a SNP-major PLINK 1 .bed file")

  # 2 bits per sample: 00 = hom a1, 01 = missing, 10 = het, 11 = hom a2.
  # Dosage counts a2, the allele PLINK lists second.
  lut <- c(0L, NA_integer_, 1L, 2L)
  dos <- matrix(NA_integer_, nrow(want), n, dimnames = list(want$rsid, ids))
  for (k in seq_len(nrow(want))) {
    seek(con, 3 + (want$idx[k] - 1) * bytes_per)
    raw <- readBin(con, "raw", bytes_per)
    bits <- as.integer(rawToBits(raw))
    code <- bits[seq(1, length(bits), 2)] + 2L * bits[seq(2, length(bits), 2)]
    dos[k, ] <- lut[code[seq_len(n)] + 1L]
  }
  message(sprintf("PLINK: %d of %d requested SNPs x %d samples, call rate %.1f%%",
                  nrow(dos), length(unique(rsids)), n, 100 * mean(!is.na(dos))))
  dos
}

# VCF (plain or gzipped) for sequencing data. Reads the file once, keeping only
# the requested sites, so pre-subset a WGS VCF if large.
#   rsids : marker rsIDs, matched against the ID column
#   sites : optional data.frame(rsid, chrom, pos) so markers can be found by
#           position when the ID column is "." (which is the norm for release
#           and pipeline VCFs). fetch_positions.py builds this table.
# Dosage counts non-reference alleles, so at a multi-allelic site 1/2 scores 2.
# Allele orientation is learned downstream, so REF/ALT conventions do not
# matter here, but a multi-allelic marker is noisier and is reported as such.
read_vcf_dosage <- function(path, rsids = NULL, sites = NULL, min_dp = 0) {
  site_key <- NULL
  if (!is.null(sites)) {
    sites <- as.data.frame(sites)
    site_key <- setNames(as.character(sites$rsid),
                         paste0(sub("^chr", "", sites$chrom), ":", sites$pos))
    site_key <- site_key[!duplicated(names(site_key))]
  }
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  samples <- NULL
  keep <- list()
  repeat {
    ln <- readLines(con, n = 20000)
    if (!length(ln)) break
    hdr <- grep("^#CHROM", ln)
    if (length(hdr)) samples <- strsplit(ln[hdr], "\t")[[1]][-(1:9)]
    ln <- ln[!startsWith(ln, "#")]
    if (!length(ln)) next
    f <- tstrsplit(ln, "\t", fixed = TRUE)
    # Join on the ID column where it carries rsIDs, and on chrom:pos otherwise.
    # The position route is not a nicety: the 1000 Genomes phase 3 release VCFs
    # have "." in ID for every variant, and pipeline output often does too.
    by_id <- if (is.null(rsids)) rep(TRUE, length(ln)) else f[[3]] %in% rsids
    by_pos <- rep(FALSE, length(ln))
    if (!is.null(site_key)) {
      k <- paste0(sub("^chr", "", f[[1]]), ":", f[[2]])
      by_pos <- k %in% names(site_key)
    }
    sel <- by_id | by_pos
    if (any(sel)) {
      g <- lapply(f, `[`, sel)
      # row name: the rsID if the file gave one, else the one we looked up
      nm <- g[[3]]
      if (!is.null(site_key)) {
        lk <- site_key[paste0(sub("^chr", "", g[[1]]), ":", g[[2]])]
        nm <- ifelse(is.na(nm) | nm == "." | !nm %in% rsids, lk, nm)
      }
      g[[3]] <- nm
      keep[[length(keep) + 1L]] <- g
    }
  }
  if (!length(keep)) stop("no requested sites found in ", path)
  f <- lapply(seq_along(keep[[1]]), function(i) unlist(lapply(keep, `[[`, i)))
  ok <- !is.na(f[[3]])
  f <- lapply(f, `[`, ok)
  fmt <- f[[9]]
  gt_field <- vapply(strsplit(fmt, ":", fixed = TRUE), function(x) match("GT", x), 1L)
  dp_field <- vapply(strsplit(fmt, ":", fixed = TRUE), function(x) match("DP", x), 1L)

  parse_gt <- function(col) {
    parts <- strsplit(col, ":", fixed = TRUE)
    gt <- mapply(function(p, i) if (length(p) >= i) p[i] else NA_character_, parts, gt_field)
    # count non-reference alleles; handles / and |, multi-allelic codes, and
    # haploid calls ("0"/"1"), which appear on chrX and chrY in male samples
    a <- strsplit(gt, "[/|]")
    d <- vapply(a, function(x) {
      x <- x[x != "."]
      if (!length(x)) return(NA_integer_)
      n <- suppressWarnings(as.integer(x))
      if (anyNA(n)) return(NA_integer_)
      if (length(n) == 1L) return(if (n[1] > 0L) 2L else 0L)
      sum(n > 0L)
    }, 1L)
    if (min_dp > 0) {
      dp <- suppressWarnings(as.integer(mapply(
        function(p, i) if (!is.na(i) && length(p) >= i) p[i] else NA, parts, dp_field)))
      d[!is.na(dp) & dp < min_dp] <- NA_integer_
    }
    d
  }
  dos <- do.call(cbind, lapply(f[-(1:9)], parse_gt))
  dimnames(dos) <- list(f[[3]], samples)
  dos <- dos[!duplicated(rownames(dos)), , drop = FALSE]
  multi <- sum(grepl(",", f[[5]], fixed = TRUE))
  message(sprintf("VCF: %d SNPs x %d samples, call rate %.1f%%%s",
                  nrow(dos), ncol(dos), 100 * mean(!is.na(dos)),
                  if (multi) sprintf(" (%d multi-allelic)", multi) else ""))
  dos
}

# Affymetrix SNP6 Birdseed CHPs, via read_calvin_chp.py + the pd annotation.
read_snp6_dosage <- function(calls_files, probeset_ann) {
  calls <- rbindlist(lapply(calls_files, fread))
  calls[, gsm := sub("_.*", "", sample)]
  calls[, dosage := c(`6` = 0, `8` = 1, `7` = 2)[as.character(call)]]
  calls <- merge(calls, probeset_ann, by = "probeset")
  calls <- calls[order(confidence)][!duplicated(calls[order(confidence), .(gsm, rsid)])]
  as.matrix(dcast(calls, rsid ~ gsm, value.var = "dosage"), rownames = "rsid")
}
