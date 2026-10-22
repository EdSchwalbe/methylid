# Likelihood-ratio matching of methylation samples to SNP-array samples.
#
# Each feature is a methylation probe whose beta depends on a SNP genotyped on
# the SNP array (the rs probes, plus CpG probes with a common SNP at the CpG or
# single-base-extension position). For each feature, beta ~ dosage is fitted on
# training pairs (known or confidently matched), and features with R^2 below
# min_r2 are dropped, keeping the best feature per SNP. A pair (i, j) is scored
#   log10 LR = sum_f log10 [ p(beta_if | g_jf) / sum_g freq_f(g) p(beta_if | g) ]
# with p a normal mixed with a small uniform outlier component, so that LOH or
# copy-number driven allelic imbalance in tumours cannot veto a true match.

fit_features <- function(beta, dos, feat, pairs, min_r2 = 0.5, min_pairs = 30) {
  fits <- lapply(seq_len(nrow(feat)), function(k) {
    b <- beta[feat$feature[k], pairs$meth]
    g <- dos[feat$rsid[k], pairs$snp]
    ok <- !is.na(b) & !is.na(g)
    if (sum(ok) < min_pairs || length(unique(g[ok])) < 2) return(NULL)
    f <- lm(b[ok] ~ g[ok])
    r2 <- summary(f)$r.squared
    data.frame(feature = feat$feature[k], rsid = feat$rsid[k], a = coef(f)[1], b = coef(f)[2],
               sd = max(sigma(f), 0.02), r2 = r2)
  })
  fits <- do.call(rbind, fits)
  fits <- fits[fits$r2 >= min_r2, ]
  fits <- fits[order(-fits$r2), ]
  fits[!duplicated(fits$rsid), ]
}

score_lr <- function(beta, dos, fits, outlier = 0.02) {
  B <- beta[fits$feature, , drop = FALSE]
  G <- dos[fits$rsid, , drop = FALSE]
  freq <- sapply(0:2, function(g) rowMeans(G == g, na.rm = TRUE))
  dens <- lapply(0:2, function(g) {
    (1 - outlier) * dnorm(B, fits$a + fits$b * g, fits$sd) + outlier
  })
  denom <- log10(Reduce(`+`, Map(function(d, k) d * freq[, k], dens, 1:3)))
  llr <- 0
  for (g in 0:2) {
    A <- log10(dens[[g + 1]]) - denom
    A[is.na(A)] <- 0
    I <- (G == g) * 1
    I[is.na(I)] <- 0
    llr <- llr + crossprod(A, I)
  }
  n <- crossprod(!is.na(B) * 1, !is.na(G) * 1)
  list(llr = llr, n = n)
}

summarise_best <- function(llr, n) {
  ord <- t(apply(llr, 1, order, decreasing = TRUE))
  i <- seq_len(nrow(llr))
  # The total LR of a true match grows with the number of SNPs compared, so an
  # absolute threshold does not port across platforms. Dividing by the number of
  # SNPs does: a true pair averages about +0.27 per SNP whatever the platform or
  # marker count, while an unrelated pair is negative once a few dozen SNPs are
  # compared. That per-SNP mean is the statistic to threshold on.
  lower <- t(apply(llr, 1, function(x) {
    x <- sort(x[is.finite(x)])
    x[seq_len(max(1, floor(length(x) * 0.8)))]
  }))
  mu <- apply(lower, 1, median)
  sdev <- pmax(apply(lower, 1, mad), 1e-6)
  data.frame(meth_sample = rownames(llr),
             best_snp_sample = colnames(llr)[ord[, 1]],
             best_llr = llr[cbind(i, ord[, 1])],
             second_snp_sample = colnames(llr)[ord[, 2]],
             second_llr = llr[cbind(i, ord[, 2])],
             n_snps = n[cbind(i, ord[, 1])],
             # the portable statistic: mean log10 LR per SNP compared
             per_snp = llr[cbind(i, ord[, 1])] / pmax(n[cbind(i, ord[, 1])], 1),
             second_per_snp = llr[cbind(i, ord[, 2])] /
               pmax(n[cbind(i, ord[, 2])], 1),
             # secondary diagnostic: standing out from this sample's own null.
             # Inflated by taking a maximum over many candidates, so it is
             # reported but not used to call matches.
             best_z = (llr[cbind(i, ord[, 1])] - mu) / sdev,
             second_z = (llr[cbind(i, ord[, 2])] - mu) / sdev)
}

# Cross-fitted matching: seeds are split into folds; each seed sample is scored
# with features fitted on the other folds; non-seed samples use all seeds.
match_lr <- function(beta, dos, feat, seeds, n_folds = 2, ...) {
  set.seed(1)
  fold <- sample(rep_len(seq_len(n_folds), nrow(seeds)))
  llr <- matrix(NA_real_, ncol(beta), ncol(dos), dimnames = list(colnames(beta), colnames(dos)))
  n <- llr
  fits_all <- fit_features(beta, dos, feat, seeds, ...)
  s <- score_lr(beta, dos, fits_all)
  llr[] <- s$llr
  n[] <- s$n
  for (k in seq_len(n_folds)) {
    fk <- fit_features(beta, dos, feat, seeds[fold != k, ], ...)
    test <- seeds$meth[fold == k]
    s <- score_lr(beta[, test, drop = FALSE], dos, fk)
    llr[test, ] <- s$llr
    n[test, ] <- s$n
    message(sprintf("fold %d: %d features from %d training pairs", k, nrow(fk), sum(fold != k)))
  }
  list(matches = summarise_best(llr, n), llr = llr, n = n, features = fits_all)
}
