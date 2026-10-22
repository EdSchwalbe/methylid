# Match methylation-array samples to SNP-array samples using the rs probes.
#
# Inputs
#   meth_beta : rs-probe beta matrix (probes x samples), from extract_meth_snps.R
#   snp_calls : long table sample/rsid/dosage (count of allele B; NA = no call)
# Score per pair = mean |beta - dosage/2| over shared autosomal SNPs (0 = identical).
# Opposite homozygotes (beta < 0.2 vs BB, or > 0.8 vs AA) are counted separately:
# LOH and allelic imbalance can turn a het into a hom, but cannot produce these.
#
# Allele orientation (which allele the probe's high-beta signal reports) is learned
# from the data: an initial guess from allele frequencies, then refined from the
# signs of beta~dosage correlations across confident matches until stable.

match_meth_snp <- function(beta, dos, min_snps = 10, n_iter = 10, conf_gap = 0.1) {
  rs <- intersect(rownames(beta), rownames(dos))
  beta <- beta[rs, , drop = FALSE]
  dos <- dos[rs, , drop = FALSE]

  p_m <- rowMeans(beta, na.rm = TRUE)
  p_s <- rowMeans(dos, na.rm = TRUE) / 2
  flip <- abs(p_m - (1 - p_s)) < abs(p_m - p_s)

  score_all <- function(flip) {
    b <- beta
    b[flip, ] <- 1 - b[flip, ]
    g <- dos / 2
    ok_b <- !is.na(b)
    ok_g <- !is.na(g)
    b0 <- b; b0[!ok_b] <- 0
    g0 <- g; g0[!ok_g] <- 0
    n <- crossprod(ok_b * 1, ok_g * 1)
    # sum |b - g| expanded over the three genotype classes of g
    absdiff <- crossprod(b0, (ok_g & g0 == 0) * 1) +
      crossprod(abs(b0 - 0.5) * ok_b, (ok_g & g0 == 0.5) * 1) +
      crossprod((1 - b0) * ok_b, (ok_g & g0 == 1) * 1)
    opp <- crossprod((ok_b & b0 < 0.2) * 1, (ok_g & g0 == 1) * 1) +
      crossprod((ok_b & b0 > 0.8) * 1, (ok_g & g0 == 0) * 1)
    s <- absdiff / n
    s[n < min_snps] <- NA
    list(score = s, n = n, opp = opp)
  }

  for (it in seq_len(n_iter)) {
    sc <- score_all(flip)
    best <- apply(sc$score, 1, function(x) order(x)[1:2])
    gap <- sc$score[cbind(seq_len(ncol(beta)), best[2, ])] -
      sc$score[cbind(seq_len(ncol(beta)), best[1, ])]
    conf <- which(gap > conf_gap)
    if (length(conf) < 20) break
    r <- sapply(rs, function(k) cor(beta[k, conf], dos[k, best[1, conf]], use = "complete.obs"))
    new_flip <- r < 0
    message(sprintf("iter %d: %d confident matches, %d orientation changes",
                    it, length(conf), sum(new_flip != flip)))
    if (all(new_flip == flip)) break
    flip <- new_flip
  }

  s <- sc$score
  ord <- t(apply(s, 1, order))
  i <- seq_len(nrow(s))
  res <- data.frame(
    meth_sample = rownames(s),
    best_snp_sample = colnames(s)[ord[, 1]],
    best_score = s[cbind(i, ord[, 1])],
    second_snp_sample = colnames(s)[ord[, 2]],
    second_score = s[cbind(i, ord[, 2])],
    n_snps = sc$n[cbind(i, ord[, 1])],
    opp_hom = sc$opp[cbind(i, ord[, 1])]
  )
  res$gap <- res$second_score - res$best_score
  list(matches = res, score = s, n = sc$n, opp = sc$opp,
       orientation = data.frame(rsid = rs, flipped = flip))
}
