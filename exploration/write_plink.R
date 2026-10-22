# Write a dosage matrix as PLINK 1 binary, for testing the reader round-trip.
write_plink_dosage <- function(dos, prefix, chr = 1, a1 = "A", a2 = "G") {
  n <- ncol(dos)
  write.table(data.frame(chr, rownames(dos), 0, seq_len(nrow(dos)), a1, a2),
              paste0(prefix, ".bim"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  write.table(data.frame(colnames(dos), colnames(dos), 0, 0, 0, -9),
              paste0(prefix, ".fam"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  code <- c(`0` = 0L, `1` = 2L, `2` = 3L)  # hom a1, het, hom a2
  for (k in seq_len(nrow(dos))) {
    g <- dos[k, ]
    cc <- ifelse(is.na(g), 1L, code[as.character(g)])
    bits <- as.integer(rbind(cc %% 2L, cc %/% 2L))
    pad <- (-length(bits)) %% 8
    if (pad) bits <- c(bits, integer(pad))
    writeBin(packBits(as.integer(bits), "raw"), con)
  }
  invisible(prefix)
}
