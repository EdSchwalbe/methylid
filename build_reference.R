# Build the probe -> SNP feature table for a methylation platform.
#
# Features are methylation probes whose beta tracks a genotype:
#   rs   : the dedicated SNP probes (rs* on 450K/EPIC, nv*/rs* on EPICv2)
#   CpG  : a common SNP at the interrogated C or G
#   SBE  : a common SNP at the single-base-extension position
# Only the probe -> rsID mapping is needed: which allele a probe's high beta
# reports is learned from the data later (the sign of the beta ~ dosage slope),
# so strand and A/B allele conventions never have to be reconciled here.
#
# Usage: build_reference.R PLATFORM OUT.tsv [MIN_MAF] [DBSNP_TABLE]
#   PLATFORM     450k | epic | epicv2
#   DBSNP_TABLE  optional: restrict to SNPs on your genotyping platform, as a
#                file with one rsID per line (e.g. cut -f2 your.bim)
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
platform <- tolower(args[1])
out_path <- args[2]
min_maf <- if (length(args) > 2) as.numeric(args[3]) else 0.05
keep_rs <- if (length(args) > 3) unique(readLines(args[4])) else NULL

snp_probe_ids <- function(platform) {
  if (platform == "epicv2") {
    suppressMessages(library(sesameData))
    addr <- sesameData::sesameDataGet("EPICv2.address")
    ids <- addr$ordering$Probe_ID
    return(grep("^(rs|nv)", ids, value = TRUE))
  }
  suppressMessages(library(minfi))
  man <- if (platform == "450k") {
    suppressMessages(library(IlluminaHumanMethylation450kmanifest))
    IlluminaHumanMethylation450kmanifest
  } else {
    suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
    IlluminaHumanMethylationEPICmanifest
  }
  c(minfi::getProbeInfo(man, type = "SnpI")$Name, minfi::getProbeInfo(man, type = "SnpII")$Name)
}

cpg_sbe_snps <- function(platform, min_maf) {
  anno <- switch(platform,
    "450k" = {
      suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
      IlluminaHumanMethylation450kanno.ilmn12.hg19
    },
    "epic" = {
      suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
      IlluminaHumanMethylationEPICanno.ilm10b4.hg19
    },
    "epicv2" = {
      suppressMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
      IlluminaHumanMethylationEPICv2anno.20a1.hg38
    })
  # newest dbSNP build shipped by the annotation package
  tbl <- grep("^SNPs\\.[0-9]+CommonSingle$", names(minfi::getAnnotationObject(anno)@data), value = TRUE)
  tbl <- tbl[order(as.integer(gsub("\\D", "", tbl)))]
  si <- as.data.frame(minfi::getAnnotation(anno, what = tail(tbl, 1)))
  message("using ", tail(tbl, 1))
  out <- rbind(
    data.table(feature = rownames(si), rsid = si$CpG_rs, maf = si$CpG_maf, kind = "CpG"),
    data.table(feature = rownames(si), rsid = si$SBE_rs, maf = si$SBE_maf, kind = "SBE")
  )[!is.na(rsid) & maf >= min_maf]

  # Note on EPICv2: it carries over 83% of EPIC v1 probes in general, but only
  # 0.4% of those with a common SNP at the interrogated CpG -- Illumina removed
  # the SNP-affected probes. So the EPIC v1 markers cannot be mapped across, and
  # EPICv2's usable set comes from its own (dbSNP 151, hg38) annotation, mostly
  # at different probes. Still thousands of SNPs, far more than the ~50 needed,
  # and the 65 rs probes of the 450K are restored. Which markers actually track
  # genotype is settled empirically by the R^2 filter during calibration.
  out
}

rs <- snp_probe_ids(platform)
message(sprintf("%s: %d dedicated SNP probes", platform, length(rs)))
# EPICv2 probe IDs carry a replicate suffix (rs1019916_BC21); the rsID is the
# part before it, and several probes may report the same SNP.
feat <- data.table(feature = rs, rsid = sub("_.*$", "", rs),
                   maf = NA_real_, kind = "rs")
feat <- feat[grepl("^rs", rsid)]

cs <- cpg_sbe_snps(platform, min_maf)
message(sprintf("%s: %d CpG/SBE probe-SNP pairs at MAF >= %.2f (%d SNPs)",
                platform, nrow(cs), min_maf, uniqueN(cs$rsid)))
feat <- rbind(feat, cs)

if (!is.null(keep_rs)) {
  before <- nrow(feat)
  feat <- feat[rsid %in% keep_rs]
  message(sprintf("restricted to genotyping platform: %d of %d pairs, %d SNPs",
                  nrow(feat), before, uniqueN(feat$rsid)))
}

fwrite(unique(feat), out_path, sep = "\t")
message("wrote ", out_path, ": ", nrow(unique(feat)), " probe-SNP pairs over ",
        uniqueN(feat$rsid), " SNPs")
