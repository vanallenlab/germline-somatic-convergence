#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Prioritize a list of candidate negative control phenotypes for permutation testing

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)

# Inputs
args <- c("/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/data/negative_control_phenos/data-2025-01-10.csv",
          "/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/data/negative_control_phenos/ukb31063_h2_topline.02Oct2019.tsv.gz")

# Load genetic correlations
rg <- read.table(args[1], sep=",", header=T)

# Helper function for collapsing rg results across a vector of codes
collapse.rg <- function(cancer, sub.codes){
  rg.sub <- rg[which(rg$ID1 %in% sub.codes | rg$ID2 %in% sub.codes), ]
  other.codes <- setdiff(as.character(unlist(rg.sub[, c("ID1", "ID2")])), sub.codes)
  as.data.frame(do.call("rbind", lapply(other.codes, function(oc){
    oc.sub <- rg.sub[which(rg.sub$ID1 %in% oc | rg.sub$ID2 %in% oc),
                     c("ID1", "ID2", "rpheno", "rg", "rg.SE", "Z", "p")]
    oc.sub <- oc.sub[!duplicated(oc.sub[, -c(1:2)]), ]
    best.row <- as.vector(oc.sub[which.min(oc.sub$p), ])
    # Find and retain most significant rg
    c("cancer" = cancer,
      "cancer_code" = as.character(setdiff(best.row[c("ID1", "ID2")], oc))[1],
      "other_pheno" = oc,
      "rg" = best.row$rg,
      "p" = best.row$p,
      "rpheno" = best.row$rpheno)
  })))
}

# Only keep pairs of traits that are genetically uncorrelated with any cancer type
pancan.codes <- c("134", "2453", "40009")
pancan.rg <- collapse.rg("pancan", pancan.codes)
rg.numeric.cols <- c("rg", "p", "rpheno")
pancan.rg[, rg.numeric.cols] <- apply(pancan.rg[, rg.numeric.cols], 2, as.numeric)
pancan.rg <- pancan.rg[which(pancan.rg$p > 0.25 & abs(pancan.rg$rg) < 0.05 & abs(pancan.rg$rpheno) < 0.05), ]
elig.phenos <- unique(pancan.rg$other_pheno)
rg <- rg[which(rg$ID1 %in% elig.phenos | rg$ID2 %in% elig.phenos), ]

# Subset to cancer-related pairs (but not cancer-cancer pairs)
prostate.codes <- c("20001_1044", "40001_C61", "C3_PROSTATE", "C61", "C_PROSTATE")
cancer.codes <- list("breast" = c("20001_1002", "C3_BREAST_3", "C50", "C_BREAST_3"),
                     "colorectal" = c("20001_1020", "20001_1022", "20001_1023", "20002_1460", "40001_C19", "40001_C20", "40001_C189", "C3_COLON", "C3_RECTOISIGMOID_JUNCTION", "C3_RECTUM", "C18", "C19", "C20", "C_COLON", "C_RECTOISIGMOID_JUNCTION", "C_RECTUM", "D12", "C3_DIGESTIVE_ORGANS"),
                     "lung" = c("20001_1001", "20107_3", "20110_3", "22140", "C3_BRONCHUS_LUNG", "C34", "C_BRONCHUS_LUNG", "LUNG_CANCER", "LUNG_CANCER_MESOT"),
                     "prostate" = prostate.codes,
                     "renal" = c("20001_1034", "20001_1035", "20002_1192", "20002_1193", "20002_1194", "20002_1405", "20002_1427", "40001_C64", "C3_BLADDER", "C3_KIDNEY_NOTRENALPELVIS", "C64", "C67", "C_BLADDER", "C_KIDNEY_NOTRENALPELVIS", "E27", "I12", "N17", "N18", "N19", "N28", prostate.codes))
all.cancer.codes <- unique(unlist(cancer.codes))
cancer.pairs <- which(rg$ID1 %in% all.cancer.codes | rg$ID2 %in% all.cancer.codes)
intercancer.pairs <- which(rg$ID1 %in% all.cancer.codes & rg$ID2 %in% all.cancer.codes)
rg <- rg[setdiff(cancer.pairs, intercancer.pairs), ]

# Simplify rg results by collapsing per cancer type
rg <- as.data.frame(do.call("rbind", lapply(names(cancer.codes), function(cancer){
  sub.codes <- cancer.codes[[cancer]]
  collapse.rg(cancer, sub.codes)
})))
rg[, rg.numeric.cols] <- apply(rg[, rg.numeric.cols], 2, as.numeric)

# Filter based on genetic correlations
rg <- rg[which(rg$p > 0.5 & abs(rg$rg) < 0.05 & abs(rg$rpheno) < 0.05), ]

# Load heritability estimates
h2g <- read.table(args[2], sep="\t", comment.char="", header=T, quote="")

# Filter heritability data
h2g <- h2g[which(h2g$h2_liability > 0.05 & h2g$h2_p < 0.05
                 & h2g$confidence %in% c("high", "medium")
                 & !h2g$isBadPower), ]

# Subset heritability data to columns of interest
h2g <- h2g[, c("phenotype", "description", "h2_liability", "h2_z", "confidence", "Neff")]

# Intersect rg & h2g data
dat <- merge(rg, h2g, by.x="other_pheno", by.y="phenotype", all=F, sort=F)
numeric.cols <- c("rg", "p", "rpheno", "h2_liability", "h2_z", "Neff")
dat[, numeric.cols] <- apply(dat[, numeric.cols], 2, as.numeric)

# For each cancer, rank-order candidate phenotypes based on minimizing |rg|, P(rg), and |rpheno|
candidates <- lapply(names(cancer.codes), function(cancer){
  dat.sub <- dat[which(dat$cancer == cancer), ]
  dat.sub$rg_rank <- order(abs(dat.sub$rg))
  dat.sub$p_rank <- order(-dat.sub$p)
  dat.sub$rpheno_rank <- order(abs(dat.sub$rpheno))
  dat.sub$mean_rank <- apply(dat.sub[, grep("_rank$", colnames(dat.sub))], 1, mean)
  dat.sub[order(dat.sub$mean_rank), ]
})
names(candidates) <- names(cancer.codes)

# candidates$breast$description
# candidates$colorectal$description
# candidates$lung$description
candidates$prostate$description
