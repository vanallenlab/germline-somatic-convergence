#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Postprocess a .tsv of permuted convergent pair counts

# Load libraries
require(RLCtools)

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
tsv.in <- as.character(args[1])
tsv.out <- as.character(args[2])

# # DEV:
# tsv.in <- "~/scratch/perm_test.results.tsv"
# tsv.out <- "~/scratch/perm_test.results.postprocessed.tsv"

# Make vector for expected column names
tiers <- paste("tier", 1:4, sep="")
origins <- c("germline", "somatic")
contexts <- c("coding", "noncoding")
strata.names <- as.vector(unlist(sapply(contexts, function(germline_context){
  sapply(contexts, function(somatic_context){
    paste("germline_", germline_context, ".somatic_", somatic_context, ".",
          tiers, sep="")
  })
})))

# Load permuted data as data.frame
df <- read.table(tsv.in, sep="\t", header=F,
                      col.names=c("cancer", "perm_idx.major", "perm_idx.minor",
                                  strata.names))

# Make summary columns for all contexts
for(tier in tiers){
  ac <- apply(df[, grep(tier, colnames(df))], 1, sum, na.rm=T)
  df[, paste("germline_any.somatic_any", tier, sep=".")] <- ac
}

# Make summary columns for multi-tier criteria
strata.bases <- gsub(".tier1", "", colnames(df)[grep("tier1", colnames(df))])
t2p <- df[, grep("tier1", colnames(df))] + df[, grep("tier2", colnames(df))]
df[, paste(strata.bases, "tier2plus", sep=".")] <- t2p
t3p <- df[, grep("tier2plus", colnames(df))] + df[, grep("tier3", colnames(df))]
df[, paste(strata.bases, "tier3plus", sep=".")] <- t3p
t4p <- df[, grep("tier3plus", colnames(df))] + df[, grep("tier4", colnames(df))]
df[, paste(strata.bases, "any", sep=".")] <- t4p

# Make summary rows for pooled result across all cancers
# and, separately, all negative control phenotype
prim.cancers <- sort(unique(df$cancer[grep("_", df$cancer, invert=T)]))
nc.phenos <- sort(unique(sapply(setdiff(df$cancer, prim.cancers), function(ps){
  paste(unlist(strsplit(ps, split="_"))[-1], collapse="_")
})))
perm.infos <- unique(df[, grep("perm_idx", colnames(df))])
append.df <- do.call("rbind", lapply(1:nrow(perm.infos), function(info.ridx){
  major.idx <- perm.infos[info.ridx, 1]
  minor.idx <- perm.infos[info.ridx, 2]
  sub.df <- df[which(df$perm_idx.major == major.idx
                     & df$perm_idx.minor == minor.idx), ]
  ac <- as.numeric(apply(sub.df[which(sub.df$cancer %in% prim.cancers), -(1:3)],
                         2, sum, na.rm=T))
  append.df <- as.data.frame(t(data.frame(c("any", major.idx, minor.idx, ac))))
  colnames(append.df) <- colnames(sub.df)
  for(pheno in nc.phenos){
    if(length(unique(sub.df$cancer[grep(pheno, sub.df$cancer)])) > 1){
      ac <- as.numeric(apply(sub.df[grep(pheno, sub.df$cancer), -(1:3)],
                             2, function(v){sum(as.numeric(v), na.rm=T)}))
      append.df <- as.data.frame(rbind(append.df, c(paste("any", pheno, sep="_"),
                                                    major.idx, minor.idx, ac)))
    }
  }
  any.nc.idx <- intersect(grep(paste(prim.cancers, collapse="|"), sub.df$cancer),
                          grep(paste(nc.phenos, collapse="|"), sub.df$cancer))
  ac <- as.numeric(apply(sub.df[any.nc.idx, -(1:3)], 2,
                         function(v){sum(as.numeric(v), na.rm=T)}))
  append.df <- as.data.frame(rbind(append.df, c("any_negative_control",
                                                major.idx, minor.idx, ac)))
  rownames(append.df) <- NULL
  return(append.df)
}))
df <- as.data.frame(rbind(df, as.data.frame(append.df)))

# Format permuted data and write to outfile
df[, -c(1:3)] <- apply(df[, -c(1:3)], 2, as.numeric)
write.table(df, tsv.out, col.names=T, row.names=F, quote=F, sep="\t")
