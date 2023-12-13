#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Code to curate GTEx gene expression matrixes for permutation weights

# Set parameters
options(scipen=1000, stringsAsFactors=F)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023")
gtex.in <- "other_data/gtex_median_expression_gencode_filtered_SH_20231205.csv"
outdir <- "other_data/permutation_weights"


# Read GTEx data
gtex <- read.table(gtex.in, sep=",", header=T, check.names=F)


# Get organ-level averages for organs with multiple biopsy sites
gtex$Colon <- apply(gtex[, grep("Colon", colnames(gtex))], 1, mean, na.rm=T)
gtex <- gtex[, -grep("Colon - ", colnames(gtex), fixed=T)]
gtex$Kidney <- apply(gtex[, grep("Kidney", colnames(gtex))], 1, mean, na.rm=T)
gtex <- gtex[, -grep("Kidney - ", colnames(gtex), fixed=T)]


# Deduplicate genes
gtex <- gtex[-which(duplicated(gtex$Gene)), ]


# Scale each tissue by the fraction of total RNA contributed by each gene
gtex[, -1] <- apply(gtex[, -1], 2, function(vals){
  # First, assign all genes a non-zero value
  zero.idxs <- which(vals == 0)
  if(length(zero.idxs) > 0){
    vals[zero.idxs] <- min(vals[-zero.idxs], na.rm=T)
  }

  # Second, scale all genes
  vals / sum(vals)
})

# Make matrix of weights for all tissues
t2c.map <- c("Breast - Mammary Tissue" = "breast",
             "Lung" = "lung",
             "Prostate" = "prostate",
             "Colon" = "colorectal",
             "Kidney" = "renal")
expression.weights <- do.call("rbind", lapply(colnames(gtex)[-1], function(tissue){
  sub.df <- gtex[, c(1, which(colnames(gtex) == tissue))]
  sub.df$cancer <- t2c.map[tissue]
  colnames(sub.df) <- c("gene", "weight", "cancer")
  return(sub.df)
}))

# Write expression weights to file
write.table(expression.weights,
            paste(outdir, "gene_weights.expression.cancer_specific.tsv", sep="/"),
            col.names=T, row.names=F, sep="\t", quote=F)
