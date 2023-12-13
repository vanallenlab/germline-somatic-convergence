#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Code to compute composite permutation weights by multiplying one set of
# Bayesian priors with GTEx gene expression values

# Set parameters
options(scipen=1000, stringsAsFactors=F)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023")
args <- commandArgs(trailingOnly=TRUE)

# Takes three positional arguments: weights, expression weights, and output path
# # DEV:
# args <- c("other_data/permutation_weights/gene_weights.hrc_noncoding_snps.tsv",
#           "other_data/permutation_weights/gene_weights.expression.cancer_specific.tsv",
#           "other_data/permutation_weights/gene_weights.hrc_noncoding_snps.composite.cancer_specific.tsv")

# Read old weights
old.weights <- read.table(args[1], header=T, sep="\t", comment.char="")
old.genes <- unique(old.weights[, 1])

# Read expression data
gtex <- read.table(args[2], sep="\t", header=T, check.names=F, comment.char="")

# Iterate over each cancer type and multiply genes
weights.out <- do.call("rbind", lapply(unique(gtex$cancer), function(cancer){
  gtex.sub <- gtex[which(gtex$cancer == cancer), ]

  # Fill missing expression values for genes as the minimum of all non-zero genes
  missing.genes <- setdiff(old.genes, gtex.sub[, 1])
  n.missing <- length(missing.genes)
  gtex.missing <- data.frame(missing.genes,
                             rep(min(gtex.sub[, 2]), n.missing),
                             rep(cancer, n.missing))
  colnames(gtex.missing) <- colnames(gtex.sub)
  gtex.sub <- rbind(gtex.sub, gtex.missing)
  rownames(gtex.sub) <- gtex.sub[, 1]

  # Multiply weights
  new.sub <- gtex.sub[old.weights[, 1], ]
  new.sub[, 2] <- new.sub[, 2] * old.weights[, 2]
  rownames(new.sub) <- NULL
  new.sub
}))
colnames(weights.out) <- c("gene", "weight", "cancer")

# Write composite weights to file
write.table(weights.out, args[3], col.names=T, row.names=F, sep="\t", quote=F)
