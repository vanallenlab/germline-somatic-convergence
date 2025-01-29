#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Compute composite permutation weights by multiplying two sets of gene weights
# Normalizes each set of weights first to ensure equivalent scaling

# Set parameters
options(scipen=1000, stringsAsFactors=F)

# Takes three positional arguments: weights_a, weights_b, and output path
args <- commandArgs(trailingOnly=TRUE)

# # DEV:
# args <- c("/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.composite_germline_noncoding.tsv",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.ppi_connections.tsv",
#           "~/scratch/composite_weights.test.tsv")

# Read & normalize weights A
weights.a <- read.table(args[1], header=T, sep="\t", comment.char="")
weights.a$prob <- as.numeric(weights.a[, 2]) / sum(as.numeric(weights.a[, 2]))

# Read & normalize weights A
weights.b <- read.table(args[2], header=T, sep="\t", comment.char="")
weights.b$prob <- as.numeric(weights.b[, 2]) / sum(as.numeric(weights.b[, 2]))

# Composite weights = P(A_g) x P(B_g)
w.m <- merge(weights.a, weights.b, by=1, all=F, sort=F, suffixes=c(".a", ".b"))
w.m$prob <- w.m$prob.a * w.m$prob.b
w.m$prob <- w.m$prob / sum(w.m$prob)
colnames(w.m)[1] <- "#gene"

# Write composite weights to file
write.table(w.m[, c("#gene", "prob")], args[3], col.names=T, row.names=F, sep="\t", quote=F)

