#!/usr/bin/env Rscript

# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Curate gene-specific nonsynonymous mutation rates from gnomAD

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)

# Read command-line arguments
args <- commandArgs(trailingOnly=TRUE)
in.tsv <- as.character(args[1])
out.tsv <- as.character(args[2])

# Read & clean gnomAD mutation rate data
x <- read.table(in.tsv, header=T, sep="\t", comment.char="")
x$mu <- x$mu_lof + x$mu_mis
x <- x[which(!is.na(x$mu)), c("gene", "mu")]
x <- x[order(x$gene), ]

# Write to out.tsv
colnames(x)[1] <- "#gene"
write.table(x, out.tsv, sep="\t", quote=F, col.names=T, row.names=F)
