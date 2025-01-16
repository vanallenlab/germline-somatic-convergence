#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Curate gene-specific nonsynonymous mutation rates from gnomAD

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)

# Read command-line arguments
args <- commandArgs(trailingOnly=TRUE)
# # Dev arguments
# args <- c("~/Desktop/Collins/VanAllen/germline_somatic_convergence/misc/gwas_perm_weight_investigation_jan2025/gnomad_metrics.tsv.gz",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gencode.v47.autosomal.protein_coding.genes.list",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.coding_nonsynonymous.tsv")
in.tsv <- as.character(args[1])
elig.genes.in <- as.character(args[2])
out.tsv <- as.character(args[3])

# Read & clean gnomAD mutation rate data
x <- read.table(in.tsv, header=T, sep="\t", comment.char="")
x$mu <- x$mu_lof + x$mu_mis
x <- x[which(!is.na(x$mu)), c("gene", "mu")]
x <- x[!duplicated(x$gene), ]

# Load eligible genes, filter gnomAD data, and impute missing genes as median
elig <- read.table(args[2], header=F)[, 1]
x <- x[which(x$gene %in% elig), ]
x <- as.data.frame(rbind(x, data.frame("gene"=elig[which(!elig %in% x$gene)], "mu"=median(x$mu, na.rm=T))))

# Write to out.tsv
x <- x[order(x$gene), ]
colnames(x)[1] <- "#gene"
write.table(x, out.tsv, sep="\t", quote=F, col.names=T, row.names=F)
