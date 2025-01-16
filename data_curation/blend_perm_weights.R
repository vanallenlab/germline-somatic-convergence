#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2024
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Scale & average two raw gene weight files

# Set parameters
options(scipen=1000, stringsAsFactors=F)

# Takes exactly two positional arguments:
# 1. First gene weight .tsv
# 2. Second gene weight .tsv
# 3. Output path to scaled & averaged composite gene weight .tsv
args <- commandArgs(trailingOnly=TRUE)

# # Dev parameters
# args <- c("/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.gwas_catalog_noncoding.tsv",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.sudmant_noncoding_snps.tsv",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.composite_germline_noncoding.tsv")

# Load raw weights
w1 <- read.table(args[1], header=T, sep="\t", comment.char="", check.names=F)
w2 <- read.table(args[2], header=T, sep="\t", comment.char="", check.names=F)

# Scale raw weights
w1$weight <- w1[, 2] / sum(w1[, 2], na.rm=T)
w2$weight <- w2[, 2] / sum(w2[, 2], na.rm=T)

# Merge scaled weights
m <- merge(w1[, c(colnames(w1)[1], "weight")],
           w2[, c(colnames(w2)[1], "weight")],
           by=colnames(w1)[1], sort=F, all=F, suffixes=c(".w1", ".w2"))
m$weight <- apply(m[, c("weight.w1", "weight.w2")], 1, mean, na.rm=T)

# Write scaled composite weights to output file
write.table(m[, c(colnames(m)[1], "weight")], args[3],
            col.names=T, row.names=F, sep="\t", quote=F)

