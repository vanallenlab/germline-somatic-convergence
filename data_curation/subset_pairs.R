#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Subset convergent pairs to only pairs present in an input gene list

# Load libraries
require(RLCtools)

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
# # DEV
# args <- c("/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.union.tsv",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gencode.v47.autosomal.protein_coding.ebi_ppi_universe.genes.list",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.union.ebi_subset.tsv")
pairs.in <- as.character(args[1])
glist.in <- as.character(args[2])
pairs.out <- as.character(args[3])


# Read original pairs
pairs <- read.table(pairs.in, check.names=F, sep="\t", comment.char="", header=T)

# Read restriction gene list
glist <- read.table(glist.in, header=F)[, 1]

# Subset pairs to restriction list
pairs <- pairs[which(pairs$germline_gene %in% glist
                     & pairs$somatic_gene %in% glist), ]

# Write to output file
write.table(pairs, pairs.out, col.names=T, row.names=F, sep="\t", quote=F)
