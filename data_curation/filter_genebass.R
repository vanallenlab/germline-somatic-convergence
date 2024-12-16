#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Filter GeneBass coding variant associations for a single trait

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)
p.thresh <- 10e-6

# # DEV:
# args <- c("/Users/ryan/Partners\ HealthCare\ Dropbox/Ryan\ Collins/VanAllen/VALab_germline_somatic_2023/other_data/genebass/atrial_fibrilation.lof.genebass.csv",
#           "/Users/ryan/Partners\ HealthCare\ Dropbox/Ryan\ Collins/VanAllen/VALab_germline_somatic_2023/other_data/genebass/atrial_fibrilation.mislc.genebass.csv",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/gene_lists/germline_coding/atrial_fibrilation.germline.coding.genes.list")

# Read command-line arguments
args <- commandArgs(trailingOnly=TRUE)
lof.csv <- as.character(args[1])
mislc.csv <- as.character(args[2])
outfile <- as.character(args[3])

# Read & combine LoF + mis|LC associations
lof <- read.table(lof.csv, header=T, sep=",")
mislc <- read.table(mislc.csv, header=T, sep=",")
ss <- rbind(lof, mislc)

# Filter to significant genes and write to outfile
write.table(sort(unique(ss$Gene.Name[which(ss$P.Value.SKATO < p.thresh)])),
            outfile, quote=F, col.names=F, row.names=F)
