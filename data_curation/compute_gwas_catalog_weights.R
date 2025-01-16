#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2024
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Compute gene weights for noncoding and coding GWAS hits

# Set parameters
options(scipen=1000, stringsAsFactors=F)

# Takes exactly two positional arguments:
# 1. Filtered & annotated GWAS catalog loci
# 2. List of eligible gene symbols
# 3. Output prefix for gene weights; will be suffixed by "_coding.tsv" and "_noncoding.tsv"
args <- commandArgs(trailingOnly=TRUE)

# # Dev parameters
# args <- c("/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.filtered.annotated.tsv.gz",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gencode.v47.autosomal.protein_coding.genes.list",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.gwas_catalog")

# Read GWAS loci
loci <- read.table(args[1], header=T, sep="\t", comment.char="", quote="")

# Read list of eligible genes
all.genes <- read.table(args[2], header=F)[, 1]

# Function to compute weights for a subset of loci
calc.weights <- function(loci, all.genes){
  g <- unlist(strsplit(loci$MAPPED_GENE, split="(, | - |; | x )"))
  g <- g[which(g %in% all.genes)]
  g <- c(g, all.genes)
  weights <- as.data.frame(table(g))
  colnames(weights) <- c("#gene", "gwas_hits")
  return(weights)
}

# Compute coding & noncoding weights
w.c <- calc.weights(loci[which(loci$gene_context == "coding_exon"), ], all.genes)
w.nc <- calc.weights(loci[which(loci$gene_context != "coding_exon"), ], all.genes)

# Write weights to output files
write.table(w.c, paste(args[3], "coding.tsv", sep="_"),
            col.names=T, row.names=F, quote=F, sep="\t")
write.table(w.nc, paste(args[3], "noncoding.tsv", sep="_"),
            col.names=T, row.names=F, quote=F, sep="\t")
