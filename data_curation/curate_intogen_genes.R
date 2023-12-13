#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Code to curate IntOGen somatic driver gene lists

# Set parameters
options(scipen=1000, stringsAsFactors=F)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023")
intogen.in <- "other_data/2023-05-31_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0){
  outdir <- as.character(args[1])
}else{
  outdir <- "other_data/intogen_genes/"
}
if(!dir.exists(outdir)){
  dir.create(outdir)
}


# Read IntOGen driver compendium
idat <- read.table(intogen.in, sep="\t", header=T)

# Only retain genes with FDR Q<0.05
idat <- idat[which(idat$QVALUE_COMBINATION < 0.05), ]

# Build map of eligible cancer terms
# This must be done manually from IntoGen metadata
cancer.map <- list("Breast" = c("BRCA"),
                   "Colorectal" = c("COADREAD", "COAD"),
                   "Lung" = c("NSCLC", "SCLC", "LUAD", "LUSC"),
                   "Melanoma" = c("MEL", "UM"),
                   "Prostate" = c("PRAD"),
                   "Renal" = c("CCRCC", "RCC", "PRCC", "CHRCC"))

# Annotate IntOGen data with simple true/false columns for each cancer type
for(cancer in names(cancer.map)){
  cterms <- cancer.map[[cancer]]
  somatic.labs <- idat$CANCER_TYPE %in% cterms
  idat[, paste(cancer, "somatic", sep=".")] <- somatic.labs
}

# Write lists of germline & somatic drivers to gene list files
for(cancer in names(cancer.map)){
    outdir.sub <- paste(outdir, "/somatic_coding", sep="")
    if(!dir.exists(outdir.sub)){
      dir.create(outdir.sub)
    }
    genes <- sort(unique(idat$SYMBOL[which(idat[, paste(cancer, "somatic", sep=".")])]))
    f.out <- paste(outdir.sub, "/", tolower(cancer),
                   ".somatic.coding.genes.list", sep="")
    write.table(genes, f.out, col.names=F, row.names=F, quote=F)
}

