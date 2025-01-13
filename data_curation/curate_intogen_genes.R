#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2023
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Code to curate IntOGen somatic driver gene lists

# Set parameters
options(scipen=1000, stringsAsFactors=F)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023")
intogen.in <- "other_data/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
cohorts.in <- "other_data/2024-06-18_IntOGen-Cohorts/cohorts.tsv"
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
idat.all <- read.table(intogen.in, sep="\t", header=T, check.names=F)

# Load total sample size for each cancer type
cohorts <- read.table(cohorts.in, header=T, sep="\t", check.names=F, quote="")
cancer.n <- sapply(sort(unique(cohorts$CANCER)), function(cancer){
  sum(cohorts[which(cohorts$CANCER == cancer), "SAMPLES"], na.rm=T)
})

# Collapse driver compendium across studies into unique (gene, cancer) pairs
idat <- as.data.frame(do.call("rbind", lapply(sort(unique(idat.all$SYMBOL)), function(gene){
  gdat <- idat.all[which(idat.all$SYMBOL == gene), ]
  do.call("rbind", lapply(sort(unique(gdat$CANCER_TYPE)), function(cancer){
    cdat <- gdat[which(gdat$CANCER_TYPE == cancer), ]
    driver.total.n <- sum(cdat$TOTAL_SAMPLES, na.rm=T)
    driver.freq <- sum(cdat$SAMPLES, na.rm=T) / driver.total.n
    driver.prop <- driver.total.n / as.numeric(cancer.n[cancer])
    c("SYMBOL" = gene,
      "CANCER_TYPE" = cancer,
      "METHODS" = paste(sort(unique(unlist(strsplit(cdat$METHODS, split=",")))), collapse=","),
      "DRIVER_FREQ" = driver.freq,
      "DRIVER_PROP" = driver.prop,
      "FREQ_LOW_BOUND" = driver.freq * driver.prop,
      "BEST_Q" = min(cdat$QVALUE_COMBINATION, na.rm=T))
  }))
})))

# Only retain genes with FDR Q<0.05 in at least one cohort
idat <- idat[which(idat$BEST_Q < 0.05), ]

# Only retain genes with either:
# A. >=5% frequency in cohorts where driver status was reported, and those cohorts
#    comprise >10% of all samples in that cancer type, or
# B. Lower bound of >=1% frequency in cancer type
a.idx <- which(idat$DRIVER_FREQ >= 0.05 & idat$DRIVER_PROP >= 0.1)
b.idx <- which(idat$FREQ_LOW_BOUND >= 0.01)
idat <- idat[union(a.idx, b.idx), ]

# Only retain genes implicated by intOGen combination or multiple methods
idat <- idat[union(which(idat$METHODS == "combination"),
                   grep(",", idat$METHODS)), ]

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

