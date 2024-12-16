#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2024
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Code to filter GWAS catalog loci for a single cancer type

# Set parameters
options(scipen=1000, stringsAsFactors=F)

# Takes exactly two positional arguments:
# 1. Raw GWAS loci as .tsv
# 2. Desired output .tsv for filtered GWAS loci
args <- commandArgs(trailingOnly=TRUE)

# # Dev parameters
# args <- c("/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gwas_catalog/breast.gwas_catalog.12_05_24.unfiltered.tsv",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gwas_catalog/breast.gwas_catalog.12_05_24.filtered.tsv")

# Read raw loci
loci <- read.table(args[1], header=T, sep="\t", check.names=F, quote="")

# Only retain loci with at least one mapped gene
loci <- loci[which(loci$MAPPED_GENE != "" & !is.na(loci$MAPPED_GENE)), ]

# Only retain loci with P <= 5x10-8
loci <- loci[which(loci$PVALUE_MLOG >= -log10(5*10^-8)), ]

# Only retain loci with a specific & defined alternate allele
loci <- loci[which(!grepl("?", loci$`STRONGEST SNP-RISK ALLELE`, fixed=T)), ]

# Drop epistatic associations (marked with an "x" in their coordinates)
loci <- loci[!grepl("x", loci$CHR_POS, fixed=T), ]

# Only keep studies from the last decade
loci <- loci[which(as.numeric(format(as.Date(loci$DATE), "%Y")) > 2014), ]

# Only keep case:control studies (not case:case)
loci <- loci[which(grepl("cases", loci$`INITIAL SAMPLE SIZE`) &
                     grepl("controls", loci$`INITIAL SAMPLE SIZE`)), ]

# Only keep studies with >=10k total cases and >=10k total controls
loci$case_n <- sapply(loci$`INITIAL SAMPLE SIZE`, function(s){
  cohort.comp <- unlist(strsplit(s, split=", "))
  sum(as.numeric(unlist(strsplit(gsub(",", "", cohort.comp[grep("cases", cohort.comp)]),
                                 split="[ ]+"))), na.rm=T)
})
loci$control_n <- sapply(loci$`INITIAL SAMPLE SIZE`, function(s){
  cohort.comp <- unlist(strsplit(s, split=", "))
  sum(as.numeric(unlist(strsplit(gsub(",", "", cohort.comp[grep("controls", cohort.comp)]),
                                 split="[ ]+"))), na.rm=T)
})
loci <- loci[which(loci$case_n >= 10000 & loci$control_n >= 10000), ]

# Write filtered loci to output .tsv
write.table(loci, args[2], sep="\t", quote=F, col.names=T, row.names=F)
