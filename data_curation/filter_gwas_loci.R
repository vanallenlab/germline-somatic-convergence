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
# args <- c("/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gwas_catalog/menarche.gwas_catalog.01_13_25.unfiltered.tsv",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gwas_catalog/menarche.gwas_catalog.01_13_25.filtered.tsv")

# Read raw loci
loci <- read.table(args[1], header=T, sep="\t", check.names=F, quote="")

# Only retain autosomal variants
loci <- loci[which(loci$CHR_ID %in% 1:22), ]

# Only retain loci with at least one mapped gene
loci <- loci[which(loci$MAPPED_GENE != "" & !is.na(loci$MAPPED_GENE)), ]

# Only retain loci with P <= 5x10-8
loci <- loci[which(loci$PVALUE_MLOG >= -log10(5*10^-8)), ]

# Only retain loci with a defined alternate allele (or, at very least, an rsID)
has.alt.idx <- which(!grepl("?", loci$`STRONGEST SNP-RISK ALLELE`, fixed=T))
has.rsid <- grep("^rs[0-9]+", loci$`STRONGEST SNP-RISK ALLELE`)
loci <- loci[union(has.alt.idx, has.rsid), ]

# Drop epistatic associations (marked with an "x" in their coordinates)
loci <- loci[!grepl("x", loci$CHR_POS, fixed=T), ]

# Only keep studies from the last decade
loci <- loci[which(as.numeric(format(as.Date(loci$DATE), "%Y")) > 2014), ]

# Only include case-control or population association studies
case.control.idx <- intersect(grep("cases", loci$`INITIAL SAMPLE SIZE`),
                              grep("controls", loci$`INITIAL SAMPLE SIZE`))
genpop.idx <- grep("individuals", loci$`INITIAL SAMPLE SIZE`)
loci <- loci[union(case.control.idx, genpop.idx), ]

# Only keep studies with >=10k total cases and >=10k total controls
# or, for population studies, at least >=20k total individuals
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
loci$individual_n <- sapply(loci$`INITIAL SAMPLE SIZE`, function(s){
  cohort.comp <- unlist(strsplit(s, split=", "))
  sum(as.numeric(unlist(strsplit(gsub(",", "", cohort.comp[grep("individuals", cohort.comp)]),
                                 split="[ ]+"))), na.rm=T)
})
cc.n.pass <- (loci$case_n >= 10000 & loci$control_n >= 10000)
total.n.pass <- loci$individual_n >= 20000
loci <- loci[which(cc.n.pass | total.n.pass), ]

# Write filtered loci to output .tsv
write.table(loci, args[2], sep="\t", quote=F, col.names=T, row.names=F)
