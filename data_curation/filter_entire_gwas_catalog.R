#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2024
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Filter entire GWAS catalog as preprocessing for permutation weighting

# Set parameters
options(scipen=1000, stringsAsFactors=F)

# Takes exactly two positional arguments:
# 1. Raw GWAS catalog results as .tsv
# 2. Desired output .tsv for filtered GWAS loci
args <- commandArgs(trailingOnly=TRUE)

# # Dev parameters
# args <- c("/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.tsv.gz",
#           "/Users/ryan/Partners HealthCare Dropbox/Ryan Collins/VanAllen/VALab_germline_somatic_2023/other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.filtered.tsv")

# Read raw loci
loci <- read.table(args[1], header=T, sep="\t", check.names=F, quote="")

# Only retain loci with at least one mapped gene
loci <- loci[which(loci$MAPPED_GENE != "" & !is.na(loci$MAPPED_GENE)), ]

# Only retain loci with P <= 5x10-8
loci <- loci[which(loci$PVALUE_MLOG >= -log10(5*10^-8)), ]

# Drop epistatic associations (marked with an "x" in their coordinates)
loci <- loci[!grepl("x", loci$CHR_POS, fixed=T), ]

# Only keep studies from the last decade
loci <- loci[which(as.numeric(format(as.Date(loci$DATE), "%Y")) > 2014), ]

# Exclude any cancer phenotype
drop.terms <- c("cancer", "carcinoma", "neoplasm", "tumor", "tumour",
                "malignancy", "lymphoma", "sarcoma", "leukemia")
keep.idxs <- intersect(grep(paste(drop.terms, collapse="|"),
                            tolower(loci$`DISEASE/TRAIT`), invert=TRUE),
                       grep(paste(drop.terms, collapse="|"),
                            tolower(loci$MAPPED_TRAIT), invert=TRUE))
loci <- loci[keep.idxs, ]

# Deduplicate the same variant mapped to the same gene for the same trait
unique.rows <- which(!duplicated(loci[, c("CHR_ID", "CHR_POS", "MAPPED_GENE", "DISEASE/TRAIT")]))
loci <- loci[unique.rows, ]

# Write filtered loci to output .tsv
write.table(loci, args[2], sep="\t", quote=F, col.names=T, row.names=F)
