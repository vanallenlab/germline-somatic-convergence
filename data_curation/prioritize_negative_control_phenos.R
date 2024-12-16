#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Prioritize a list of candidate negative control phenotypes for permutation testing

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)

# # DEV:
# args <- c("/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/data/negative_control_phenos/neale_ukbb.candidate_rg_data.csv",
#           "/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/data/negative_control_phenos/ukb31063_h2_topline.02Oct2019.tsv.gz")

# Read command-line arguments
args <- commandArgs(trailingOnly=TRUE)
rg.csv <- as.character(args[1])
h2g.tsv <- as.character(args[2])


# Load genetic correlations
rg <- read.table(rg.csv, sep=",", header=T)

# Filter genetic correlations
rg <- rg[which(rg$p > 0.25 & abs(rg$rg) < 0.02 & abs(rg$rpheno) < 0.02), ]

# Load heritability estimates
h2g <- read.table(h2g.tsv, sep="\t", comment.char="", header=T, quote="")

# Filter heritability data
h2g <- h2g[which(h2g$h2_liability > 0.1 & h2g$h2_p < 0.05 / nrow(h2g)
                 & h2g$n_cases > 5000 & h2g$n_controls > 20000), ]

# Find the intersection of qualified rg + h2g data
candidates <- intersect(h2g$phenotype, rg$ID2)
cat("Candidate phenotypes:\n")
cat(paste(h2g$description[which(h2g$phenotype %in% candidates)], collapse="\n"))

