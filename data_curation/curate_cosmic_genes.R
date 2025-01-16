#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Code to curate COSMIC gene lists of germline and somatic drivers

# Set parameters
options(scipen=1000, stringsAsFactors=F)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023/")
cosmic.in <- "other_data/COSMIC.CGC.12_05_24.tsv"
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0){
  outdir <- as.character(args[1])
}else{
  outdir <- "gene_lists"
}


# Read COSMIC gene census
cdat <- read.table(cosmic.in, sep="\t", header=T)

# Build map of eligible cancer terms
# This must be done manually from COSMIC metadata as there is quite a variety
all.ctypes <- unlist(strsplit(unlist(as.vector(cdat[, c("Tumour.Types.Somatic.",
                                                        "Tumour.Types.Germline.")])),
                              split=","))
all.ctypes <- sort(unique(sub("^[ ]+|[ ]+$", "", all.ctypes)))
write.table(all.ctypes, "other_data/all_COSMIC_cancer_types.list", quote=F,
            row.names=F, col.names=F)
cancer.map <- list("Breast" = c("breast", "breast cancer", "breast carcinoma",
                                "lobular breast", "luminal A breast",
                                "phyllodes tumour", "phyllodes tumour of the breast",
                                "secretory breast"),
                   "Colorectal" = c("colon", "colon adenocarcinoma", "colon cancer",
                                    "colon carcinoma", "colorectal",
                                    "colorectal adenocarcinoma", "colorectal cancer",
                                    "colorectal cancer susceptibility", "CRC",
                                    "large intestine", "large intestine carcinoma",
                                    "rectal cancer"),
                   "Lung" = c("lung", "lung adenocarcinoma", "lung cancer",
                              "lung carcinoma", "lung SCC", "Lung SCC", "NSCLC",
                              "SCLC", "small cell lung carcinoma"),
                   "Melanoma" = c("cutaneous melanoma", "desmoplastic melanoma",
                                  "malignant melanoma of soft parts", "melanoma",
                                  "mucosal melanoma", "skin and uveal melanoma",
                                  "uveal melanoma", "malignant melanoma"),
                   "Prostate" = c("prostae adenocarcinoma", "prostate",
                                  "prostate cancer", "prostate carcinoma"),
                   "Renal" = c("CCRCC", "clear cell renal carcinoma",
                               "clear cell renal cell carcinoma", "kidney cancer",
                               "papillary renal", "RCC", "renal",
                               "renal angiomyolipoma", "renal cell carcinoma",
                               "renal cell carcinoma (childhood epithelioid)"))

# Annotate COSMIC data with simple true/false columns for each cancer type
for(cancer in names(cancer.map)){
  # Somatic annotation
  somatic.labs <- sapply(cdat$Tumour.Types.Somatic., function(cstr){
    cterms <- sub("^[ ]+|[ ]+$", "", unlist(strsplit(cstr, split=",")))
    length(intersect(cancer.map[[cancer]], cterms)) > 0
  })
  cdat[, paste(cancer, "somatic", sep=".")] <- somatic.labs

  # Germline annotation
  germline.labs <- sapply(cdat$Tumour.Types.Germline., function(cstr){
    cterms <- sub("^[ ]+|[ ]+$", "", unlist(strsplit(cstr, split=",")))
    length(intersect(cancer.map[[cancer]], cterms)) > 0
  })
  cdat[, paste(cancer, "germline", sep=".")] <- germline.labs
}

# Write lists of germline & somatic drivers to gene list files
for(cancer in names(cancer.map)){
  for(context in c("germline", "somatic")){
    outdir.sub <- paste(outdir, "/", context, "_coding", sep="")
    if(!dir.exists(outdir.sub)){
      dir.create(outdir.sub)
    }
    genes <- sort(unique(cdat[which(cdat[, paste(cancer, context, sep=".")]), 1]))
    f.out <- paste(outdir.sub, "/", tolower(cancer),
                   ".", context, ".coding.genes.list", sep="")
    write.table(genes, f.out, col.names=F, row.names=F, quote=F)
  }
}

