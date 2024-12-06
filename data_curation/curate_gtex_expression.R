#!/usr/bin/env Rscript

# Germline:Somatic Convergence Analysis 2023
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Code to curate GTEx gene expression matrixes for permutation weights

# Set parameters
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=T)
PedSV::load.constants("scales")
require(viridis, quietly=T)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023")
gtex.in <- "other_data/gtex_median_expression_gencode_filtered_SH_20231205.csv"
outdir <- "other_data/permutation_weights"


# Read GTEx data
gtex <- read.table(gtex.in, sep=",", header=T, check.names=F)


# Get organ-level averages for organs with multiple biopsy sites
gtex$Colon <- apply(gtex[, grep("Colon", colnames(gtex))], 1, mean, na.rm=T)
gtex <- gtex[, -grep("Colon - ", colnames(gtex), fixed=T)]
gtex$Kidney <- apply(gtex[, grep("Kidney", colnames(gtex))], 1, mean, na.rm=T)
gtex <- gtex[, -grep("Kidney - ", colnames(gtex), fixed=T)]


# Deduplicate genes
gtex <- gtex[-which(duplicated(gtex$Gene)), ]


# Divide all genes for each tissue into quintiles by expression in tissue
# Genes with strictly zero expression in a tissue are assigned to a sixth "quintile"
gtex.raw <- gtex
# Add a fifth bin for all unexpressed genes
gtex[, -1] <- apply(gtex[, -1], 2, function(vals){
  # First, hold out all genes with strictly zero expression
  # These will be assigned to their own "quintile"
  zero.idxs <- which(vals == 0)
  assignments <- rep(NA, length(vals))
  assignments[zero.idxs] <- 1

  # Second, divide non-zero genes into quintiles per tissue
  assignments[-zero.idxs] <- ceiling(5 * rank(vals[-zero.idxs]) / length(vals[-zero.idxs])) + 1
  return(assignments)
})

# Make matrix of weights for all tissues
t2c.map <- c("Breast - Mammary Tissue" = "breast",
             "Lung" = "lung",
             "Prostate" = "prostate",
             "Colon" = "colorectal",
             "Kidney" = "renal")
expression.weights <- do.call("rbind", lapply(colnames(gtex)[-1], function(tissue){
  sub.df <- gtex[, c(1, which(colnames(gtex) == tissue))]
  sub.df$cancer <- t2c.map[tissue]
  colnames(sub.df) <- c("gene", "weight", "cancer")
  return(sub.df)
}))

# Write expression weights to file
write.table(expression.weights,
            paste(outdir, "gene_weights.expression.cancer_specific.tsv", sep="/"),
            col.names=T, row.names=F, sep="\t", quote=F)

# Generate a histogram of expression values for each tissue colored by quintile
pdf(paste(outdir, "/../../results/expression_bins_per_tissue.pdf", sep=""),
    height=4.5, width=10)
par(mfrow=c(2, 3), bty="n")
sapply(1:length(t2c.map), function(k){
  vals <- log10(gtex.raw[, names(t2c.map)[k]] + 10e-3)
  h <- hist(vals, breaks=100, plot=F)
  prep.plot.area(xlims=range(h$breaks), ylims=c(0, max(h$counts)),
                 parmar=c(3, 3, 2, 1))
  qb <- c(-Inf, sapply(2:6, function(q){min(vals[which(gtex[, names(t2c.map)[k]] == q)])}), Inf)
  sapply(1:6, function(q){
    lower <- qb[q]; upper <- qb[q+1]
    h.idxs <- which(h$mids > lower & h$mids <= upper)
    rect(xleft=h$breaks[h.idxs], xright=h$breaks[h.idxs+1],
         ybottom=0, ytop=h$counts[h.idxs],
         col=viridis(6)[q])
  })
  clean.axis(1, at=log10(logscale.major), labels=logscale.major, title="GTEx TPM")
  clean.axis(2, title="Genes")
  mtext(3, line=0, font=2, 
        text=paste(toupper(substr(t2c.map[k], 1, 1)), 
                   substr(t2c.map[k], 2, 100), sep=""))
})
prep.plot.area(c(0, 1), c(0, 1), rep(1, 4))
legend("center", legend=rev(c("Not expressed", paste("Q", 1:5, sep=""))),
       fill=rev(viridis(6)), bty="n", cex=1.3)
dev.off()
