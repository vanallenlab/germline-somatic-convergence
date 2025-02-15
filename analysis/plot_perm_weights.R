# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Plot distributions of gene sampling probabilities for permutations

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)
require(OncoModR)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023/")

# Read weights
load.weights <- function(path){
  df <- read.table(path, header=T, sep="\t", comment.char="")
  colnames(df) <- c("gene", "weight")
  df$weight <- as.numeric(df$weight)
  df$prob <- df$weight / sum(df$weight)
  return(df)
}
weights <- list("coding_gwas" = load.weights("other_data/permutation_weights/gene_weights.composite_germline_coding.tsv"),
                "noncoding_gwas" = load.weights("other_data/permutation_weights/gene_weights.composite_germline_noncoding.tsv"),
                "coding_mu" = load.weights("other_data/permutation_weights/gene_weights.coding_nonsynonymous.tsv"),
                "somatic_coding" = load.weights("other_data/permutation_weights/gene_weights.coding_nonsynonymous.somatic_ppi_adjustment.tsv"),
                "somatic_noncoding" = load.weights("other_data/permutation_weights/gene_weights.genome_territory.tsv"))

# Plot weights
plot.weights <- function(df, title=NULL, units="Sampling probability", breaks=1000){
  par(mar=c(4, 3.5, 1.5, 0.5))
  xlims <- c(0, quantile(df$prob, 0.99))
  hist(df$prob, col="dodgerblue4", ylab="", xlab="", yaxt="n", xaxt="n", yaxs="i",
       breaks=breaks, xlim=xlims, main=title, xpd=T)
  clean.axis(1, title=units, infinite=T, labels=paste(round(100*axTicks(1), 3), "%", sep=""))
  clean.axis(2, title="Genes", title.line=1.5, infinite=T)
}
pdf("results/permutation_weight_distributions.pdf",
    height=3, width=10)
par(mfrow=c(1, 5))
plot.weights(weights[[1]], title="GWAS (coding)", breaks=200)
plot.weights(weights[[2]], title="GWAS (noncoding)", breaks=500)
plot.weights(weights[[3]], title="Coding rare variants", breaks=500)
plot.weights(weights[[4]], title="Somatic coding", breaks=1000)
plot.weights(weights[[5]], title="Somatic noncoding", breaks=500)
dev.off()
