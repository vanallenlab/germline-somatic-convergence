#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Find candidate "hub" genes involved in several independent germline:somatic pairs

# Note that this code is meant for local analysis and references hardcoded paths
# on RLC's machine. This will not generalize to all systems.

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)
setwd("~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023/")

# Load all gene pairs
pairs <- read.table("results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.union.tsv",
                    header=T, sep="\t", comment.char="", check.names=F)
colnames(pairs)[1] <- gsub("^#", "", colnames(pairs)[1])

# Exclude gene pairs where coding variation implicates the same gene (not as interesting)
drop.idx <- intersect(which(pairs$germline_context == "coding"
                            & pairs$somatic_context == "coding"),
                      grep("same_gene", pairs$criteria))
pairs <- pairs[-drop.idx, ]

# Count number of times each gene appears in a germline:somatic pair
all.genes <- sort(unique(c(pairs$germline_gene, pairs$somatic_gene)))
gene.counts <- sort(sapply(all.genes, function(g){length(which(pairs$germline_gene == g | pairs$somatic_gene == g))}), decreasing=TRUE)
gene.counts.plot <- rev(gene.counts[which(gene.counts>=10)])
pdf("results/pairs_per_gene.barplot.pdf", height=4.5, width=8)
par(mar=c(0.25, 4, 2.5, 0.25))
h <- barplot(gene.counts.plot, horiz=T, xaxs="i", yaxs="i", xaxt="n",
        names=NA, col="#02679A", border="#003354")
abline(v=0:1000, col="white")
RLCtools::clean.axis(3, title="Convergent pair memberships per gene",
                     label.units="count", cex.axis=5/6, tck=-0.015)
sapply(1:length(gene.counts.plot), function(i){
  axis(2, at=h[i, 1], tick=F, line=-0.9, cex.axis=3/6,
       font=if(gene.counts.plot[i] >= 20){4}else{3},
       labels=names(gene.counts.plot)[i], las=2)
})
dev.off()

# # Print list of pairs involving top genes
# for(gene in names(gene.counts[which(gene.counts >= 10)])){
#   cat(paste("\n\n", gene, ":", sep=""))
#   print(pairs[which(pairs$germline_gene == gene | pairs$somatic_gene == gene), ])
# }
