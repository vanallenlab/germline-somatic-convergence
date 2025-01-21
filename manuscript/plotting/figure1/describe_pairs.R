#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Descriptive analysis of convergent gene pairs and hub genes

# Set parameters and load libraries
options(scipen=1000, stringsAsFactors=FALSE)
require(argparse)
require(RLCtools)

# Parse command line arguments and options
parser <- ArgumentParser(description="Summarize convergent pairs")
parser$add_argument("--config", metavar=".R", type="character", required=TRUE,
                    help="path to config.R file")
parser$add_argument("--project-root", metavar="path", type="character",
                    help="root directory for all project data [default: pwd]")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    default="./pair_analysis", help="file prefix for all plots")
args <- parser$parse_args()

# # DEV arguments
# args <- list("config" = "~/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence/manuscript/plotting/config.R",
#              "project_root" = "~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023",
#              "out_prefix" = "~/scratch/gsc_dev")

# Check --project-root
if(is.null(args$project_root)){
  args$project_root <- getwd()
}else{
  args$project_root <- gsub("/$", "", args$project_root)
}

# Source config.R
source(args$config)

# Read pairs
pairs <- read.table(paste(args$project_root, "results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.union.tsv", sep="/"),
                    header=T, sep="\t", comment.char="", check.names=F)
colnames(pairs)[1] <- gsub("^#", "", colnames(pairs)[1])

# Collect data for main pair count barplot
bar.dat <- as.data.frame(do.call("rbind", lapply(contexts, function(germ.context){
  sub.df <- as.data.frame(do.call("rbind", lapply(contexts, function(som.context){
    sub.df <- as.data.frame(do.call("rbind", lapply(cancers, function(cancer){
      v <- sapply(1:4, function(tier){
        length(which(pairs$cancer == cancer
                     & pairs$germline_context == germ.context
                     & pairs$somatic_context == som.context
                     & pairs$tier == tier))
      })
      c(cancer, v)
    })))
    as.data.frame(cbind(rep(som.context, nrow(sub.df)), sub.df))
  })))
  as.data.frame(cbind(rep(germ.context, nrow(sub.df)), sub.df))
})))
colnames(bar.dat) <- c("germline", "somatic", "cancer", names(tier.pch))
bar.dat[, names(tier.pch)] <- apply(bar.dat[, names(tier.pch)], 2, as.numeric)
bar.dat$total <- apply(bar.dat[, grep("^tier", colnames(bar.dat))], 1, sum)

# Barplot of convergent pairs; grouped by context; colored by cancer; shaded by tier
group.spacer <- 1
bar.xbuf <- 0.15
pdf(paste(args$out_prefix, "pair_count_bars.pdf", sep="."),
    height=2.3, width=4.75)
prep.plot.area(xlims=c(-1, (4*n.cancers)+(3*group.spacer)),
               ylims=c(0, 1.1*max(bar.dat$total)),
               parmar=c(2, 3.75, 0.1, 0))
sapply(0:3, function(i){
  ax.at <- c(0, n.cancers)+(i*n.cancers)+(i*group.spacer)
  axis(1, at=ax.at, tck=0, labels=NA)
  axis(1, at=mean(ax.at), tick=F, cex.axis=5/6, line=-0.9,
       labels=if(i<2){"Coding"}else{"Noncoding"},
       col.axis=if(i<2){"black"}else{"gray70"})
  axis(1, at=mean(ax.at), tick=F, cex.axis=5/6, line=0,
       labels=if(i %in% c(0, 2)){"Coding"}else{"Noncoding"},
       col.axis=if(i %in% c(0, 2)){"black"}else{"gray70"})
})
axis(1, at=-1, tick=F, xpd=T, cex.axis=5.5/6, line=-0.9, label="Germline:", hadj=1)
axis(1, at=-1, tick=F, xpd=T, cex.axis=5.5/6, line=-0, label="Somatic:", hadj=1)
clean.axis(2, label.units="count", title="Convergent\ngene pairs", infinite.positive=TRUE)
sapply(1:n.cancers, function(ki){
  cidxs <- ki+((0:3)*n.cancers)+((0:3)*group.spacer)
  y.at <- t(apply(bar.dat[which(bar.dat$cancer == cancers[ki]),
                          paste("tier", 1:4, sep="")], 1, cumsum))
  rect.cols <- sapply(2:5, function(ci){rep(cancer.palettes[[cancers[ki]]][[ci]], 4)})
  rect(xleft=rep(cidxs-1+bar.xbuf, 4), xright=rep(cidxs-bar.xbuf, 4),
       ybottom=as.vector(cbind(rep(0, 4), y.at[, 1:3])),
       ytop=as.vector(y.at), border=rect.cols, col=rect.cols, lwd=0.5)
  rect(xleft=cidxs-1+bar.xbuf, xright=cidxs-bar.xbuf,
       ybottom=rep(0, 4), ytop=y.at[, 4])
  text(x=cidxs-0.5, y=y.at[, 4]+(0.04*diff(par("usr")[3:4])),
       labels=prettyNum(y.at[, 4], big.mark=","),
       col=unlist(cancer.colors[[cancers[ki]]]), cex=4/6, xpd=T)
})
dev.off()

# Prep hub gene plots
hub.cutoff <- 4
somatic.hubs <- names(which(sort(table(pairs$somatic_gene), decreasing=TRUE) >= hub.cutoff))
germline.hubs <- names(which(sort(table(pairs$germline_gene), decreasing=TRUE) >= hub.cutoff))
som.hub.prop <- length(somatic.hubs) / length(c(somatic.hubs, germline.hubs))
total.hub.width <- 7
hub.height <- 2.65
plot.hubs <- function(hub.genes, origin, pt.cex=0.775){
  n.genes <- length(hub.genes)
  max.hits <- max(unlist(apply(pairs[, c("germline_gene", "somatic_gene")], 2, table)))
  prep.plot.area(xlims=c(0, n.genes-0.25), ylims=c(-0.5, max.hits),
                 parmar=c(2.4, 3.5, 0.1, 0.25))
  clean.axis(2, label.units="count", title="Convergent pairs\ncontaining gene",
             infinite.positive=TRUE)
  sapply(1:n.genes, function(i){
    gene <- hub.genes[i]
    hub.df <- pairs[which(pairs[, paste(origin, "gene", sep="_")] == gene), ]
    n.points <- nrow(hub.df)
    v.cancers <- hub.df$cancer
    v.tier <- paste("tier", hub.df$tier, sep="")
    v.germ <- hub.df$germline_context
    v.som <- hub.df$somatic_context
    sapply(n.points:1, function(k){
      points(x=i-0.5, y=k-0.5, pch=tier.pch[[v.tier[k]]], xpd=T, cex=pt.cex,
             bg=cancer.palettes[[v.cancers[k]]][[context.shades[[v.germ[k]]]]],
             col=cancer.palettes[[v.cancers[k]]][[context.shades[[v.som[k]]]]],
             lwd=1.5*pt.cex)
    })
    axis(1, at=i-c(0.25, 0.75), tck=0, labels=NA, xpd=T)
    text(x=i-0.45, y=-0.05*diff(par("usr")[3:4]), srt=45, adj=1,
         labels=hub.genes[i], font=3, cex=4/6, xpd=T)
  })
}

# Somatic hub genes
pdf(paste(args$out_prefix, "hub_genes.somatic.pdf", sep="."),
    height=hub.height, width=total.hub.width*som.hub.prop)
plot.hubs(somatic.hubs, "somatic")
dev.off()

# Germline hub genes
pdf(paste(args$out_prefix, "hub_genes.germline.pdf", sep="."),
    height=hub.height, width=total.hub.width*som.hub.prop)
plot.hubs(germline.hubs, "germline")
dev.off()
