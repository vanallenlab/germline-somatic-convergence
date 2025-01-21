#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Summarize raw results from permutation testing

# Set parameters and load libraries
options(scipen=1000, stringsAsFactors=FALSE)
require(argparse)
require(beeswarm)
require(RLCtools)

# Parse command line arguments and options
parser <- ArgumentParser(description="Summarize permutation statistics")
parser$add_argument("--config", metavar=".R", type="character", required=TRUE,
                    help="path to config.R file")
parser$add_argument("--project-root", metavar="path", type="character",
                    help="root directory for all project data [default: pwd]")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    default="./perm_stats", help="file prefix for all plots")
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

# Read raw permutation stats
stats <- read.table(paste(args$project_root, "results",
                          "VALab_germline_somatic_2024.v2.all_permutation_results",
                          "cosmic_intogen_union_results",
                          "cosmic_intogen_union.composite_sampling_results",
                          "cosmic_intogen_union.composite_sampling.all_stats.tsv", sep="/"),
                    header=T, sep="\t", comment.char="", check.names=F)

# Sort stats to ensure negative control always appears directly after corresponding cancer
stats <- stats[with(stats, order(cancer, germline, somatic, criteria)), ]

# Add helper flag for whether the row is a negative control point
stats$nc <- grepl("_", stats$cancer)

# Convert fold-enrichments to log2 and winsorize lower bound
stats$lfc <- log2(stats$fold)
min.lfc <- -2
stats$lfc[which(stats$lfc < min.lfc)] <- min.lfc

# Function to plot paired & connected swarms of log2(fold) between empirical and negative control points
paired.swarm <- function(stats, column="lfc", at=c(0.5, 1.5), pt.cex=0.5){
  # Because stats has been sorted above, we can assume that melting the dataframe is safe
  emp.df <- stats[which(!stats$nc), ]
  nc.df <- stats[which(stats$nc), ]

  # Get pointwise graphical parameters
  pwpch <- unlist(tier.pch[emp.df$criteria])
  pwbg <- sapply(1:nrow(emp.df), function(i){
    cancer.palettes[[emp.df$cancer[i]]][[context.shades[[emp.df$germline[i]]]]]
  })
  pwcol <- sapply(1:nrow(emp.df), function(i){
    cancer.palettes[[emp.df$cancer[i]]][[context.shades[[emp.df$somatic[i]]]]]
  })
  pwsym <- c("TRUE" = "+", "FALSE" = "")[as.character(grepl("plus|any", emp.df$criteria))]
  # linecol <- sapply(1:nrow(emp.df), function(i){
  #   hex2grey(cancer.palettes[[emp.df$cancer[i]]][["light2"]])
  # })
  linecol <- rep("gray90", length(pwpch))

  # Get coordinates of each swarm
  emp.xy <- beeswarm(emp.df[, column], cex=pt.cex, at=at[1], pwpch=pwpch, add=T,
                     col=NA, bg=NA, corral="wrap", corralWidth=0.75)[, c("x", "y")]
  nc.xy <- beeswarm(nc.df[, column], cex=pt.cex, at=at[2], pwpch=pwpch, add=T,
                    col=NA, bg=NA, corral="wrap", corralWidth=0.75)[, c("x", "y")]

  # Add connecting lines
  segments(x0=emp.xy$x, x1=nc.xy$x, y0=emp.xy$y, y1=nc.xy$y, col=linecol)

  # Plot points
  points(emp.xy, pch=pwpch, cex=pt.cex, bg=pwbg, col=pwcol, xpd=T)
  points(emp.xy, pch=pwsym, cex=0.7*pt.cex, col="white", xpd=T)
  points(nc.xy, pch=pwpch, cex=pt.cex, bg=hex2grey(pwbg), col=hex2grey(pwcol), xpd=T)
  points(nc.xy, pch=pwsym, cex=0.7*pt.cex, col="white", xpd=T)
}

# LFC by tier, all cancers
pdf(paste(args$out_prefix, "lfc_by_tier.emp_vs_negctrl.pdf", sep="."), height=2, width=5)
max.lfc <- ceiling(max(stats$lfc))
prep.plot.area(xlims=c(0, 11), ylims=c(min.lfc, max.lfc),
               parmar=c(1, 3.5, 0.5, 0))
abline(h=0)
paired.swarm(stats[which(stats$criteria == "tier1"), ], at=(1:2)-0.5)
paired.swarm(stats[which(stats$criteria == "tier2"), ], at=(4:5)-0.5)
paired.swarm(stats[which(stats$criteria == "tier3"), ], at=(7:8)-0.5)
paired.swarm(stats[which(stats$criteria == "tier4"), ], at=(10:11)-0.5)
axis(1, at=c(1, 4, 7, 10), tick=F, line=-0.9, labels=paste("Tier", 1:4), cex.axis=5/6)
y.ax.at <- min.lfc:max.lfc
clean.axis(2, at=y.ax.at, labels=2^y.ax.at, infinite.positive=TRUE,
           title="Fold-enrichment vs.\nrandom permutations")
dev.off()

