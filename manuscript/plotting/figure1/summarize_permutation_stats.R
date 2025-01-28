#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Summarize raw results from permutation testing

# Set parameters and load libraries
options(scipen=1000, stringsAsFactors=FALSE)
require(argparse)
require(RLCtools)

# Declare constants
models <- c("uniform", "bayesian", "expression", "composite")
n.effective.tests <- 120

# Helper function to get pointwise graphical parameters
get.pw.params <- function(stats){
  pwpch <- unlist(tier.pch[stats$criteria])
  pwbg <- unlist(sapply(1:nrow(stats), function(i){
    cancer <- stats$cancer[[i]]
    if(cancer %in% names(cancer.palettes)){
      cancer.palettes[[cancer]][[context.shades[[stats$germline[i]]]]]
    }else{
      hex2grey(cancer.palettes[[unlist(strsplit(cancer, split="_"))[1]]][[context.shades[[stats$germline[i]]]]])
    }
  }))
  pwcol <- unlist(sapply(1:nrow(stats), function(i){
    cancer <- stats$cancer[[i]]
    if(cancer %in% names(cancer.palettes)){
      cancer.palettes[[stats$cancer[i]]][[context.shades[[stats$somatic[i]]]]]
    }else{
      hex2grey(cancer.palettes[[unlist(strsplit(cancer, split="_"))[1]]][[context.shades[[stats$somatic[i]]]]])
    }
  }))
  pwsym <- c("TRUE" = "+", "FALSE" = "")[as.character(grepl("plus|any", stats$criteria))]
  # linecol <- sapply(1:nrow(stats), function(i){
  #   hex2grey(cancer.palettes[[stats$cancer[i]]][["light2"]])
  # })
  linecol <- rep("gray90", length(pwpch))
  return(list("pwpch" = pwpch, "pwbg" = pwbg, "pwcol" = pwcol,
              "pwsym" = pwsym, linecol=linecol))
}

# Parse command line arguments and options
parser <- ArgumentParser(description="Summarize permutation statistics")
parser$add_argument("--config", metavar=".R", type="character", required=TRUE,
                    help="path to config.R file")
parser$add_argument("--model", metavar="string", type="character", required=TRUE,
                    help=paste("specify which model results should be analyzed ",
                               "(options: ", paste(models, collapse="|"), ")", sep=""))
parser$add_argument("--project-root", metavar="path", type="character",
                    help="root directory for all project data [default: pwd]")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    default="./perm_stats", help="file prefix for all plots")
args <- parser$parse_args()

# # DEV arguments
# args <- list("config" = "~/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence/manuscript/plotting/config.R",
#              "model" = "composite",
#              "project_root" = "~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023",
#              "out_prefix" = "~/scratch/gsc_dev")

# Check --project-root
if(is.null(args$project_root)){
  args$project_root <- getwd()
}else{
  args$project_root <- gsub("/$", "", args$project_root)
}

# Check --model
if(!args$model %in% models){
  stop(paste("--model", args$model, "not a recognized option. See --help for more info."))
}

# Source config.R
source(args$config)

# Read raw permutation stats
stats <- read.table(paste(args$project_root, "/results",
                          "/VALab_germline_somatic_2024.v2.all_permutation_results",
                          "/cosmic_intogen_union_results",
                          "/cosmic_intogen_union.", args$model, "_sampling_results",
                          "/cosmic_intogen_union.", args$model,
                          "_sampling.all_stats.tsv", sep=""),
                    header=T, sep="\t", comment.char="", check.names=F)

# Sort stats to ensure negative control always appears directly after corresponding cancer
stats <- stats[with(stats, order(cancer, germline, somatic, criteria)), ]

# Add helper flag for whether the row is a negative control point
stats$nc <- grepl("_", stats$cancer)

# Convert fold-enrichments to log2 and winsorize lower bound
stats$lfc <- log2(stats$fold)
min.lfc <- -2
stats$lfc[which(stats$lfc < min.lfc)] <- min.lfc


# Plot layered QQ of negative control and true results
# Only show results for strata with at least 1 expected pair
qq.df <- stats[grep("plus", stats$criteria, invert=TRUE), ]
qq.df <- qq.df[which(qq.df$expected >= 1), ]
qq.df <- qq.df[grep("^any", qq.df$cancer, invert=T), ]
nc.qq.df <- qq.df[which(qq.df$nc), ]
pos.qq.df <- qq.df[which(!qq.df$nc), ]
qq.pt.cex <- 0.75
qq.ymax <- 15
pdf(paste(args$out_prefix, args$model, "qq.pdf", sep="."), height=2.2, width=2.2)
plot.qq(pos.qq.df$p, ymax=qq.ymax, cap.pvals=TRUE,
        pt.cex=qq.pt.cex, do.fdr=FALSE, cutoff=0.05/n.effective.tests, label.cex=5/6,
        pt.color=unlist(cancer.colors[pos.qq.df$cancer]),
        ci.color="gray95", oe.line.color="gray60",
        y.title=expression("Obs." ~ -log[10] ~ italic(P)),
        x.label.line=-1, x.title.line=0.1, y.title.line=0.15,
        parmar=c(2, 2.25, 0.25, 0.25))
plot.qq(nc.qq.df$p, pt.color="grey70", add=T, do.fdr=FALSE, pt.cex=qq.pt.cex)
dev.off()


# Volcano-style plot of P-values and fold-enrichments
volcano.logp.max <- 1.05*15
volcano.lfc.max <- log2(64)
volcano.df <- stats[grep("plus", stats$criteria, invert=TRUE), ]
volcano.df <- volcano.df[grep("^any", volcano.df$cancer, invert=T), ]
volcano.df$p[which(volcano.df$p < 10^-volcano.logp.max)] <- 10^-volcano.logp.max
volcano.df$lfc[which(volcano.df$lfc > volcano.lfc.max)] <- volcano.lfc.max
pdf(paste(args$out_prefix, args$model, "volcano.pdf", sep="."),
    width=3.3, height=2.2)
prep.plot.area(xlims=c(min.lfc, volcano.lfc.max), ylims=c(0, volcano.logp.max),
               xaxs="r", yaxs="r", parmar=c(2, 2, 0.25, 0.25))
abline(h=-log10(c(0.05/n.effective.tests, 0.05)), lty=c(5, 2),
       col=c("black", "gray50"))
pw.params <- get.pw.params(volcano.df)
sapply(1:nrow(volcano.df), function(k){
  points(x=volcano.df$lfc[k], y=-log10(volcano.df$p)[k],
         cex=0.9, pch=pw.params$pwpch[k],
         col=pw.params$pwcol[k], bg=pw.params$pwbg[k])
  points(x=volcano.df$lfc[k], y=-log10(volcano.df$p)[k],
         cex=0.7, pch=pw.params$pwsym[k], col="white")

})
x.ax.at <- axTicks(1)
x.ax.labs <- 2^x.ax.at
x.ax.labs[1] <- paste("\"\"<", x.ax.labs[1], sep="")
x.ax.labs[length(x.ax.labs)] <- paste("\"\">", x.ax.labs[length(x.ax.labs)], sep="")
clean.axis(1, at=x.ax.at, labels=x.ax.labs, parse.labels=TRUE,
           title=expression(log[2] ~ "fold-enrichment vs. permutation"),
           infinite=TRUE, label.line=-0.9, title.line=0.1)
y.ax.labs <- y.ax.at <- axTicks(2)
y.ax.labs[length(y.ax.labs)] <- paste("\"\">", y.ax.labs[length(y.ax.labs)])
clean.axis(2, at=y.ax.at, labels=y.ax.labs, title=expression(-log[10] ~ italic(P)),
           infinite=TRUE, parse.labels=TRUE, title.line=0)
dev.off()

