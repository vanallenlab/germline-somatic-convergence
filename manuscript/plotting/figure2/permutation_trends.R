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

# Function to plot paired & connected swarms of log2(fold) between empirical and negative control points
paired.swarm <- function(stats, column="lfc", at=c(0.5, 1.5), pt.cex=0.65,
                         top.keys=c("cancer", "germline", "somatic", "criteria")){
  # Because stats has been sorted above, we can assume that melting the dataframe is safe
  emp.df <- stats[which(!stats$nc), ]
  nc.df <- stats[which(stats$nc), ]

  # Get pointwise graphical parameters
  pw.params <- get.pw.params(emp.df)
  pwpch <-  pw.params$pwpch
  pwbg <- pw.params$pwbg
  pwcol <- pw.params$pwcol
  pwsym <- pw.params$pwsym
  linecol <- pw.params$linecol

  # Find point index corresponding to top-most summary; this point will be overplotted
  top.idx <- which(apply(emp.df[, top.keys], 1, function(keys){all(keys == "any")}))

  # Get coordinates of each swarm
  emp.xy <- beeswarm(emp.df[-top.idx, column], cex=pt.cex, at=at[1],
                     pwpch=pwpch[-top.idx], add=T, col=NA, bg=NA,
                     corral="wrap", corralWidth=0.75)[, c("x", "y")]
  nc.xy <- beeswarm(nc.df[-top.idx, column], cex=pt.cex, at=at[2],
                    pwpch=pwpch[-top.idx], add=T, col=NA, bg=NA,
                    corral="wrap", corralWidth=0.75)[, c("x", "y")]

  # Add connecting lines
  segments(x0=emp.xy$x, x1=nc.xy$x, y0=emp.xy$y, y1=nc.xy$y, col=linecol)

  # Plot points
  points(emp.xy, pch=pwpch[-top.idx], cex=pt.cex, bg=pwbg[-top.idx],
         col=pwcol[-top.idx], xpd=T)
  points(emp.xy, pch=pwsym[-top.idx], cex=0.7*pt.cex, col="white", xpd=T)
  points(nc.xy, pch=pwpch[-top.idx], cex=pt.cex, bg=hex2grey(pwbg)[-top.idx],
         col=hex2grey(pwcol)[-top.idx], xpd=T)
  points(nc.xy, pch=pwsym[-top.idx], cex=0.7*pt.cex, col="white", xpd=T)

  # Add overplot point & line for highest-level summary
  segments(x0=at[1], x1=at[2], y0=emp.df[top.idx, column],
           y1=nc.df[top.idx, column], col=pwcol[top.idx], lwd=1.5)
  points(x=at[1], y=emp.df[top.idx, column], pch=pwpch[top.idx], cex=2*pt.cex,
         bg=pwbg[top.idx], col=pwcol[top.idx], xpd=T)
  points(x=at[1], y=emp.df[top.idx, column], pch=pwsym[top.idx],
         cex=0.7*2*pt.cex, col="white", xpd=T)
  points(x=at[2], y=nc.df[top.idx, column], pch=pwpch[top.idx], cex=2*pt.cex,
         bg=hex2grey(pwbg[top.idx]), col=hex2grey(pwcol[top.idx]), xpd=T)
  points(x=at[2], y=nc.df[top.idx, column], pch=pwsym[top.idx],
         cex=0.7*2*pt.cex, col="white", xpd=T)
}

# Function to prepare plot area for paired swarms
prep.paired.swarm <- function(xlims=c(0, 11), ylims=c(min.lfc, max.lfc),
                              parmar=c(3, 3.75, 0.5, 0)){
  prep.plot.area(xlims=xlims, ylims=ylims, parmar=parmar)
  abline(h=0, lty=5)
  sapply(c(1, 4, 7, 10), function(x){
    axis(1, at=x-0.5, tick=F, line=-0.9, cex.axis=5/6, font=3, labels="Obs.")
    axis(1, at=x+0.5, tick=F, line=-0.9, cex.axis=5/6, font=3, labels="Neg.",
         col.axis="gray50")
  })
  y.ax.at <- min.lfc:max.lfc
  y.ax.labels <- 2^y.ax.at
  y.ax.labels[1] <- paste("\"\"<", y.ax.labels[1])
  y.ax.labels[length(y.ax.labels)] <- paste("\"\">", y.ax.labels[length(y.ax.labels)])
  clean.axis(2, at=y.ax.at, labels=y.ax.labels, infinite.positive=TRUE, parse.labels=TRUE,
             title.line=0.75, title="Fold-enrichment vs.\nrandom permutations")
}

# Function to plot connected swarms of log2(fold) between 2+ sets of points
swarm.sequence <- function(stats.list, column="lfc", pt.cex=0.65,
                         top.keys=c("cancer", "germline", "somatic", "criteria")){
  # Confirm the same number of rows appears in all stats.list
  if(length(unique(sapply(stats.list, nrow))) > 1){
    stop("Not all elements of stats.list have the same number of rows")
  }
  n.swarms <- length(stats.list)

  # Virtually swarm each column and save (x, y) coordinates
  swarm.xy <- lapply(1:length(stats.list), function(i){
    this.df <- stats.list[[i]]
    beeswarm(this.df[, column], cex=pt.cex, at=i-0.5, add=T, col=NA, bg=NA,
             corral="wrap", corralWidth=0.75)[, c("x", "y")]
  })

  # Plot one set of lines & points per element in stats.list
  for(i in 1:n.swarms){
    this.df <- stats.list[[i]]
    pw.params <- get.pw.params(this.df)
    this.xy <- swarm.xy[[i]]
    top.idx <- which(apply(this.df[, top.keys], 1, function(keys){all(keys == "any")}))
    if(i < n.swarms){
      next.xy <- swarm.xy[[i+1]]
      segments(x0=this.xy$x, x1=next.xy$x, y0=this.xy$y, y1=next.xy$y,
               col=pw.params$linecol)
      # sapply(this.df$cancer, function(cancer){cancer.palettes[[cancer]][["light2"]]})
      segments(x0=this.xy[top.idx, "x"], x1=next.xy[top.idx, "x"],
               y0=this.xy[top.idx, "y"], y1=next.xy[top.idx, "y"], lwd=2)
    }
    points(this.xy[-top.idx, ], pch=pw.params$pwpch[-top.idx], cex=pt.cex,
           bg=pw.params$pwbg[-top.idx], col=pw.params$pwcol[-top.idx], xpd=T)
    points(this.xy[-top.idx, ], pch=pw.params$pwsym[-top.idx],
           cex=0.7*pt.cex, col="white", xpd=T)
    points(this.xy[top.idx, ], pch=pw.params$pwpch[top.idx], cex=2*pt.cex,
           bg=pw.params$pwbg[top.idx], col=pw.params$pwcol[top.idx], xpd=T)
  }
}

# Function to prepare plot area for paired swarms
prep.swarm.seq <- function(stats.list, ylims=c(min.lfc, max.lfc),
                           parmar=c(2, 3.75, 0.5, 0)){
  prep.plot.area(xlims=c(0, length(stats.list)), ylims=ylims, parmar=parmar)
  abline(h=0, lty=5)
  y.ax.at <- ylims[1]:ylims[2]
  y.ax.labels <- 2^y.ax.at
  inf.neg <- TRUE
  if(y.ax.labels[1] == -2){
    y.ax.labels[1] <- paste("\"\" <", y.ax.labels[1])
    inf.neg <- FALSE
  }
  clean.axis(2, at=y.ax.at, labels=y.ax.labels, infinite.positive=TRUE,
             infinite.negative=inf.neg, title.line=0.75, parse.labels=TRUE,
             title="Fold-enrichment vs.\nnull expectation")
}

# Parse command line arguments and options
parser <- ArgumentParser(description="Visualize trends gleaned from permutations")
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
stats <- read.table(paste(args$project_root, "/results",
                          "/VALab_germline_somatic_2024.v2.all_permutation_results",
                          "/cosmic_intogen_union_results",
                          "/cosmic_intogen_union.composite_sampling_results",
                          "/cosmic_intogen_union.composite_sampling.all_stats.tsv", sep=""),
                    header=T, sep="\t", comment.char="", check.names=F)

# Sort stats to ensure negative control always appears directly after corresponding cancer
stats <- stats[with(stats, order(cancer, germline, somatic, criteria)), ]

# Add helper flag for whether the row is a negative control point
stats$nc <- grepl("_", stats$cancer)

# Convert fold-enrichments to log2 and winsorize lower bound
stats$lfc <- log2(stats$fold)
min.lfc <- -2
stats$lfc[which(stats$lfc < min.lfc)] <- min.lfc
max.lfc <- 6
stats$lfc[which(stats$lfc > max.lfc)] <- max.lfc

# Set dimensions and parameters for all paired swarm plots
ps.height <- 2.5
ps.width <- 5

# Paired swarm of LFC by tier, all cancers
pdf(paste(args$out_prefix, "composite.lfc_by_tier.emp_vs_negctrl.pdf", sep="."),
    height=ps.height, width=ps.width)
prep.paired.swarm()
paired.swarm(stats[which(stats$criteria == "tier1"), ], at=(1:2)-0.5,
             top.keys=c("cancer", "germline", "somatic"))
paired.swarm(stats[which(stats$criteria == "tier2"), ], at=(4:5)-0.5,
             top.keys=c("cancer", "germline", "somatic"))
paired.swarm(stats[which(stats$criteria == "tier3"), ], at=(7:8)-0.5,
             top.keys=c("cancer", "germline", "somatic"))
paired.swarm(stats[which(stats$criteria == "tier4"), ], at=(10:11)-0.5,
             top.keys=c("cancer", "germline", "somatic"))
axis(1, at=c(1, 4, 7, 10), tick=F, line=0.2, labels=paste("Tier", 1:4), cex.axis=5.5/6)
dev.off()

# Paired swarm of LFC by context, all cancers
pdf(paste(args$out_prefix, "composite.lfc_by_context.emp_vs_negctrl.pdf", sep="."),
    height=ps.height, width=ps.width)
prep.paired.swarm()
paired.swarm(stats[which(stats$germline == "coding"
                         & stats$somatic == "coding"
                         & !grepl("plus", stats$criteria)), ],
             at=(1:2)-0.5, top.keys=c("cancer", "criteria"))
paired.swarm(stats[which(stats$germline == "coding"
                         & stats$somatic == "noncoding"
                         & !grepl("plus", stats$criteria)), ],
             at=(4:5)-0.5, top.keys=c("cancer", "criteria"))
paired.swarm(stats[which(stats$germline == "noncoding"
                         & stats$somatic == "coding"
                         & !grepl("plus", stats$criteria)), ],
             at=(7:8)-0.5, top.keys=c("cancer", "criteria"))
paired.swarm(stats[which(stats$germline == "noncoding"
                         & stats$somatic == "noncoding"
                         & !grepl("plus", stats$criteria)), ],
             at=(10:11)-0.5, top.keys=c("cancer", "criteria"))
sapply(0:3, function(i){
  axis(1, at=1+(3*i), tick=F, cex.axis=5/6, line=0.1,
       labels=if(i<2){"Coding"}else{"Noncoding"},
       col.axis=if(i<2){"black"}else{"gray70"})
  axis(1, at=1+(3*i), tick=F, cex.axis=5/6, line=0.9,
       labels=if(i %in% c(0, 2)){"Coding"}else{"Noncoding"},
       col.axis=if(i %in% c(0, 2)){"black"}else{"gray70"})
})
axis(1, at=0, tick=F, xpd=T, cex.axis=5.5/6, line=0.1, label="Germline:", hadj=1)
axis(1, at=0, tick=F, xpd=T, cex.axis=5.5/6, line=0.9, label="Somatic:", hadj=1)
dev.off()

# Set device params for main effect decay plots
main.width <- 3.5
main.height <- 2.25

# Main effect size decay by tier
stats.list <- lapply(c("tier1", "tier2", "tier3", "tier4"), function(tier){
  stats[which(!stats$nc & stats$criteria == tier & stats$germline == "any" & stats$somatic == "any"), ]
})
main.lfc.range <- range(do.call("rbind", stats.list)[, "lfc"])
pdf(paste(args$out_prefix, "composite.lfc_decay_by_tier.main.pdf", sep="."),
    width=main.width, height=main.height)
prep.swarm.seq(stats.list, ylims=c(floor(main.lfc.range[1]), ceiling(main.lfc.range[2])))
swarm.sequence(stats.list, top.keys=c("cancer", "germline", "somatic"))
axis(1, at=0.5, tick=F, line=0.1, labels="Tier 1", col.axis="white", cex.axis=5.5/6)
axis(1, at=(1:3)+0.5, tick=F, line=0.1, labels=paste("Tier", 2:4), cex.axis=5.5/6)
dev.off()

# Main effect size decay by context
stats.list <- list(stats[which(!stats$nc & stats$criteria == "any" & stats$germline == "coding" & stats$somatic == "coding"), ],
                   stats[which(!stats$nc & stats$criteria == "any" & stats$germline == "coding" & stats$somatic == "noncoding"), ],
                   stats[which(!stats$nc & stats$criteria == "any" & stats$germline == "noncoding" & stats$somatic == "coding"), ],
                   stats[which(!stats$nc & stats$criteria == "any" & stats$germline == "noncoding" & stats$somatic == "noncoding"), ])
main.lfc.range <- range(do.call("rbind", stats.list)[, "lfc"])
pdf(paste(args$out_prefix, "composite.lfc_decay_by_context.main.pdf", sep="."),
    width=main.width, height=main.height)
prep.swarm.seq(stats.list, ylims=c(floor(main.lfc.range[1]), ceiling(main.lfc.range[2])))
swarm.sequence(stats.list, top.keys=c("cancer", "criteria"))
sapply(0:3, function(i){
  axis(1, at=i+0.5, tick=F, cex.axis=4.5/6, line=-0.9,
       labels=if(i<2){"Coding"}else{"Noncoding"},
       col.axis=if(i<2){"black"}else{"gray70"})
  axis(1, at=i+0.5, tick=F, cex.axis=4.5/6, line=0.1,
       labels=if(i %in% c(0, 2)){"Coding"}else{"Noncoding"},
       col.axis=if(i %in% c(0, 2)){"black"}else{"gray70"})
})
axis(1, at=0, tick=F, xpd=T, cex.axis=5.5/6, line=-0.9, label="Germline:", hadj=1)
axis(1, at=0, tick=F, xpd=T, cex.axis=5.5/6, line=0.1, label="Somatic:", hadj=1)
dev.off()

# Supplementary effect size decay by tier
# TODO: implement this

# Supplementary effect size decay by context
# TODO: implement this


