#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Barplots of gene set sizes by cancer type and context

# Set parameters and load libraries
options(scipen=1000, stringsAsFactors=FALSE)
require(argparse)
require(RLCtools)

# Parse command line arguments and options
parser <- ArgumentParser(description="Summarize input gene sets")
parser$add_argument("--config", metavar=".R", type="character", required=TRUE,
                    help="path to config.R file")
parser$add_argument("--project-root", metavar="path", type="character",
                    help="root directory for all project data [default: pwd]")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    default="./gene_set_sizes", help="file prefix for all plots")
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

# Read all gene lists into memory
gsets <- lapply(origins, function(origin){
  x1 <- lapply(contexts, function(context){
    x2 <- lapply(cancers, function(cancer){
      read.table(paste(args$project_root, "/gene_lists/", origin, "_", context, "/",
                       cancer, ".", origin, ".", context, ".genes.list", sep=""),
                 header=F)[, 1]
    })
    names(x2) <- cancers
    return(x2)
  })
  names(x1) <- contexts
  return(x1)
})
names(gsets) <- origins

# Simple barplot of genes per cancer by context
bar.xbuf <- 0.15
ymax <- max(sapply(gsets, function(x1){sapply(x1, function(x2){sapply(x2, length)})}))
sapply(origins, function(origin){
  sapply(contexts, function(context){
    v <- sapply(gsets[[origin]][[context]], length)
    pdf(paste(args$out_prefix, "gene_count_bars", origin, context, "pdf", sep="."),
        height=0.9, width=1.3)
    prep.plot.area(xlims=c(-bar.xbuf, n.cancers), ylims=c(0, 1.25*ymax),
                   parmar=c(0.3, 2.3, 0.3, 0))
    axis(1, at=c(0, n.cancers), tck=0, labels=NA)
    # v.col <- sapply(cancers, function(cancer){cancer.palettes[[cancer]][[context.shades[[context]]]]})
    v.col <- unlist(cancer.colors[cancers])
    rect(xleft=(1:n.cancers)-1+bar.xbuf, xright=(1:n.cancers)-bar.xbuf,
         ybottom=rep(0, n.cancers), ytop=v, col=v.col)
    text(x=(1:n.cancers)-0.5, y=v+(0.075*diff(par("usr")[3:4])),
         labels=prettyNum(v, big.mark=","), cex=4/6, col=v.col)
    clean.axis(2, infinite=TRUE, label.units="count", title="Genes",
               title.line=0.4, tck=-0.035)
    dev.off()
  })
})

# Upset of # of cancer types per gene; one plot per context
# TODO: implement this

# Upset of # of contexts per gene; one plot per cancer type
# TODO: implement this

# Read filtered COSMIC and intOGen lists into memory
# TODO: implement this

# Upset of intOGen vs. COSMIC
# TODO: implement this
