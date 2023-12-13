#!/usr/bin/env Rscript

# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Script to shuffle a gene list based on gene-specific weights

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)
require(argparse, quietly=TRUE)

# Helper function to convert any character into a seed
# Taken from the "TeachingDemos" library
# (but not imported to avoid needing a custom Docker)
char2seed <- function(x, set=TRUE, ...){

  tmp <- c(0:9,0:25,0:25)
  names(tmp) <- c(0:9,letters,LETTERS)

  x <- gsub("[^0-9a-zA-Z]","",as.character(x))

  xsplit <- tmp[ strsplit(x,'')[[1]] ]

  seed <- sum(rev( 7^(seq(along=xsplit)-1) ) * xsplit)
  seed <- as.integer( seed %% (2^31-1) )

  if(set){
    set.seed(seed, ...)
    return(invisible(seed))
  } else {
    return(seed)
  }
}


# Read arguments
parser <- ArgumentParser(description="Weighted shuffling of a gene list")
parser$add_argument("--tsv-in", metavar=".tsv", type="character", required=TRUE,
                    help="two-column .tsv of gene symbols and weights")
parser$add_argument("--seed", metavar="string or file", type="character", required=TRUE,
                    help=paste("seed string for shuffling. Can also provide a ",
                               "file containing one or more seeds, which will",
                               "be used to create multiple output files"))
parser$add_argument("--cancer", metavar="string", type="character",
                    help=paste("(Optional) cancer type to extract from --tsv-in. ",
                               "If specified, --tsv-in must have a third column ",
                               "for cancer type."))
parser$add_argument("--eligible-genes", metavar="file", type="character",
                    help=paste("(Optional) file with gene symbols to include"))
parser$add_argument("--outfile", metavar=".txt", type="character", default="stdout",
                    help=paste("path to output file [default: stdout]. If",
                               "a file is passed as --seed, the value of --outfile",
                               "will be ignored and one outfile per line in --seed",
                               "will be written to the execution directory."))
args <- parser$parse_args()

# # DEV (single):
# args <- list("tsv_in" = "~/scratch/gene_weights.uniform.tsv",
#              "seed" = "0btgecc1",
#              "eligible_genes" = "~/scratch/gencode.v44.autosomal.protein_coding.genes.list",
#              "outfile" = "~/scratch/shuffled_genes.dev.list")

# # DEV (multi):
# args <- list("tsv_in" = "~/Dropbox (Partners HealthCare)/VanAllen/VALab_germline_somatic_2023/other_data/permutation_weights/gene_weights.coding_nonsynonymous.tsv",
#              "seed" = "~/scratch/test_seeds.list",
#              "outfile" = "~/scratch/shuffled_genes.dev.list")

# Read genes and weights
x <- read.table(args$tsv_in, header=T, sep="\t", comment.char="")
colnames(x) <- c("gene", "weight", "cancer")[1:ncol(x)]

# Subset to a specific cancer type, if optioned
if(!is.null(args$cancer)){
  x <- x[which(x$cancer == args$cancer), 1:2]
}

# Subset to --eligible-genes, if optioned
if(!is.null(args$eligible_genes)){
  elig <- read.table(args$eligible_genes, header=F)[, 1]
  x <- x[which(x[, 1] %in% elig), ]
}

# Ensure genes begin sorted alphabetically and no weights are missing
x <- x[which(!is.na(x$weight)), ]
x <- x[order(x$gene), ]

# Normalize weights to reflect sampling probabilities
x$prob <- x$weight / sum(x$weight, na.rm=T)

# Read seeds
if(file.exists(args$seed)){
  seeds <- read.table(args$seed, header=F)[, 1]
}else{
  seeds <- args$seed
}

# Shuffle once per seed
sapply(seeds, function(seed){

  # Set seed
  char2seed(seed)

  # Sample genes based on their weights
  new.order <- sample(1:nrow(x), size=nrow(x), replace=FALSE, prob=x$prob)

  # Write shuffled gene list to --outfile unless --seed was a file of seeds
  if(length(seeds) > 1){
    outfile <- paste("shuffled", seed, "txt", sep=".")
  }else{
    if(args$outfile %in% c("stdout", "/dev/stdout", "-")){
      outfile <- stdout()
    }else{
      outfile <- args$outfile
    }
  }
  write.table(x[new.order, "gene"], outfile, col.names=F, row.names=F, quote=F)
})

