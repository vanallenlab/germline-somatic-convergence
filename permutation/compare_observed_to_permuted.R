#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Compare observed germline:somatic convergence to expected convergence based on permutation

# Load libraries
require(RLCtools)

# Set parameters
options(scipen=1000, stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
observed.in <- as.character(args[1])
permuted.in <- as.character(args[2])
out.prefix <- as.character(args[3])
strata.labels <- c("any" = "Any Convergence Criterion",
                   "same" = "Identical Genes",
                   "ligand" = "Ligand/Receptor Pairs",
                   "known" = "Protein-Protein Interactions",
                   "protein" = "Same Protein Complex")
cancer.colors <- c("all" = "#404040",
                   "breast" = "#FF7A7A",
                   "colorectal" = "#BA71D6",
                   "lung" = "#FF7C1F",
                   "prostate" = "#66B9E8",
                   "renal" = "#3CB333")

# # DEV:
# observed.in <- "~/scratch/test.observed_counts.tsv"
# permuted.in <- "~/scratch/germ_som_convergence.permutation_results.uniform_weighting.n100.txt.gz"
# out.prefix <- "~/scratch/convergence_test"

# Make vector for expected column names
cancers <- c("breast", "colorectal", "lung", "prostate", "renal")
criteria <- c("any", "same_gene", "ligand_receptor", "known_ppi", "protein_complex")
origins <- c("germline", "somatic")
contexts <- c("coding", "noncoding")
strata.names <- c(
  paste("germline_any.somatic_any", criteria, sep="."),
  as.vector(unlist(sapply(contexts, function(germline_context){
    sapply(contexts, function(somatic_context){
      paste("germline_", germline_context, ".somatic_", somatic_context, ".",
            criteria, sep="")
    })
  })))
)

# Load observed data as data.frame
obs.df <- read.table(observed.in, sep="\t", header=F, col.names=c("cancer", strata.names))

# Load permuted data as data.frame
perm.df <- read.table(permuted.in, sep="\t", header=F,
                      col.names=c("cancer", "perm_idx.major", "perm_idx.minor", strata.names))

# Run one comparison for each cancer type and strata
res <- do.call("rbind", lapply(c("all", cancers), function(cancer){
  as.data.frame(do.call("rbind", lapply(strata.names, function(strata){
    # Permutation test
    obs.val <- as.numeric(obs.df[which(obs.df$cancer == cancer), strata])
    perm.vals <- as.numeric(perm.df[which(perm.df$cancer == cancer), strata])
    exp.mean <- mean(perm.vals, na.rm=T)
    exp.se <- sd(perm.vals, na.rm=T) / sqrt(length(perm.vals))
    exp.ci <- quantile(perm.vals, c(0.025, 0.975))
    fold <- obs.val / exp.mean
    fold.ci <- sort(obs.val / (exp.mean + (qnorm(c(0.025, 0.975)) * exp.se)))
    fold.ci[1] <- max(c(0, fold.ci[1]))
    perms.gt <- length(which(perm.vals >= obs.val))
    n.perms <- length(perm.vals)
    p <- (perms.gt + 1) / (n.perms + 1)

    # Plot obs/exp histogram
    highlight.color <- cancer.colors[cancer]
    pdf(paste(out.prefix, cancer, strata, "pdf", sep="."),
        height=3.5, width=3.5)
    RLCtools::density.w.outliers(perm.vals, style="hist", min.bin.width=1,
                                 xlims=range(c(perm.vals, obs.val)),
                                 x.title="Germline:Somatic Gene Pairs",
                                 x.title.line=0.25, x.label.units="counts",
                                 add.y.axis=FALSE, color="gray70", border="gray70",
                                 parmar=c(2.25, 0.5, 7, 0.5))
    text(x=ceiling(mean(perm.vals)), y=0.9*par("usr")[4], pos=4, cex=5/6, font=3,
         label="Permuted null", col="gray70")
    segments(x0=obs.val, x1=obs.val, y0=0, y1=0.5*par("usr")[4],
             col=highlight.color, lwd=6)
    points(x=obs.val, y=0.5*par("usr")[4], pch=25, col=highlight.color,
           bg=highlight.color, cex=2, xpd=T)
    text(x=obs.val, y=0.5*par("usr")[4], pos=3, labels="Obs.",
         col=highlight.color, cex=5/6, xpd=T)
    strata.parts <- lapply(unlist(strsplit(strata, split=".", fixed=T)),
                           function(s){unlist(strsplit(s, split="_"))})
    l6 <- paste(cancer, "Cancer")
    mtext(3, line=6, text=title.case(l6, case="sentence"), font=2)
    l5 <- sub("any", "all", paste(c(strata.parts[[1]][2:1], "Variants"), collapse=" "))
    mtext(3, line=5, text=title.case(l5, case="sentence"))
    l4 <- sub("any", "all", paste(c(strata.parts[[2]][2:1], "Variants"), collapse=" "))
    mtext(3, line=4, text=title.case(l4, case="sentence"))
    l3 <- strata.labels[strata.parts[[3]][1]]
    mtext(3, line=3, text=title.case(l3, case="sentence"))
    l2 <- paste(obs.val, "observed vs.", round(exp.mean, 1), "expected")
    mtext(3, line=2, text=l2)
    if(is.na(fold) | is.infinite(fold)){
      l1 <- "Undefined fold-change"
    }else{
      l1 <- paste(round(fold, 1), "-fold ", if(fold>=1){"enriched"}else{"depleted"},
                  " (95% CI: ", paste(round(fold.ci, 1), collapse=" to "),
                  ")", sep="")
    }
    mtext(3, line=1, text=l1)
    l0 <- if(p > 1/(n.perms+1)){
      RLCtools::format.pval(p)
    }else{
      RLCtools::format.pval(1/n.perms, equality="<")
    }
    mtext(3, line=0, text=l0)
    dev.off()

    # Return permutation res
    c(cancer,
      paste(setdiff(strata.parts[[1]], "germline"), sep="_"),
      paste(setdiff(strata.parts[[2]], "somatic"), sep="_"),
      paste(strata.parts[[3]], collapse="_"),
      obs.val, exp.mean, exp.ci, fold, fold.ci, p)
  })))
}))
colnames(res) <- c("cancer", "germline", "somatic", "criteria", "observed",
                   "expected", "expected_lower_95pct", "expected_upper_95pct",
                   "fold", "fold_lower_95pct", "fold_upper_95pct", "p")
write.table(res, paste(out.prefix, "all_stats.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
