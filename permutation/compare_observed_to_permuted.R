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
strata.labels <- c("tier1" = "Tier 1 (same gene)",
                   "tier2" = "Tier 2 (PPI known func.)",
                   "tier2plus" = "Tiers 1 & 2",
                   "tier3" = "Tier 3 (PPI unknown func.)",
                   "tier3plus" = "Tiers 1, 2, and 3",
                   "tier4" = "Tier 4 (same complex)",
                   "any" = "Any convergence criterion")
cancer.colors <- c("any" = "#404040",
                   "breast" = "#E67F88",
                   "colorectal" = "#C097CF",
                   "lung" = "#E8AE6B",
                   "prostate" = "#84BAD9",
                   "renal" = "#65B360")

# # DEV:
# observed.in <- "~/scratch/test.observed_counts.tsv"
# permuted.in <- "~/scratch/perm_test.results.postprocessed.tsv"
# out.prefix <- "~/scratch/convergence_test"

# Make vector for expected column names
tiers <- paste("tier", 1:4, sep="")
origins <- c("germline", "somatic")
contexts <- c("coding", "noncoding")
strata.names <- as.vector(unlist(sapply(contexts, function(germline_context){
  sapply(contexts, function(somatic_context){
    paste("germline_", germline_context, ".somatic_", somatic_context, ".",
          tiers, sep="")
  })
})))

# Helper function for filling in missing summary data from observed counts
fill.summary.data <- function(df){
  # Make summary columns for all contexts
  for(tier in tiers){
    ac <- apply(df[, grep(tier, colnames(df))], 1, sum, na.rm=T)
    df[, paste("germline_any.somatic_any", tier, sep=".")] <- ac
  }

  # Make summary columns for multi-tier criteria
  strata.bases <- gsub(".tier1", "", colnames(df)[grep("tier1", colnames(df))])
  t2p <- df[, grep("tier1", colnames(df))] + df[, grep("tier2", colnames(df))]
  df[, paste(strata.bases, "tier2plus", sep=".")] <- t2p
  t3p <- df[, grep("tier2plus", colnames(df))] + df[, grep("tier3", colnames(df))]
  df[, paste(strata.bases, "tier3plus", sep=".")] <- t3p
  t4p <- df[, grep("tier3plus", colnames(df))] + df[, grep("tier4", colnames(df))]
  df[, paste(strata.bases, "any", sep=".")] <- t4p

  # Make summary rows for pooled result across all cancers
  # and, separately, all negative control phenotype
  prim.cancers <- sort(unique(df$cancer[grep("_", df$cancer, invert=T)]))
  nc.phenos <- sort(unique(sapply(setdiff(df$cancer, prim.cancers), function(ps){
    paste(unlist(strsplit(ps, split="_"))[-1], collapse="_")
  })))
  ac <- as.numeric(apply(df[which(df$cancer %in% prim.cancers), -1],
                         2, sum, na.rm=T))
  df <- as.data.frame(rbind(df, c("any", ac)))
  for(pheno in nc.phenos){
    if(length(unique(df$cancer[grep(pheno, df$cancer)])) > 1){
      ac <- as.numeric(apply(df[grep(pheno, df$cancer), -1], 2,
                             function(v){sum(as.numeric(v), na.rm=T)}))
      df <- as.data.frame(rbind(df, c(paste("any", pheno, sep="_"), ac)))
    }
  }
  any.nc.idx <- intersect(grep(paste(prim.cancers, collapse="|"), df$cancer),
                          grep(paste(nc.phenos, collapse="|"), df$cancer))
  ac <- as.numeric(apply(df[any.nc.idx, -1], 2,
                         function(v){sum(as.numeric(v), na.rm=T)}))
  df <- as.data.frame(rbind(df, c("any_negative_control", ac)))
  df[, -1] <- apply(df[, -1], 2, as.numeric)

  return(df)
}

# Load observed data as data.frame and fill missing summary counts
obs.df <- read.table(observed.in, sep="\t", header=F, col.names=c("cancer", strata.names))
obs.df <- fill.summary.data(obs.df)
cancers <- sort(unique(obs.df$cancer))
all.strata.names <- colnames(obs.df)[-1]

# Load permuted data as data.frame
perm.df <- read.table(permuted.in, sep="\t", header=T)

# Ensure strata match between observed and permuted inputs
if(length(intersect(all.strata.names, colnames(perm.df))) != length(all.strata.names)){
  cat(paste("Some strata (column names) differ between observed and permuted",
            "inputs. Please check your input files and try running again."))
  quit(save=F, status=1)
}

# Run one comparison for each cancer type and strata
res <- do.call("rbind", lapply(cancers, function(cancer){
  as.data.frame(do.call("rbind", lapply(all.strata.names, function(strata){
    # Permutation test
    obs.val <- as.numeric(obs.df[which(obs.df$cancer == cancer), strata])
    perm.vals <- as.numeric(perm.df[which(perm.df$cancer == cancer), strata])
    exp.mean <- mean(perm.vals, na.rm=T)
    exp.sd <- sd(perm.vals, na.rm=T)
    exp.ci <- quantile(perm.vals, c(0.025, 0.975))
    fold <- obs.val / exp.mean
    if(all(exp.ci == 0)){
      fold.ci <- c(0, Inf)
    }else if(min(exp.ci) <= 0){
      fold.ci <- c(obs.val / max(exp.ci), Inf)
    }else{
      fold.ci <- sort(obs.val / exp.ci)
    }
    perms.gt <- length(which(perm.vals >= obs.val))
    n.perms <- length(perm.vals)
    p.perm <- (perms.gt + 1) / (n.perms + 1)
    p.pois <- ppois(obs.val, exp.mean, lower.tail=FALSE)
    # To estimate P-value beyond the number of permutations, we use a Poisson null
    if(perms.gt == 0){
      p <- p.pois
    }else{
      # Otherwise, we take the more conservative P-value between Poisson & permuted nulls
      p <- max(p.perm, p.pois, na.rm=T)
    }

    # Plot obs/exp histogram
    highlight.color <- if(cancer %in% names(cancer.colors)){cancer.colors[cancer]}else{"gray40"}
    pdf(paste(out.prefix, cancer, strata, "pdf", sep="."),
        height=2.9, width=2.3)
    RLCtools::density.w.outliers(perm.vals, style="hist", min.bin.width=1,
                                 xlims=range(c(perm.vals, obs.val)),
                                 x.title="Germline:somatic gene pairs",
                                 x.title.line=0, x.label.units="counts",
                                 add.y.axis=FALSE, color="gray70", border="gray70",
                                 parmar=c(2, 0.5, 6, 0.5))
    if(mean(perm.vals) <= mean(par("usr")[1:2])){
      pn.lab.at <- ceiling(mean(perm.vals))
      pn.lab.pos <- 4
    }else{
      pn.lab.at <- floor(mean(perm.vals))
      pn.lab.pos <- 2
    }
    text(x=pn.lab.at, y=0.8*par("usr")[4], pos=pn.lab.pos, cex=5/6, font=3,
         label="Permuted null", col="gray70")
    segments(x0=obs.val, x1=obs.val, y0=0, y1=0.5*par("usr")[4],
             col=highlight.color, lwd=6)
    points(x=obs.val, y=0.5*par("usr")[4], pch=25, col=highlight.color,
           bg=highlight.color, cex=2, xpd=T)
    text(x=obs.val, y=0.5*par("usr")[4], pos=3, labels="Obs.",
         col=highlight.color, cex=5/6, xpd=T)
    strata.parts <- lapply(unlist(strsplit(strata, split=".", fixed=T)),
                           function(s){unlist(strsplit(s, split="_"))})
    if(grepl("_", cancer)){
      cparts <- unlist(strsplit(cancer, split="_"))
      l6 <- paste(cparts[1], "vs.", paste(cparts[-1], collapse=" "))
    }else{
      l6 <- paste(cancer, "cancer")
    }
    mtext(3, line=6-0.9, text=title.case(l6, case="sentence"), cex=5/6)
    l5 <- sub("any", "all", paste(c(strata.parts[[1]][2:1], "Variants"), collapse=" "))
    mtext(3, line=5-0.75, text=title.case(l5, case="sentence"), cex=5/6)
    l4 <- sub("any", "all", paste(c(strata.parts[[2]][2:1], "Alterations"), collapse=" "))
    mtext(3, line=4-0.6, text=title.case(l4, case="sentence"), cex=5/6)
    l3 <- strata.labels[strata.parts[[3]][1]]
    mtext(3, line=3-0.45, text=title.case(l3, case="sentence"), cex=5/6)
    l2 <- paste(obs.val, "observed vs.", round(exp.mean, 1), "expected")
    mtext(3, line=2-0.3, text=l2, cex=5/6)
    if(is.na(fold) | is.infinite(fold)){
      l1 <- "Undefined fold-change"
    }else{
      l1 <- paste(round(fold, 1), "-fold ", if(fold>=1){"enriched"}else{"depleted"},
                  " (95% CI: ", paste(round(fold.ci, 1), collapse="-"),
                  ")", sep="")
    }
    mtext(3, line=1-0.15, text=l1, cex=5/6)
    l0 <- RLCtools::format.pval(p, min.neg.log10.p=10)
    mtext(3, line=0, text=l0, cex=5/6)
    dev.off()

    # Return permutation res
    if(is.nan(fold)){fold <- 1}
    c(cancer,
      paste(setdiff(strata.parts[[1]], "germline"), sep="_"),
      paste(setdiff(strata.parts[[2]], "somatic"), sep="_"),
      paste(strata.parts[[3]], collapse="_"),
      obs.val, exp.mean, exp.sd, as.numeric(exp.ci), fold, as.numeric(fold.ci),
      p.perm, p.pois, p)
  })))
}))
colnames(res) <- c("cancer", "germline", "somatic", "criteria", "observed",
                   "expected", "expected_sd", "expected_lower_95pct",
                   "expected_upper_95pct", "fold", "fold_lower_95pct",
                   "fold_upper_95pct", "p_perm", "p_poisson", "p")
write.table(res, paste(out.prefix, "all_stats.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
