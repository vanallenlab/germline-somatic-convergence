#!/usr/bin/env Rscript

# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Compare observed germline:somatic convergence to expected convergence based on permutation

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

# Helper function to coerce all words in a character vector to title case
titleCase <- function(x){
  as.character(sapply(x, function(s){
    paste(toupper(substr(s, 1, 1)),
          tolower(substr(s, 2, nchar(s))),
          sep="")
  }))
}

# Helper function to format P-values
format.pval <- function(p, nsmall=2, max.decimal=3, equality="=", min.neg.log10.p=100){
  if(-log10(p)>min.neg.log10.p){
    bquote(italic(P) %~~% 0)
  }else if(ceiling(-log10(p)) > max.decimal){
    parts <- unlist(strsplit(format(p, scientific=T), split="e"))
    base <- gsub(" ", "", formatC(round(as.numeric(parts[1]), nsmall), digits=1+nsmall), fixed=T)
    exp <- gsub(" ", "", as.numeric(parts[2]), fixed=T)
    bquote(italic(P) ~ .(equality) ~ .(base) ~ "x" ~ 10 ^ .(exp))
  }else{
    bquote(italic(P) ~ .(equality) ~ .(formatC(round(p, max.decimal), digits=max.decimal)))
  }
}

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
    exp.ci <- quantile(perm.vals, c(0.025, 0.975))
    fold <- obs.val / exp.mean
    perms.gt <- length(which(perm.vals >= obs.val))
    n.perms <- length(perm.vals)
    p <- (perms.gt + 1) / (n.perms + 1)

    # Plot obs/exp histogram
    pdf(paste(out.prefix, cancer, strata, "pdf", sep="."),
        height=3.5, width=3.5)
    par(mar=c(2.25, 0.5, 6, 0.5), bty="n")
    highlight.color <- "#F89820"
    h <- hist(perm.vals, breaks=seq(0, max(perm.vals)+1), plot=F)$density
    h <- h / max(h, na.rm=T)
    xmax <- 1.1 * max(c(obs.val, perm.vals), na.rm=T)
    plot(NA, NA, xlim=c(0, xmax), ylim=c(0, 1),
         xaxt="n", xlab="", yaxt="n", ylab="", yaxs="i")
    rect(xleft=(1:length(h))-1, xright=1:length(h), ybottom=0, ytop=h,
         col="#02679A", border="#003354")
    text(x=mean(perm.vals), y=0.9, pos=4, cex=5/6, font=3,
         label="Permuted Null", col="#003354")
    segments(x0=obs.val, x1=obs.val, y0=0, y1=2/3, col=highlight.color, lwd=4)
    points(x=obs.val, y=2/3, pch=25, col=highlight.color, bg=highlight.color, cex=2)
    text(x=obs.val, y=2/3, pos=3, labels="Obs.", col=highlight.color, cex=5/6)
    axis(1, at=c(0, 10e10), tck=0, labels=NA)
    axis(1, at=axTicks(1), tck=-0.025, labels=NA)
    axis(1, at=axTicks(1), tick=F, line=-0.65, cex.axis=5/6)
    mtext(1, line=1.25, text="Germline:Somatic Gene Pairs")
    strata.parts <- lapply(unlist(strsplit(strata, split=".", fixed=T)), function(s){unlist(strsplit(s, split="_"))})
    mtext(3, line=5, text=paste(titleCase(cancer), "Cancer"), font=2)
    mtext(3, line=4, text=sub("Any", "All", paste(titleCase(c(strata.parts[[1]][2:1], "Variants")), collapse=" ")))
    mtext(3, line=3, text=sub("Any", "All", paste(titleCase(c(strata.parts[[2]][2:1], "Variants")), collapse=" ")))
    mtext(3, line=2, text=strata.labels[strata.parts[[3]][1]])
    mtext(3, line=1, text=paste(obs.val, "Observed vs.", round(exp.mean, 1), "Expected"))
    p.label <- if(p > 1/(n.perms+1)){format.pval(p)}else{format.pval(1/n.perms, equality="<")}
    mtext(3, line=0, text=p.label)
    dev.off()

    # Return permutation res
    c(cancer,
      paste(setdiff(strata.parts[[1]], "germline"), sep="_"),
      paste(setdiff(strata.parts[[2]], "somatic"), sep="_"),
      paste(strata.parts[[3]], collapse="_"),
      obs.val, exp.mean, exp.ci, fold, p)
  })))
}))
colnames(res) <- c("cancer", "germline", "somatic", "criteria", "observed", "expected",
                   "expected_lower_95pct", "expected_upper_95pct", "fold","p")
write.table(res, paste(out.prefix, "all_stats.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
