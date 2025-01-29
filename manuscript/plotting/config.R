#!/usr/bin/env Rscript

# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Configuration file for germline:somatic convergence figures

# Recommended implementation: source this file within any plotting Rscript

# Disable undesired variable formatting
options(stringsAsFactors=FALSE, scipen=1000)

# Tier shapes and shading
tier.pch <- list("tier1" = 22,
                 "tier2" = 24,
                 "tier2plus" = 24,
                 "tier3" = 23,
                 "tier3plus" = 23,
                 "tier4" = 21,
                 "any" = 21)
tier.shades <- list("tier1" = "dark1",
                    "tier2" = "main",
                    "tier2plus" = "main",
                    "tier3" = "light1",
                    "tier3plus" = "light1",
                    "tier4" = "light2",
                    "any" = "light2")

# Context shading
context.shades <- list("coding" = "dark1",
                       "any" = "main",
                       "noncoding" = "light2")

# Cancer color palettes
pancan.colors <- list("dark2" = "#202020",
                      "dark1" = "#303030",
                      "main" = "#404040",
                      "light1" = "#666666",
                      "light2" = "#8C8C8C")
breast.colors <- list("dark2" = "#734044",
                      "dark1" = "#AD5F66",
                      "main" = "#E67F88",
                      "light1" = "#EB99A0",
                      "light2" = "#F0B2B8")
colorectal.colors <- list("dark2" = "#604C68",
                          "dark1" = "#90719B",
                          "main" = "#C097CF",
                          "light1" = "#CDACD9",
                          "light2" = "#D9C1E2")
lung.colors <- list("dark2" = "#745736",
                    "dark1" = "#AE8350",
                    "main" = "#E8AE6B",
                    "light1" = "#EDBE89",
                    "light2" = "#F1CEA6")
prostate.colors <- list("dark2" = "#425D6D",
                        "dark1" = "#638BA3",
                        "main" = "#84BAD9",
                        "light1" = "#9DC8E1",
                        "light2" = "#B5D6E8")
renal.colors <- list("dark2" = "#335A30",
                     "dark1" = "#4C8648",
                     "main" = "#65B360",
                     "light1" = "#84C280",
                     "light2" = "#A3D1A0")
cancer.palettes <- list("any" = pancan.colors,
                        "breast" = breast.colors,
                        "colorectal" = colorectal.colors,
                        "lung" = lung.colors,
                        "prostate" = prostate.colors,
                        "renal" = renal.colors)
cancer.colors <- lapply(cancer.palettes, function(p){p[["main"]]})

# General declarative variables and variable names
# Note that we exclude pan-cancer from the main `cancers` vector because it will
# be treated separately in most analyses
cancers <- setdiff(names(cancer.colors), "any")
n.cancers <- length(cancers)
origins <- c("germline", "somatic")
contexts <- c("coding", "noncoding")
