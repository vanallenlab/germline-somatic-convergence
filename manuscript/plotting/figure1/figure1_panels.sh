# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Generate all plots for figure 1 and its related supplements

# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence
export CONFIG=$CODEDIR/manuscript/plotting/config.R

# Make output directory
mkdir convergence_paper/figures/summary_figure/assets

# Barplots summarizing input gene lists
$CODEDIR/manuscript/plotting/figure1/gene_set_sizes.R \
  --config $CONFIG \
  --out-prefix convergence_paper/figures/summary_figure/assets/gene_set_sizes

# Barplots summarizing empirical pairs and hub genes
$CODEDIR/manuscript/plotting/figure1/describe_pairs.R \
  --config $CONFIG \
  --out-prefix convergence_paper/figures/summary_figure/assets/gene_set_sizes

# Plot GSEA results
# TODO: implement this

# Switch back to original execution directory
# (Only necessary if this script is called directly from the command line)
cd -
