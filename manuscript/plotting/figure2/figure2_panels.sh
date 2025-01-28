# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Generate all plots for figure 2 (permutation testing) and its related supplements

# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence
export CONFIG=$CODEDIR/manuscript/plotting/config.R
export OUTDIR=convergence_paper/figures/trends_figure/assets

# Make output directory
if ! [ -e $OUTDIR ]; then mkdir $OUTDIR; fi

# Plot decay of enrichment vs. tier & context
$CODEDIR/manuscript/plotting/figure2/permutation_trends.R \
  --config $CONFIG \
  --out-prefix $OUTDIR/perm_trends

# Plot GSEA results
# TODO: implement this

# Switch back to original execution directory
# (Only necessary if this script is called directly from the command line)
cd -

