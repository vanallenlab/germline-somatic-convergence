# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Generate all plots for figure 2 (permutation testing) and its related supplements

# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence
export CONFIG=$CODEDIR/manuscript/plotting/config.R
export OUTDIR=convergence_paper/figures/permutation_figure/assets

# Make output directory
if ! [ -e $OUTDIR ]; then mkdir $OUTDIR; fi

# Summarize raw results from permutation tests
$CODEDIR/manuscript/plotting/figure2/summarize_permutation_stats.R \
  --config $CONFIG \
  --out-prefix $OUTDIR/perm_stats

# Plot permutation statistics across different null models
# (For supplement, but good to show that composite method is most conservative)
# TODO: implement this



# Switch back to original execution directory
# (Only necessary if this script is called directly from the command line)
cd -

