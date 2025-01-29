# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Generate all plots for figure 1 (data description) and its related supplements

# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence
export CONFIG=$CODEDIR/manuscript/plotting/config.R
export OUTDIR=convergence_paper/figures/summary_figure/assets

# Make output directory
if ! [ -e $OUTDIR ]; then mkdir $OUTDIR; fi

# Barplots summarizing input gene lists
$CODEDIR/manuscript/plotting/figure1/gene_set_sizes.R \
  --config $CONFIG \
  --out-prefix $OUTDIR/gene_set_sizes

# Barplots summarizing empirical pairs and hub genes
$CODEDIR/manuscript/plotting/figure1/describe_pairs.R \
  --config $CONFIG \
  --out-prefix $OUTDIR/pair_analysis

# Summarize raw results from permutation tests across different null models
# (3/4 of these will be for supplement, but good to show that composite method 
#  is best calibrated)
for model in uniform bayesian expression composite; do
  $CODEDIR/manuscript/plotting/figure1/summarize_permutation_stats.R \
    --config $CONFIG \
    --model $model \
    --out-prefix $OUTDIR/perm_stats
done

# Switch back to original execution directory
# (Only necessary if this script is called directly from the command line)
cd -
