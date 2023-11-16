# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Bash code to identify all germline:somatic overlaps among candidate genes


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germ_som_pilot_2023/germline-somatic-exploration-2023


# Prep output directory
for subdir in results results/find_pairs; do
  if [ ! -e $subdir ]; then
    mkdir $subdir
  fi
done


# Run find_pairs.py for all combination of coding|noncoding and germline|somatic
for cancer in breast colorectal lung prostate renal; do
  for germ_context in coding noncoding; do
    for som_context in coding noncoding; do
      $CODEDIR/analysis/find_pairs.py \
        --germline gene_lists/germline_$germ_context/$cancer.germline.$germ_context.genes.list \
        --somatic gene_lists/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
        --cellchat-db other_data/cellchat_db_formatted.csv \
        --ppi-db other_data/ebi_intact.all_interactions.tsv \
        --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
        --out-tsv results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.pairs.tsv
    done
  done
done


# Concatenate all results into a single table
for cancer in breast colorectal lung prostate renal; do
  for germ_context in coding noncoding; do
    for som_context in coding noncoding; do
      fgrep -v "#" \
        results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.pairs.tsv \
      | sed '/^$/d' \
      | awk -v cancer=$cancer -v gc=$germ_context -v sc=$som_context -v OFS="\t" \
        '{ print cancer, $1, gc, $2, sc, $3 }'
    done
  done
done \
| sort -Vk1,1 -k2,2V -k4,4V -k3,3V -k5,5V \
| cat \
  <( echo -e "#cancer\tgermline_gene\tgermline_context\tsomatic_gene\tsomatic_context\tcriteria" ) \
  - \
> results/VALab_germline_somatic_2023.pilot.gene_pairs.annotated.tsv
