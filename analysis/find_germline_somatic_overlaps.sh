# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Bash code to identify all germline:somatic overlaps among candidate genes


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence


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
      if [ $som_context != "coding" ]; then
        # Only one source of somatic noncoding hits
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$cancer.germline.$germ_context.genes.list \
          --somatic gene_lists/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.pairs.tsv
      else
        # For somatic coding hits, we are evaluating four possible definitions, as follows:
        # 1. Union of intOGen + COSMIC (this is the default)
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$cancer.germline.$germ_context.genes.list \
          --somatic gene_lists/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.union.pairs.tsv
        # 2. Intersection of intOGen + COSMIC
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$cancer.germline.$germ_context.genes.list \
          --somatic other_data/cosmic_intogen_intersection/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.intersection.pairs.tsv
        # 3. COSMIC only
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$cancer.germline.$germ_context.genes.list \
          --somatic other_data/cosmic_genes/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.cosmic_only.pairs.tsv
        # 4. intOGen only
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$cancer.germline.$germ_context.genes.list \
          --somatic other_data/intogen_genes/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.intogen_only.pairs.tsv
      fi
    done
  done
done


# Concatenate all results into a single table
for som_coding_def in union intersection cosmic_only intogen_only; do
  for cancer in breast colorectal lung prostate renal; do
    for germ_context in coding noncoding; do
      for som_context in coding noncoding; do
        if [ $som_context != "coding" ]; then
          res_tsv=results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.pairs.tsv
        else
          res_tsv=results/find_pairs/$cancer.germline_$germ_context.somatic_$som_context.$som_coding_def.pairs.tsv
        fi
        fgrep -v "#" $res_tsv | sed '/^$/d' \
        | awk -v cancer=$cancer -v gc=$germ_context -v sc=$som_context -v OFS="\t" \
          '{ print cancer, $1, gc, $2, sc, $3, $4 }'
      done
    done
  done | fgrep -v "germline_gene" \
  | sort -Vk1,1 -k2,2V -k4,4V -k3,3V -k5,5V \
  | cat \
    <( echo -e "#cancer\tgermline_gene\tgermline_context\tsomatic_gene\tsomatic_context\tcriteria\ttier" ) \
    - \
  > results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.$som_coding_def.tsv
done


# # Print summary of results for tables in slides
# res_tsv=results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.union.tsv
# for germ_context in coding noncoding; do
#   for som_context in coding noncoding; do
#     echo -e "\n\n\n==========\nGermline: $germ_context\nSomatic: $som_context\n=========="
#     awk -v gc=$germ_context -v sc=$som_context \
#       '{ if ($3==gc && $5==sc) print }' $res_tsv | wc -l \
#     | awk '{ print "Total pairs:", $1, "pairs\n" }'
#     for tier in 1 2 3; do
#       awk -v gc=$germ_context -v sc=$som_context -v tr=$tier \
#         '{ if ($3==gc && $5==sc && $NF==tr) print }' $res_tsv | wc -l \
#       | awk -v tr=$tier '{ print "\n* Tier", tr, ":", $1, "pairs" }'
#       if [ $tier == 1 ]; then
#         awk -v gc=$germ_context -v sc=$som_context -v tr=$tier -v OFS="\t" \
#           '{ if ($3==gc && $5==sc && $NF==tr) print $1, $2 }' $res_tsv \
#         | sort -Vk1,1 -k2,2V
#       else
#         awk -v gc=$germ_context -v sc=$som_context -v tr=$tier -v OFS="\t" \
#           '{ if ($3==gc && $5==sc && $NF==tr) print $1, $2, $4 }' $res_tsv \
#         | sort -Vk1,1 -k2,2V -k3,3V
#       fi
#     done
#   done
# done


# Run find_pairs.py for all combination of cancers and negative control phenotypes
while read cancer nc_pheno; do
  for germ_context in coding noncoding; do
    for som_context in coding noncoding; do
      echo -e "$cancer\t$nc_pheno\t$germ_context\t$som_context"
      $CODEDIR/analysis/find_pairs.py \
        --germline gene_lists/germline_$germ_context/$nc_pheno.germline.$germ_context.genes.list \
        --somatic gene_lists/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
        --cellchat-db other_data/cellchat_db_formatted.csv \
        --ppi-db other_data/ebi_intact.all_interactions.tsv \
        --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
        --out-tsv results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.pairs.tsv

      if [ $som_context != "coding" ]; then
        # Only one source of somatic noncoding hits
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$nc_pheno.germline.$germ_context.genes.list \
          --somatic gene_lists/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.pairs.tsv
      else
        # For somatic coding hits, we are evaluating four possible definitions, as follows:
        # 1. Union of intOGen + COSMIC (this is the default)
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$nc_pheno.germline.$germ_context.genes.list \
          --somatic gene_lists/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.union.pairs.tsv
        # 2. Intersection of intOGen + COSMIC
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$nc_pheno.germline.$germ_context.genes.list \
          --somatic other_data/cosmic_intogen_intersection/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.intersection.pairs.tsv
        # 3. COSMIC only
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$nc_pheno.germline.$germ_context.genes.list \
          --somatic other_data/cosmic_genes/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.cosmic_only.pairs.tsv
        # 4. intOGen only
        $CODEDIR/analysis/find_pairs.py \
          --germline gene_lists/germline_$germ_context/$nc_pheno.germline.$germ_context.genes.list \
          --somatic other_data/intogen_genes/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
          --cellchat-db other_data/cellchat_db_formatted.csv \
          --ppi-db other_data/ebi_intact.all_interactions.tsv \
          --protein-complexes other_data/ebi_complex_portal.all_complexes.tsv \
          --out-tsv results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.intogen_only.pairs.tsv
      fi
    done
  done
done < other_data/negative_control_phenotype_pairs.tsv


# Concatenate all results into a single table
for som_coding_def in union intersection cosmic_only intogen_only; do
  while read cancer nc_pheno; do
    for germ_context in coding noncoding; do
      for som_context in coding noncoding; do
        if [ $som_context != "coding" ]; then
          res_tsv=results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.pairs.tsv
        else
          res_tsv=results/find_pairs/${cancer}_${nc_pheno}.germline_$germ_context.somatic_$som_context.$som_coding_def.pairs.tsv
        fi
        fgrep -v "#" $res_tsv | grep -v '^germline_gene' | sed '/^$/d' \
        | awk -v cancer=${cancer}_${nc_pheno} -v gc=$germ_context -v sc=$som_context -v OFS="\t" \
          '{ print cancer, $1, gc, $2, sc, $3, $4 }'
      done
    done | fgrep -v "germline_gene"
  done < other_data/negative_control_phenotype_pairs.tsv \
  | sort -Vk1,1 -k2,2V -k4,4V -k3,3V -k5,5V \
  | cat \
    <( echo -e "#cancer\tgermline_gene\tgermline_context\tsomatic_gene\tsomatic_context\tcriteria\ttier" ) \
    - \
  > results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.$som_coding_def.negative_controls.tsv
done


# Merge primary results and negative control results
for som_coding_def in union intersection cosmic_only intogen_only; do
  cat \
    results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.$som_coding_def.tsv \
    <( sed '1d' results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.$som_coding_def.negative_controls.tsv ) \
  > results/VALab_germline_somatic_2024.v2.gene_pairs.annotated.$som_coding_def.primary_plus_neg_control_merged.tsv
done

