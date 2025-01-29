# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate germline genes for negative control phenotypes for permutation testing


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence


# Enumerate negative control pairings
TAB="$( printf '\t' )"
cat << EOF > other_data/negative_control_phenotype_pairs.tsv
breast${TAB}inguinal_hernia
colorectal${TAB}hair_loss
lung${TAB}menarche
prostate${TAB}anxiety
renal${TAB}reflux
EOF


# Filter GeneBass associations
while read cancer pheno; do
  $CODEDIR/data_curation/filter_genebass.R \
    other_data/genebass/$pheno.lof.genebass.csv \
    other_data/genebass/$pheno.mislc.genebass.csv \
    gene_lists/germline_coding/$pheno.germline.coding.genes.list
done < other_data/negative_control_phenotype_pairs.tsv


# Filter germline loci from GWAS catalog
while read cancer pheno; do
  $CODEDIR/data_curation/filter_gwas_loci.R \
    other_data/gwas_catalog/$pheno.gwas_catalog.01_17_25.unfiltered.tsv \
    other_data/gwas_catalog/$pheno.gwas_catalog.01_17_25.filtered.tsv
done < other_data/negative_control_phenotype_pairs.tsv


# Annotate GWAS catalog for gene overlap
while read cancer pheno; do
  echo $pheno
  $CODEDIR/data_curation/annotate_gwas_catalog.py \
    --tsv-in other_data/gwas_catalog/$pheno.gwas_catalog.01_17_25.filtered.tsv \
    --gtf ~/Desktop/Collins/VanAllen/germline_somatic_convergence/data/gencode/gencode.v47.annotation.gtf.gz \
    --eligible-genes other_data/gencode.v47.autosomal.protein_coding.genes.list \
    --tsv-out other_data/gwas_catalog/$pheno.gwas_catalog.01_17_25.filtered.annotated.tsv
done < other_data/negative_control_phenotype_pairs.tsv


# Divide GWAS catalog into coding & noncoding subsets
while read cancer pheno; do
  gwas_tsv=other_data/gwas_catalog/$pheno.gwas_catalog.01_17_25.filtered.annotated.tsv
  idx=$( head -n1 $gwas_tsv | sed 's/\t/\n/g' | awk '{ if ($1=="MAPPED_GENE") print NR }' )

  # Combine coding genes with GeneBass genes
  awk -v FS="\t" -v idx=$idx '{ if ($NF=="CDS") print $idx }' $gwas_tsv \
  | sed 's/,/\n/g' \
  | cat - gene_lists/germline_coding/$pheno.germline.coding.genes.list \
  | sort -V | uniq \
  > gene_lists/germline_coding/$pheno.germline.coding.genes.list2
  mv \
    gene_lists/germline_coding/$pheno.germline.coding.genes.list2 \
    gene_lists/germline_coding/$pheno.germline.coding.genes.list

  # Extract noncoding genes
  awk -v FS="\t" -v idx=$idx '{ if ($NF!="CDS") print $idx }' $gwas_tsv \
  | sed 's/, /\n/g' | sed 's/ - /\n/g' | fgrep -v MAPPED_GENE | sort -V | uniq \
  > gene_lists/germline_noncoding/$pheno.germline.noncoding.genes.list
done < other_data/negative_control_phenotype_pairs.tsv


# Filter all gene lists to protein-coding symbols defined in Gencode
while read cancer pheno; do
  for csq in coding noncoding; do
    fgrep -xf \
      other_data/gencode.v47.autosomal.protein_coding.genes.list \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list \
    > gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list2
    mv \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list2 \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list
  done
done < other_data/negative_control_phenotype_pairs.tsv


# Exclude pleiotropic genes matched for cancer type & germline context
while read cancer pheno; do
  for csq in coding noncoding; do
    fgrep -xvf \
      gene_lists/germline_${csq}/$cancer.germline.$csq.genes.list \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list \
    > gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list2
    mv \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list2 \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list
  done
done < other_data/negative_control_phenotype_pairs.tsv


# Get table of counts
while read cancer pheno; do
  cat gene_lists/germline_coding/$pheno.germline.coding.genes.list | wc -l
  cat gene_lists/germline_noncoding/$pheno.germline.noncoding.genes.list | wc -l
  cat \
    gene_lists/germline_coding/$pheno.germline.coding.genes.list \
    gene_lists/germline_noncoding/$pheno.germline.noncoding.genes.list \
  | sort -V | uniq | wc -l
done < other_data/negative_control_phenotype_pairs.tsv | paste - - -

