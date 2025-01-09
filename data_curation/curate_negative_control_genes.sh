# Germline:Somatic Convergence Project
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate germline genes for negative control phenotypes for permutation testing


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence


# Filter GeneBass associations
for pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
  $CODEDIR/data_curation/filter_genebass.R \
    other_data/genebass/$pheno.lof.genebass.csv \
    other_data/genebass/$pheno.mislc.genebass.csv \
    gene_lists/germline_coding/$pheno.germline.coding.genes.list
done


# Filter germline loci from GWAS catalog
for pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
  $CODEDIR/data_curation/filter_gwas_loci.R \
    other_data/gwas_catalog/$pheno.gwas_catalog.12_11_24.unfiltered.tsv \
    other_data/gwas_catalog/$pheno.gwas_catalog.12_11_24.filtered.tsv
done


# Annotate GWAS catalog for gene overlap
for pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
  echo $pheno
  $CODEDIR/data_curation/annotate_gwas_catalog.py \
    --tsv-in other_data/gwas_catalog/$pheno.gwas_catalog.12_11_24.filtered.tsv \
    --gtf ~/Desktop/Collins/VanAllen/germline_somatic_convergence/data/gencode/gencode.v47.annotation.gtf.gz \
    --eligible-genes other_data/gencode.v47.autosomal.protein_coding.genes.list \
    --tsv-out other_data/gwas_catalog/$pheno.gwas_catalog.12_11_24.filtered.annotated.tsv
done


# Divide GWAS catalog into coding & noncoding subsets
for pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
  gwas_tsv=other_data/gwas_catalog/$pheno.gwas_catalog.12_11_24.filtered.annotated.tsv
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
done


# Filter all gene lists to protein-coding symbols defined in Gencode
for pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
  for csq in coding noncoding; do
    fgrep -xf \
      other_data/gencode.v47.autosomal.protein_coding.genes.list \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list \
    > gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list2
    mv \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list2 \
      gene_lists/germline_${csq}/$pheno.germline.$csq.genes.list
  done
done


# Get table of counts
for pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
  cat gene_lists/germline_coding/$pheno.germline.coding.genes.list | wc -l
  cat gene_lists/germline_noncoding/$pheno.germline.noncoding.genes.list | wc -l
  cat \
    gene_lists/germline_coding/$pheno.germline.coding.genes.list \
    gene_lists/germline_noncoding/$pheno.germline.noncoding.genes.list \
  | sort -V | uniq | wc -l
done | paste - - -

