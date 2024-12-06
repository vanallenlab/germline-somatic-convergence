# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2024 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate all gene lists for germline:somatic pilot analyses

# Note that somatic noncoding gene lists must be curated manually from supplemental excel files


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence


# Curate COSMIC genes
# Starts from COSMIC gene census .tsv downloaded from COSMIC website
$CODEDIR/data_curation/curate_cosmic_genes.R other_data/cosmic_genes/
if [ -e gene_lists/germline_coding ]; then
  rm -rf gene_lists/germline_coding
fi
mv -f other_data/cosmic_genes/germline_coding gene_lists/


# Curate IntOGen genes
$CODEDIR/data_curation/curate_intogen_genes.R


# Union of COSMIC + IntOGen for main list of somatic driver genes
for cancer in breast colorectal lung prostate renal; do
  cat \
    other_data/cosmic_genes/somatic_coding/$cancer.somatic.coding.genes.list \
    other_data/intogen_genes/somatic_coding/$cancer.somatic.coding.genes.list \
  | sort -V | uniq \
  > gene_lists/somatic_coding/$cancer.somatic.coding.genes.list
done


# Intersection of COSMIC + IntOGen for high-PPV somatic driver genes
for dir in \
  other_data/cosmic_intogen_intersection \
  other_data/cosmic_intogen_intersection/somatic_coding; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done
for cancer in breast colorectal lung prostate renal; do
  fgrep -wf \
    other_data/cosmic_genes/somatic_coding/$cancer.somatic.coding.genes.list \
    other_data/intogen_genes/somatic_coding/$cancer.somatic.coding.genes.list \
  > other_data/cosmic_intogen_intersection/somatic_coding/$cancer.somatic.coding.genes.list
done


# Filter germline loci from GWAS catalog
for cancer in breast colorectal lung prostate renal; do
  $CODEDIR/data_curation/filter_gwas_loci.R \
    other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.unfiltered.tsv \
    other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.tsv
done


# Annotate GWAS catalog for gene overlap
for cancer in breast colorectal lung prostate renal; do
  $CODEDIR/data_curation/annotate_gwas_catalog.py \
    --tsv-in other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.tsv \
    --gtf ~/Desktop/Collins/VanAllen/germline_somatic_convergence/data/gencode/gencode.v47.annotation.gtf.gz \
    --tsv-out other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.annotated.tsv
done


# Divide GWAS catalog into coding & noncoding subsets
for cancer in breast colorectal lung prostate renal; do
  gwas_tsv=other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.annotated.tsv
  idx=$( head -n1 $gwas_tsv | sed 's/\t/\n/g' | awk '{ if ($1=="MAPPED_GENE") print NR }' )

  # Combine coding genes with COSMIC genes
  awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' $gwas_tsv \
  | sed 's/, /\n/g' \
  | cat - gene_lists/germline_coding/$cancer.germline.coding.genes.list \
  | sort -V | uniq \
  > gene_lists/germline_coding/$cancer.germline.coding.genes.list2
  mv \
    gene_lists/germline_coding/$cancer.germline.coding.genes.list2 \
    gene_lists/germline_coding/$cancer.germline.coding.genes.list

  # Extract noncoding genes
  awk -v FS="\t" -v idx=$idx '{ if ($NF!="coding_exon") print $idx }' $gwas_tsv \
  | sed 's/, /\n/g' | sed 's/ - /\n/g' | fgrep -v MAPPED_GENE | sort -V | uniq \
  > gene_lists/germline_noncoding/$cancer.germline.noncoding.genes.list
done


# Filter all gene lists to protein-coding symbols defined in Gencode
bedtools intersect -wa -u \
  -a ~/Desktop/Collins/VanAllen/germline_somatic_convergence/data/gencode/gencode.v47.annotation.gtf.gz \
  -b <( awk -v OFS="\t" '{ if ($1 != "chrX" && $1 != "chrY") print $1, 1, $2 }' \
          ~/Desktop/Collins/VanAllen/germline_somatic_convergence/data/gencode/hg38.genome ) \
| awk -v FS="\t" '{ if ($3=="gene") print $9 }' | fgrep -w protein_coding \
| sed 's/;/\n/g' | fgrep -w "gene_name" | awk '{ print $NF }' \
| tr -d "\"" | sort -V | uniq \
> other_data/gencode.v47.autosomal.protein_coding.genes.list
for cancer in breast colorectal lung prostate renal; do
  for context in somatic germline; do
    for csq in coding noncoding; do
      fgrep -wf \
        other_data/gencode.v44.autosomal.protein_coding.genes.list \
        gene_lists/${context}_${csq}/$cancer.$context.$csq.genes.list \
      > gene_lists/${context}_${csq}/$cancer.$context.$csq.genes.list2
      mv \
        gene_lists/${context}_${csq}/$cancer.$context.$csq.genes.list2 \
        gene_lists/${context}_${csq}/$cancer.$context.$csq.genes.list
    done
  done
done


# Get table of counts
for cancer in breast colorectal lung prostate renal; do
  cat gene_lists/germline_coding/$cancer.germline.coding.genes.list | wc -l
  cat gene_lists/germline_noncoding/$cancer.germline.noncoding.genes.list | wc -l
  cat \
    gene_lists/germline_coding/$cancer.germline.coding.genes.list \
    gene_lists/germline_noncoding/$cancer.germline.noncoding.genes.list \
  | sort -V | uniq | wc -l

  cat gene_lists/somatic_coding/$cancer.somatic.coding.genes.list | wc -l
  cat gene_lists/somatic_noncoding/$cancer.somatic.noncoding.genes.list | wc -l
  cat \
    gene_lists/somatic_coding/$cancer.somatic.coding.genes.list \
    gene_lists/somatic_noncoding/$cancer.somatic.noncoding.genes.list \
  | sort -V | uniq | wc -l  

  find gene_lists/ -name "$cancer.*.*.genes.list" \
  | xargs -I {} cat {} | sort -V | uniq | wc -l
done | paste - - - - - - -
