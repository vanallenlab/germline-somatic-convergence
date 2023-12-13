# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate all gene lists for germline:somatic pilot analyses

# Note that somatic noncoding gene lists must be curated manually from supplemental excel files


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-exploration-2023


# Curate COSMIC genes
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


# Annotate GWAS catalog for gene overlap
for cancer in breast colorectal lung prostate renal new_pca_only; do
  $CODEDIR/data_curation/annotate_gwas_catalog.py \
    --csv-in gene_lists/germline_non_organized/${cancer}_gwas_catalog_sig_filtered.csv \
    --gtf ~/scratch/gencode.v44.annotation.gtf.gz \
    --tsv-out gene_lists/germline_non_organized/${cancer}_gwas_catalog_sig_filtered.annotated.tsv
done


# Divide GWAS catalog into coding & noncoding subsets
for cancer in breast colorectal lung prostate renal; do
  idx=$( head -n1 gene_lists/germline_non_organized/${cancer}_gwas_catalog_sig_filtered.annotated.tsv \
         | sed 's/\t/\n/g' | awk '{ if ($1=="mappedGenes") print NR }' )

  # Combine coding genes with COSMIC genes
  awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' \
    gene_lists/germline_non_organized/${cancer}_gwas_catalog_sig_filtered.annotated.tsv \
  | cat - gene_lists/germline_coding/$cancer.germline.coding.genes.list \
  | sort -V | uniq \
  > gene_lists/germline_coding/$cancer.germline.coding.genes.list2
  mv \
    gene_lists/germline_coding/$cancer.germline.coding.genes.list2 \
    gene_lists/germline_coding/$cancer.germline.coding.genes.list

  # Extract noncoding genes
  awk -v FS="\t" -v idx=$idx '{ if ($NF!="coding_exon") print $idx }' \
    gene_lists/germline_non_organized/${cancer}_gwas_catalog_sig_filtered.annotated.tsv \
  | fgrep -v mappedGenes | sort -V | uniq \
  > gene_lists/germline_noncoding/$cancer.germline.noncoding.genes.list
done
# Add new prostate GWAS coding results
idx=$( head -n1 gene_lists/germline_non_organized/new_pca_only_gwas_catalog_sig_filtered.annotated.tsv \
       | sed 's/\t/\n/g' | awk '{ if ($1=="Mapped_Gene") print NR }' )
awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' \
  gene_lists/germline_non_organized/new_pca_only_gwas_catalog_sig_filtered.annotated.tsv \
| cat - gene_lists/germline_coding/prostate.germline.coding.genes.list \
| sort -V | uniq \
> gene_lists/germline_coding/prostate.germline.coding.genes.list2
mv \
  gene_lists/germline_coding/prostate.germline.coding.genes.list2 \
  gene_lists/germline_coding/prostate.germline.coding.genes.list
# Add new prostate GWAS noncoding results
awk -v FS="\t" -v idx=$idx '{ if ($NF!="coding_exon") print $idx }' \
  gene_lists/germline_non_organized/new_pca_only_gwas_catalog_sig_filtered.annotated.tsv \
| fgrep -v Mapped_Gene \
| cat - gene_lists/germline_noncoding/prostate.germline.noncoding.genes.list \
| sort -V | uniq \
> gene_lists/germline_noncoding/prostate.germline.noncoding.genes.list2
mv \
  gene_lists/germline_noncoding/prostate.germline.noncoding.genes.list2 \
  gene_lists/germline_noncoding/prostate.germline.noncoding.genes.list



# Filter all gene lists to protein-coding symbols defined in Gencode
bedtools intersect -wa -u \
  -a ~/scratch/gencode.v44.annotation.gtf.gz \
  -b <( awk -v OFS="\t" '{ if ($1 != "chrX" && $1 != "chrY") print $1, 1, $2 }' ~/scratch/hg38.genome ) \
| awk -v FS="\t" '{ if ($3=="gene") print $9 }' | fgrep -w protein_coding \
| sed 's/;/\n/g' | fgrep -w "gene_name" | awk '{ print $NF }' \
| tr -d "\"" | sort -V | uniq \
> other_data/gencode.v44.autosomal.protein_coding.genes.list
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
