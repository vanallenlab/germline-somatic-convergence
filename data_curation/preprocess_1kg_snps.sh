#!/usr/bin/env bash

# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate 1000 Genomes SNP data for permutation weights

# Note: intended to be executed on the MGB ERISOne cluster


# Set local parameters
export WRKDIR=/data/gusev/USERS/rlc47/convergence
export KGDIR=/data/gusev/1000G/phase3
cd $WRKDIR


# Set up directory trees as necessary
for SUBDIR in LSF LSF/logs LSF/scripts; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


# Normalize and filter each VCF per autosome
for contig in $( seq 1 22 ); do
  cat <<EOF > $WRKDIR/LSF/scripts/filter_1kg.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools norm -m- \
  $KGDIR/ALL.chr$contig.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
| bcftools view \
  --include 'INFO/AF >= 0.01 & INFO/AF <= 0.99' \
  --types snps -m2 -M2 \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  -Oz -o $WRKDIR/sudmant/sudmant.chr$contig.filtered.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/filter_1kg.$contig.sh
  rm $WRKDIR/LSF/logs/filter_1kg.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -n 2 -J filter_1kg_$contig \
    -o $WRKDIR/LSF/logs/filter_1kg.$contig.log \
    -e $WRKDIR/LSF/logs/filter_1kg.$contig.err \
    $WRKDIR/LSF/scripts/filter_1kg.$contig.sh
done


# Extract CDS intervals from hg19 liftover of Gencode v47 (for 1000G/hg19 compatability)
wget \
  -O $WRKDIR/gencode.v47lift37.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.annotation.gtf.gz
zcat $WRKDIR/gencode.v47lift37.annotation.gtf.gz \
| fgrep -w "transcript_type \"protein_coding\"" \
| awk -v OFS="\t" '{ if ($3=="CDS" && $1 !~ /chrX|chrY|chrM/) print $1, $4, $5, $16 }' \
| tr -d '";' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V | sed 's/^chr//g' | uniq \
| bedtools merge -i - -c 4 -o distinct | bgzip -c \
> $WRKDIR/gencode.v47lift37.autosomal.protein_coding.exons.bed.gz
tabix -p bed -f $WRKDIR/gencode.v47lift37.autosomal.protein_coding.exons.bed.gz


# Split VCF by coding vs. noncoding variants
for contig in $( seq 1 22 ); do
  cat <<EOF > $WRKDIR/LSF/scripts/split_coding.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
  # Coding regions
  bcftools view \
    -T $WRKDIR/gencode.v47lift37.autosomal.protein_coding.exons.bed.gz \
    -Oz -o $WRKDIR/sudmant/sudmant.chr$contig.filtered.coding.vcf.gz \
    $WRKDIR/sudmant/sudmant.chr$contig.filtered.vcf.gz
  # Excluding coding regions
  bcftools view \
    -T ^$WRKDIR/gencode.v47lift37.autosomal.protein_coding.exons.bed.gz \
    -Oz -o $WRKDIR/sudmant/sudmant.chr$contig.filtered.noncoding.vcf.gz \
    $WRKDIR/sudmant/sudmant.chr$contig.filtered.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/split_coding.$contig.sh
  rm $WRKDIR/LSF/logs/split_coding.$contig.*
  bsub \
    -q normal -J split_coding_$contig \
    -o $WRKDIR/LSF/logs/split_coding.$contig.log \
    -e $WRKDIR/LSF/logs/split_coding.$contig.err \
    $WRKDIR/LSF/scripts/split_coding.$contig.sh
done


# LD prune coding and noncoding SNPs, and convert to simple site BED files
for contig in $( seq 1 22 ); do
  cat <<EOF > $WRKDIR/LSF/scripts/ld_prune.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
module load Plink/1.9
cd $WRKDIR
for context in coding noncoding; do
  plink \
    --vcf $WRKDIR/sudmant/sudmant.chr$contig.filtered.\$context.vcf.gz \
    --indep-pairwise 100kb 1 0.5 \
    --threads 4 \
    --memory 6000 \
    --out $WRKDIR/sudmant/sudmant.chr$contig.filtered.\$context.pruned
  bcftools view \
    --no-version \
    --include "ID=@$WRKDIR/sudmant/sudmant.chr$contig.filtered.\$context.pruned.prune.in" \
    $WRKDIR/sudmant/sudmant.chr$contig.filtered.\$context.vcf.gz \
  | bcftools query -f '%CHROM\t%POS\n' \
  | awk -v OFS="\t" '{ print \$1, \$2, \$2+1 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bgzip -c \
  > $WRKDIR/sudmant/sudmant.chr$contig.filtered.\$context.pruned.sites.bed.gz
  tabix -f $WRKDIR/sudmant/sudmant.chr$contig.filtered.\$context.pruned.sites.bed.gz
done
EOF
  chmod a+x $WRKDIR/LSF/scripts/ld_prune.$contig.sh
  rm $WRKDIR/LSF/logs/ld_prune.$contig.*
  bsub \
    -q normal -J ld_prune_$contig \
    -o $WRKDIR/LSF/logs/ld_prune.$contig.log \
    -e $WRKDIR/LSF/logs/ld_prune.$contig.err \
    $WRKDIR/LSF/scripts/ld_prune.$contig.sh
done


# Merge filtered & pruned 1000G variants
for context in coding noncoding; do
  zcat $WRKDIR/sudmant/sudmant.chr*.filtered.$context.pruned.sites.bed.gz \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | bgzip -c \
  > $WRKDIR/sudmant/sudmant.filtered.$context.pruned.sites.bed.gz
  tabix -f $WRKDIR/sudmant/sudmant.filtered.$context.pruned.sites.bed.gz
done

