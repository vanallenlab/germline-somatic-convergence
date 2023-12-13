# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate all data needed for permutation testing

# Note that this requires several sources of upstream data to already have been curated
# This is also helper code that references hard paths on RLCs local machine
# Thus, it won't generalize to other systems


##################
# 0. Basic Setup #
##################
# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023/
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-exploration-2023
export WRKDIR=`mktemp -d`

# Download Gencode v44 GTF
wget \
  -O $WRKDIR/gencode.v44.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz 


######################
# 1. Uniform weights #
######################
awk -v OFS="\t" '{ print $1, "1" }' \
  other_data/gencode.v44.autosomal.protein_coding.genes.list \
| sort -Vk1,1 | cat <( echo -e "#gene\tweight" ) - \
> other_data/permutation_weights/gene_weights.uniform.tsv


################################
# 2. Noncoding somatic weights #
################################
# Subset GTF to gene bodies of autosomal protein-coding genes
zcat $WRKDIR/gencode.v44.annotation.gtf.gz \
| fgrep -wf other_data/gencode.v44.autosomal.protein_coding.genes.list \
| awk -v OFS="\t" '{ if ($3=="gene" && $1 !~ /chrX|chrY|chrM/) print $1, $4, $5, $14 }' \
| tr -d '";' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V | bgzip -c \
> $WRKDIR/gencode.v44.autosomal.protein_coding.genes.bed.gz
tabix -p bed -f $WRKDIR/gencode.v44.autosomal.protein_coding.genes.bed.gz

# Divide all intergenic regions in half based on their nearest upstream/downstream genes
awk -v OFS="\t" '{ print $1, "1", $2 }' other_data/hg38.genome | grep -ve 'X\|Y' \
| bedtools subtract -a - -b $WRKDIR/gencode.v44.autosomal.protein_coding.genes.bed.gz \
| bgzip -c > $WRKDIR/intergenic_blocks.bed.gz
bedtools closest \
  -a $WRKDIR/intergenic_blocks.bed.gz \
  -b $WRKDIR/gencode.v44.autosomal.protein_coding.genes.bed.gz \
  -D ref -id \
| awk -v OFS="\t" '{ if ($(NF-1) != ".") print $1, $2, int(($2+$3)/2), $(NF-1) }' \
| bgzip -c > $WRKDIR/intergenic_blocks.upstream.bed.gz
bedtools closest \
  -a $WRKDIR/intergenic_blocks.bed.gz \
  -b $WRKDIR/gencode.v44.autosomal.protein_coding.genes.bed.gz \
  -D ref -iu \
| awk -v OFS="\t" '{ if ($(NF-1) != ".") print $1, int(($2+$3)/2), $3, $(NF-1) }' \
| bgzip -c > $WRKDIR/intergenic_blocks.downstream.bed.gz
zcat \
  $WRKDIR/gencode.v44.autosomal.protein_coding.genes.bed.gz \
  $WRKDIR/intergenic_blocks.upstream.bed.gz \
  $WRKDIR/intergenic_blocks.downstream.bed.gz \
| fgrep -wf other_data/gencode.v44.autosomal.protein_coding.genes.list \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V | bgzip -c \
> $WRKDIR/all_gene_territories.bed.gz
$CODEDIR/data_curation/collect_gene_territory_weights.py \
  --bed-in $WRKDIR/all_gene_territories.bed.gz \
  --genes-list other_data/gencode.v44.autosomal.protein_coding.genes.list \
  --tsv-out other_data/permutation_weights/gene_weights.genome_territory.tsv


##########################
# 3. Coding GWAS weights #
##########################
# Extract CDS intervals from hg19 liftover of Gencode v44 (for HRC compatability)
wget \
  -O $WRKDIR/gencode.v44lift37.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz
zcat $WRKDIR/gencode.v44lift37.annotation.gtf.gz \
| fgrep -wf other_data/gencode.v44.autosomal.protein_coding.genes.list \
| awk -v OFS="\t" '{ if ($3=="CDS" && $1 !~ /chrX|chrY|chrM/) print $1, $4, $5, $16 }' \
| tr -d '";' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V | sed 's/^chr//g' | uniq \
| bedtools merge -i - -c 4 -o distinct | bgzip -c \
> $WRKDIR/gencode.v44lift37.autosomal.protein_coding.exons.bed.gz
tabix -p bed -f $WRKDIR/gencode.v44lift37.autosomal.protein_coding.exons.bed.gz

# Count number of overlaps between HRC SNPs and exons from genes of interest
wget \
  -O $WRKDIR/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz \
  ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
tabix -p vcf -f $WRKDIR/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
tabix \
  --regions $WRKDIR/gencode.v44lift37.autosomal.protein_coding.exons.bed.gz \
  $WRKDIR/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz \
| awk -v OFS="\t" '{ print $1, $2, $2+1 }' | bgzip -c \
> $WRKDIR/HRC.coding_regions.bed.gz
bedtools intersect -wa \
  -a <( zcat $WRKDIR/gencode.v44lift37.autosomal.protein_coding.exons.bed.gz ) \
  -b $WRKDIR/HRC.coding_regions.bed.gz \
| cut -f4 | sed 's/,/\n/g' | sort -Vk1,1 | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' \
| cat <( echo -e "#gene\tcoding_SNPs" ) - \
> other_data/permutation_weights/gene_weights.hrc_coding_snps.tsv


###########################
# 4. Other coding weights #
###########################
# Curate gene-specific nonsynonymous mutation rates from gnomAD
wget \
  -O $WRKDIR/gnomad_metrics.tsv.gz \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
$CODEDIR/data_curation/curate_coding_mutation_rates.R \
  $WRKDIR/gnomad_metrics.tsv.gz \
  other_data/permutation_weights/gene_weights.coding_nonsynonymous.tsv


#############################
# 5. Noncoding GWAS weights #
#############################
# Note: relies on data downloaded & processed in section 3, above

# Subset GTF to gene bodies of autosomal protein-coding genes
zcat $WRKDIR/gencode.v44lift37.annotation.gtf.gz \
| fgrep -wf other_data/gencode.v44.autosomal.protein_coding.genes.list \
| awk -v OFS="\t" '{ if ($3=="gene" && $1 !~ /chrX|chrY|chrM/) print $1, $4, $5, $14 }' \
| tr -d '";' | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V  | bgzip -c \
> $WRKDIR/gencode.v44lift37.autosomal.protein_coding.genes.bed.gz
tabix -p bed -f $WRKDIR/gencode.v44lift37.autosomal.protein_coding.genes.bed.gz

# Find gene nearest to each noncoding HRC SNP
bedtools intersect -v \
  -a $WRKDIR/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz \
  -b $WRKDIR/gencode.v44lift37.autosomal.protein_coding.exons.bed.gz \
| awk -v OFS="\t" '{ print $1, $2, $2+1 }' \
| bedtools closest -a - -b $WRKDIR/gencode.v44lift37.autosomal.protein_coding.genes.bed.gz \
| awk '{ print $NF }' | sed 's/,/\n/g' | sort -Vk1,1 | uniq -c \
| awk -v OFS="\t" '{ print $2, $1 }' | cat <( echo -e "#gene\tSNPs_nearby" ) - \
> other_data/permutation_weights/gene_weights.hrc_noncoding_snps.tsv


#####################################
# 6. Curate GTEx expression weights #
#####################################
$CODEDIR/data_curation/curate_gtex_expression.R


###################################################################################
# 7. Make composite weights by multiplying Bayesian priors with expression values #
###################################################################################
for suffix in coding_nonsynonymous hrc_coding_snps \
              hrc_noncoding_snps genome_territory; do
  $CODEDIR/data_curation/multiply_weights_by_expression.R \
    other_data/permutation_weights/gene_weights.$suffix.tsv \
    other_data/permutation_weights/gene_weights.expression.cancer_specific.tsv \
    other_data/permutation_weights/gene_weights.$suffix.composite.cancer_specific.tsv
done


###################################################
# 8. Count of genes to be sampled for each strata #
###################################################
# Re-curate COSMIC-only gene lists to identify how many *new* genes are implicated
# by coding GWAS hits
mkdir $WRKDIR/cosmic_tmp/
$CODEDIR/data_curation/curate_cosmic_genes.R $WRKDIR/cosmic_tmp/

# Loop over all strata and count genes
for cancer in breast colorectal lung prostate renal; do
  # COSMIC germline coding
  cat $WRKDIR/cosmic_tmp/germline_coding/$cancer.germline.coding.genes.list | wc -l \
  | paste <( echo -e "$cancer\tgermline\tcoding_cosmic" ) -

  # Other germline coding contributed by GWAS (after excluding COSMIC)
  idx=$( head -n1 gene_lists/germline_non_organized/${cancer}_gwas_catalog_sig_filtered.annotated.tsv \
         | sed 's/\t/\n/g' | awk '{ if ($1=="mappedGenes") print NR }' )
  awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' \
    gene_lists/germline_non_organized/${cancer}_gwas_catalog_sig_filtered.annotated.tsv \
  | fgrep -wvf $WRKDIR/cosmic_tmp/germline_coding/$cancer.germline.coding.genes.list \
  | wc -l | paste <( echo -e "$cancer\tgermline\tcoding_gwas" ) -

  # Germline noncoding
  cat gene_lists/germline_noncoding/$cancer.germline.noncoding.genes.list | wc -l \
  | paste <( echo -e "$cancer\tgermline\tnoncoding" ) -

  # Somatic coding & noncoding
  for context in coding noncoding; do
    cat gene_lists/somatic_$context/$cancer.somatic.$context.genes.list | wc -l \
    | paste <( echo -e "$cancer\tsomatic\t${context}" ) -
  done
done > other_data/permutation_gene_counts.tsv


###############
# 9. Clean up #
###############
rm -rf $WRKDIR

