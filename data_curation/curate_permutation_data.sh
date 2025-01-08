# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
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
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence
export WRKDIR=`mktemp -d`
export OUTDIR=other_data/permutation_weights
export GTF=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/data/gencode/gencode.v47.annotation.gtf.gz

# Reset output directory
if [ -e $OUTDIR ]; then
  rm -rf $OUTDIR
fi
mkdir $OUTDIR


######################
# 1. Uniform weights #
######################
awk -v OFS="\t" '{ print $1, "1" }' \
  other_data/gencode.v47.autosomal.protein_coding.genes.list \
| sort -Vk1,1 | cat <( echo -e "#gene\tweight" ) - \
> $OUTDIR/gene_weights.uniform.tsv


################################
# 2. Noncoding somatic weights #
################################
# Subset GTF to gene bodies of autosomal protein-coding genes
zcat $GTF \
| fgrep -wf other_data/gencode.v47.autosomal.protein_coding.genes.list \
| awk -v OFS="\t" '{ if ($3=="gene" && $1 !~ /chrX|chrY|chrM/) print $1, $4, $5, $14 }' \
| tr -d '";' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V | bgzip -c \
> $WRKDIR/gencode.v47.autosomal.protein_coding.genes.bed.gz
tabix -p bed -f $WRKDIR/gencode.v47.autosomal.protein_coding.genes.bed.gz

# Divide all intergenic regions in half based on their nearest upstream/downstream genes
awk -v OFS="\t" '{ print $1, "1", $2 }' other_data/hg38.genome | grep -ve 'X\|Y' \
| bedtools subtract -a - -b $WRKDIR/gencode.v47.autosomal.protein_coding.genes.bed.gz \
| bgzip -c > $WRKDIR/intergenic_blocks.bed.gz
bedtools closest \
  -a $WRKDIR/intergenic_blocks.bed.gz \
  -b $WRKDIR/gencode.v47.autosomal.protein_coding.genes.bed.gz \
  -D ref -id \
| awk -v OFS="\t" '{ if ($(NF-1) != ".") print $1, $2, int(($2+$3)/2), $(NF-1) }' \
| bgzip -c > $WRKDIR/intergenic_blocks.upstream.bed.gz
bedtools closest \
  -a $WRKDIR/intergenic_blocks.bed.gz \
  -b $WRKDIR/gencode.v47.autosomal.protein_coding.genes.bed.gz \
  -D ref -iu \
| awk -v OFS="\t" '{ if ($(NF-1) != ".") print $1, int(($2+$3)/2), $3, $(NF-1) }' \
| bgzip -c > $WRKDIR/intergenic_blocks.downstream.bed.gz
zcat \
  $WRKDIR/gencode.v47.autosomal.protein_coding.genes.bed.gz \
  $WRKDIR/intergenic_blocks.upstream.bed.gz \
  $WRKDIR/intergenic_blocks.downstream.bed.gz \
| fgrep -wf other_data/gencode.v47.autosomal.protein_coding.genes.list \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V | bgzip -c \
> $WRKDIR/all_gene_territories.bed.gz
$CODEDIR/data_curation/collect_gene_territory_weights.py \
  --bed-in $WRKDIR/all_gene_territories.bed.gz \
  --genes-list other_data/gencode.v47.autosomal.protein_coding.genes.list \
  --tsv-out $OUTDIR/gene_weights.genome_territory.tsv


##############################################
# 3. Prepare 1000G SNP data for GWAS weights #
##############################################
# Note: for computational expediency, primary data filtering was performed on
# the MGB ERIS cluster. See preprocess_1kg_snps.sh for details

# Copy prepocessed SNP data from ERIStwo to local working directory
scp \
  rlc47@eristwo.partners.org:/data/gusev/USERS/rlc47/convergence/sudmant/sudmant.filtered.*.pruned.sites.bed.gz \
  $WRKDIR/


#################################################
# 4. Prepare GWAS catalog data for GWAS weights #
#################################################
# Filter entire GWAS catalog
$CODEDIR/data_curation/filter_entire_gwas_catalog.R \
  other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.tsv.gz \
  other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.filtered.tsv

# Annotate filtered GWAS catalog for genic overlap
$CODEDIR/data_curation/annotate_gwas_catalog.py \
  --tsv-in other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.filtered.tsv \
  --gtf ~/Desktop/Collins/VanAllen/germline_somatic_convergence/data/gencode/gencode.v47.annotation.gtf.gz \
  --tsv-out other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.filtered.annotated.tsv
gzip -f other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.filtered.annotated.tsv

# Compute gene weights for coding and noncoding variants
$CODEDIR/data_curation/compute_gwas_catalog_weights.R \
  other_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.filtered.annotated.tsv.gz \
  other_data/gencode.v47.autosomal.protein_coding.genes.list \
  $OUTDIR/gene_weights.gwas_catalog


##########################
# 5. Coding GWAS weights #
##########################
# Count number of overlaps between 1000G SNPs and exons from genes of interest
# Ensures all genes have at least one count so no genes have non-zero weights
bedtools intersect -wa \
  -a <( zcat $WRKDIR/gencode.v47lift37.autosomal.protein_coding.exons.bed.gz ) \
  -b $WRKDIR/sudmant.filtered.coding.pruned.sites.bed.gz \
| cut -f4 | sed 's/,/\n/g' \
| cat - other_data/gencode.v47.autosomal.protein_coding.genes.list \
| fgrep -xf other_data/gencode.v47.autosomal.protein_coding.genes.list \
| sort -Vk1,1 | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' \
| cat <( echo -e "#gene\tcoding_SNPs" ) - \
> $OUTDIR/gene_weights.sudmant_coding_snps.tsv

# Composite weight of 1000G SNPs and GWAS catalog hits
$CODEDIR/data_curation/blend_perm_weights.R \
  $OUTDIR/gene_weights.gwas_catalog_coding.tsv \
  $OUTDIR/gene_weights.sudmant_coding_snps.tsv \
  $OUTDIR/gene_weights.composite_germline_coding.tsv


###########################
# 6. Other coding weights #
###########################
# Curate gene-specific nonsynonymous mutation rates from gnomAD
wget \
  -O $WRKDIR/gnomad_metrics.tsv.gz \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
$CODEDIR/data_curation/curate_coding_mutation_rates.R \
  $WRKDIR/gnomad_metrics.tsv.gz \
  other_data/gencode.v47.autosomal.protein_coding.genes.list \
  $OUTDIR/gene_weights.coding_nonsynonymous.tsv


#############################
# 7. Noncoding GWAS weights #
#############################
# Note: relies on data downloaded & processed in several sections above

# Subset GTF to gene bodies of autosomal protein-coding genes
# Use hg19 liftover of Gencode v47 (for 1000G/hg19 compatability)
wget -O - \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.annotation.gtf.gz
| gunzip -c \
| fgrep -wf other_data/gencode.v47.autosomal.protein_coding.genes.list \
| fgrep -w "gene_type \"protein_coding\"" \
| awk -v OFS="\t" '{ if ($3=="gene" && $1 !~ /chrX|chrY|chrM/) print $1, $4, $5, $14 }' \
| tr -d '";' | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V  | bgzip -c \
> $WRKDIR/gencode.v47lift37.autosomal.protein_coding.genes.bed.gz
tabix -p bed -f $WRKDIR/gencode.v47lift37.autosomal.protein_coding.genes.bed.gz

# Find gene nearest to each 1000 Genomes SNP
bedtools closest \
  -a $WRKDIR/sudmant.filtered.noncoding.pruned.sites.bed.gz \
  -b $WRKDIR/gencode.v47lift37.autosomal.protein_coding.genes.bed.gz \
| awk '{ print $NF }' | sed 's/,/\n/g' \
| cat - other_data/gencode.v47.autosomal.protein_coding.genes.list \
| fgrep -xf other_data/gencode.v47.autosomal.protein_coding.genes.list \
| sort -Vk1,1 | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' \
| cat <( echo -e "#gene\tSNPs_nearby" ) - \
> $OUTDIR/gene_weights.sudmant_noncoding_snps.tsv

# Composite weight of 1000G SNPs and GWAS catalog hits
$CODEDIR/data_curation/blend_perm_weights.R \
  $OUTDIR/gene_weights.gwas_catalog_noncoding.tsv \
  $OUTDIR/gene_weights.sudmant_noncoding_snps.tsv \
  $OUTDIR/gene_weights.composite_germline_noncoding.tsv


#####################################
# 8. Curate GTEx expression weights #
#####################################
$CODEDIR/data_curation/curate_gtex_expression.R


###################################################
# 9. Count of genes to be sampled for each strata #
###################################################
# Re-curate COSMIC-only and GeneBass-only gene lists to identify how many *new* 
# germline genes are implicated by coding GWAS hits
mkdir $WRKDIR/cosmic_tmp/
$CODEDIR/data_curation/curate_cosmic_genes.R $WRKDIR/cosmic_tmp/
mkdir $WRKDIR/genebass_tmp
mkdir $WRKDIR/genebass_tmp/germline_coding
for pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
  $CODEDIR/data_curation/filter_genebass.R \
    other_data/genebass/$pheno.lof.genebass.csv \
    other_data/genebass/$pheno.mislc.genebass.csv \
    $WRKDIR/genebass_tmp/germline_coding/$pheno.germline.coding.genes.list
done

# One set of counts for each definition of somatic coding variants
for som_coding_def in union intersection cosmic_only intogen_only; do

  # Loop over all strata and count genes
  for cancer in breast colorectal lung prostate renal; do

    # COSMIC germline coding
    cat $WRKDIR/cosmic_tmp/germline_coding/$cancer.germline.coding.genes.list | wc -l \
    | paste <( echo -e "$cancer\tgermline\tcoding_cosmic" ) -

    # Other germline coding contributed by GWAS (after excluding COSMIC)
    idx=$( head -n1 other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.annotated.tsv \
           | sed 's/\t/\n/g' | awk '{ if ($1=="MAPPED_GENE") print NR }' )
    awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' \
      other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.annotated.tsv \
    | sed 's/, /\n/g' | sed 's/ - /\n/g' | sort | uniq \
    | fgrep -xvf $WRKDIR/cosmic_tmp/germline_coding/$cancer.germline.coding.genes.list \
    | wc -l | paste <( echo -e "$cancer\tgermline\tcoding_gwas" ) -

    # Add germline coding data for negative control phenotypes
    for nc_pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
      # GeneBass coding
      cat $WRKDIR/genebass_tmp/germline_coding/$nc_pheno.germline.coding.genes.list | wc -l \
      | paste <( echo -e "${cancer}_${nc_pheno}\tgermline\tcoding_cosmic" ) -

      # Other germline coding contributed by GWAS (after excluding GeneBass)
      idx=$( head -n1 other_data/gwas_catalog/$nc_pheno.gwas_catalog.12_11_24.filtered.annotated.tsv \
             | sed 's/\t/\n/g' | awk '{ if ($1=="MAPPED_GENE") print NR }' )
      awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' \
        other_data/gwas_catalog/$nc_pheno.gwas_catalog.12_11_24.filtered.annotated.tsv \
      | sed 's/, /\n/g' | sed 's/ - /\n/g' | sort | uniq \
      | fgrep -xvf $WRKDIR/genebass_tmp/germline_coding/$nc_pheno.germline.coding.genes.list \
      | wc -l | paste <( echo -e "${cancer}_${nc_pheno}\tgermline\tcoding_gwas" ) -
    done

    # Germline & somatic noncoding
    for origin in germline somatic; do
      cat gene_lists/${origin}_noncoding/$cancer.$origin.noncoding.genes.list | wc -l \
      | paste <( echo -e "$cancer\t$origin\tnoncoding" ) -
      for nc_pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
        if [ $origin == "germline" ]; then
          cat gene_lists/${origin}_noncoding/$nc_pheno.$origin.noncoding.genes.list | wc -l \
          | paste <( echo -e "${cancer}_${nc_pheno}\t$origin\tnoncoding" ) -
        else
          cat gene_lists/${origin}_noncoding/$cancer.$origin.noncoding.genes.list | wc -l \
          | paste <( echo -e "${cancer}_${nc_pheno}\t$origin\tnoncoding" ) -
        fi
      done
    done

    # Somatic coding depends on definition being used
    case $som_coding_def in
      union)
        som_coding_genes=gene_lists/somatic_coding/$cancer.somatic.coding.genes.list
        ;;
      intersection)
        som_coding_genes=other_data/cosmic_intogen_intersection/somatic_coding/$cancer.somatic.coding.genes.list
        ;;
      cosmic_only)
        som_coding_genes=other_data/cosmic_genes/somatic_coding/$cancer.somatic.coding.genes.list
        ;;
      intogen_only)
        som_coding_genes=other_data/intogen_genes/somatic_coding/$cancer.somatic.coding.genes.list
        ;;
    esac
    cat $som_coding_genes | wc -l | paste <( echo -e "$cancer\tsomatic\tcoding" ) -
    for nc_pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
      cat $som_coding_genes | wc -l | paste <( echo -e "${cancer}_${nc_pheno}\tsomatic\tcoding" ) -
    done

  done > other_data/permutation_gene_counts.$som_coding_def.tsv

done



####################################################################
# 10. Count of genes per strata per expression quantile per tissue #
####################################################################
n_gex_bins=$( sed '1d' $OUTDIR/gene_weights.expression.cancer_specific.tsv \
              | cut -f2 | sort -nr | sed -n '1p' )
mkdir $WRKDIR/gex_elig_lists_tmp
# One set of counts for each definition of somatic coding variants
for som_coding_def in union intersection cosmic_only intogen_only; do
  # Loop over all strata and count genes
  for cancer in breast colorectal lung prostate renal; do
    # Gather counts per expression quantile
    for q in $( seq 1 $n_gex_bins ); do
      # Define background list of genes in this quantile
      awk -v cancer=$cancer -v q=$q -v FS="\t" \
        '{ if ($2==q && $3==cancer) print $1 }' \
        $OUTDIR/gene_weights.expression.cancer_specific.tsv \
      > $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list
      
      # COSMIC germline coding
      fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
        $WRKDIR/cosmic_tmp/germline_coding/$cancer.germline.coding.genes.list | wc -l \
      | paste <( echo -e "$cancer\t$q\tgermline\tcoding_cosmic" ) -
      for nc_pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
        fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
          $WRKDIR/genebass_tmp/germline_coding/$nc_pheno.germline.coding.genes.list | wc -l \
        | paste <( echo -e "${cancer}_${nc_pheno}\t$q\tgermline\tcoding_cosmic" ) -
      done

      # Other germline coding contributed by GWAS (after excluding COSMIC)
      idx=$( head -n1 other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.annotated.tsv \
             | sed 's/\t/\n/g' | awk '{ if ($1=="MAPPED_GENE") print NR }' )
      awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' \
        other_data/gwas_catalog/$cancer.gwas_catalog.12_05_24.filtered.annotated.tsv \
      | sed 's/, /\n/g' | sed 's/ - /\n/g' | sort | uniq \
      | fgrep -xvf $WRKDIR/cosmic_tmp/germline_coding/$cancer.germline.coding.genes.list \
      | fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
      | wc -l | paste <( echo -e "$cancer\t$q\tgermline\tcoding_gwas" ) -
      for nc_pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
        idx=$( head -n1 other_data/gwas_catalog/$nc_pheno.gwas_catalog.12_11_24.filtered.annotated.tsv \
               | sed 's/\t/\n/g' | awk '{ if ($1=="MAPPED_GENE") print NR }' )
        awk -v FS="\t" -v idx=$idx '{ if ($NF=="coding_exon") print $idx }' \
          other_data/gwas_catalog/$nc_pheno.gwas_catalog.12_11_24.filtered.annotated.tsv \
        | sed 's/, /\n/g' | sed 's/ - /\n/g' | sort | uniq \
        | fgrep -xvf $WRKDIR/genebass_tmp/germline_coding/$nc_pheno.germline.coding.genes.list \
        | fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
        | wc -l | paste <( echo -e "${cancer}_${nc_pheno}\t$q\tgermline\tcoding_gwas" ) -
      done

      # Germline & somatic noncoding
      for origin in germline somatic; do
        fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
          gene_lists/${origin}_noncoding/$cancer.$origin.noncoding.genes.list | wc -l \
        | paste <( echo -e "$cancer\t$q\t$origin\tnoncoding" ) -
        for nc_pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
          if [ $origin == "germline" ]; then
            fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
              gene_lists/${origin}_noncoding/$nc_pheno.$origin.noncoding.genes.list | wc -l \
            | paste <( echo -e "${cancer}_${nc_pheno}\t$q\t$origin\tnoncoding" ) -
          else
            fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
              gene_lists/${origin}_noncoding/$cancer.$origin.noncoding.genes.list | wc -l \
            | paste <( echo -e "${cancer}_${nc_pheno}\t$q\t$origin\tnoncoding" ) -
          fi
        done
      done

      # Somatic coding depends on definition being used
      case $som_coding_def in
        union)
          som_coding_genes=gene_lists/somatic_coding/$cancer.somatic.coding.genes.list
          ;;
        intersection)
          som_coding_genes=other_data/cosmic_intogen_intersection/somatic_coding/$cancer.somatic.coding.genes.list
          ;;
        cosmic_only)
          som_coding_genes=other_data/cosmic_genes/somatic_coding/$cancer.somatic.coding.genes.list
          ;;
        intogen_only)
          som_coding_genes=other_data/intogen_genes/somatic_coding/$cancer.somatic.coding.genes.list
          ;;
      esac
      fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
        $som_coding_genes | wc -l \
      | paste <( echo -e "$cancer\t$q\tsomatic\tcoding" ) -
      for nc_pheno in inguinal_hernia atrial_fibrilation myocardial_infarction; do
        fgrep -xf $WRKDIR/gex_elig_lists_tmp/$cancer.gex_q$q.genes.list \
          $som_coding_genes | wc -l \
        | paste <( echo -e "${cancer}_${nc_pheno}\t$q\tsomatic\tcoding" ) -
      done
    done
  done > other_data/gene_counts_per_expression_quantile.$som_coding_def.tsv
done


################
# 11. Clean up #
################
rm -rf $WRKDIR

