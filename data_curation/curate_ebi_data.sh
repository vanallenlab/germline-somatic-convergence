# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate EBI-based gene-gene interactions maps from EBI Complex Portal and EBI IntAct

# Note that EBI database .zip files must be downloaded and extracted in advance


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VanAllen/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/germline-somatic-convergence
export ppi_data=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/data/intact_ppi_human_dec5_2024
export complex_data=/Users/ryan/Desktop/Collins/VanAllen/germline_somatic_convergence/data/ebi_complex_portal_dec5_2024


# Curate EBI IntAct PPI data
# Downloaded the full human interactome in miXML 3.0 format from:
# https://www.ebi.ac.uk/intact/interactomes
$CODEDIR/data_curation/curate_ebi_ppis.py \
  --eligible-genes other_data/gencode.v47.autosomal.protein_coding.genes.list \
  --out-tsv other_data/ebi_intact.all_interactions.tsv \
  $ppi_data/*.xml


# Curate EBI Complex Portal data
# Downloaded in miSML 3.0 format from:
# https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/psi30/human.zip
$CODEDIR/data_curation/curate_ebi_complexes.py \
  --out-tsv other_data/ebi_complex_portal.all_complexes.tsv \
  $complex_data/*.xml

