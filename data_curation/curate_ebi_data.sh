# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Helper bash to curate EBI-based gene-gene interactions maps from EBI Complex Portal and EBI IntAct

# Note that EBI database .zip files must be downloaded and extracted in advance


# Set parameters
cd ~/Dropbox\ \(Partners\ HealthCare\)/VALab_germline_somatic_2023
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/germ_som_pilot_2023/germline-somatic-exploration-2023
export ppi_data=~/Downloads/ebi_intact_db
export complex_data=~/Downloads/ebi_complex_portal


# Curate EBI IntAct PPI data
$CODEDIR/data_curation/curate_ebi_complexes.py \
  --out-tsv other_data/ebi_intact.all_interactions.tsv \
  $ppi_data/*.xml


# Curate EBI Complex Portal data
$CODEDIR/data_curation/curate_ebi_complexes.py \
  --out-tsv other_data/ebi_complex_portal.all_complexes.tsv \
  $complex_data/human/*.xml

