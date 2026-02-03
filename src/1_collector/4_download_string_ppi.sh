#!/bin/bash

# This script is used to download the STRING PPI data for Ramazzottius varieornatus (missing for Hypsibius exemplaris)
echo "Downloading STRING PPI data..."

# Source the shared logging utility
source "$(dirname "$0")/download_logger.sh"

# if none make directory for STRING data
mkdir -p data/1_raw/string_ppi

# Download STRING PPI data for Ramazzottius varieornatus (taxonomy ID: 947166)
echo "Downloading Ramazzottius varieornatus PPI data..."
curl -o data/1_raw/string_ppi/rv_ppi.txt.gz "https://stringdb-downloads.org/download/protein.links.v12.0/947166.protein.links.v12.0.txt.gz"
gunzip data/1_raw/string_ppi/rv_ppi.txt.gz

# Log the download
log_download "data/1_raw/string_ppi/rv_ppi.txt" \
             "https://stringdb-downloads.org/download/protein.links.v12.0/947166.protein.links.v12.0.txt.gz" \
             "STRING database" \
             "Ramazzottius varieornatus" \
             "947166" \
             "v12.0"

# Download STRING PPI data for Hypsibius exemplaris (taxonomy ID: 2072580)
echo "Downloading Hypsibius exemplaris PPI data..."
curl -o data/1_raw/string_ppi/he_ppi.txt.gz "https://stringdb-downloads.org/download_proteomes/protein.links.v12.0/STRG0A14BLM.protein.links.v12.0.txt.gz"
gunzip data/1_raw/string_ppi/he_ppi.txt.gz

# Log the download
log_download "data/1_raw/string_ppi/he_ppi.txt" \
             "https://stringdb-downloads.org/download_proteomes/protein.links.v12.0/STRG0A14BLM.protein.links.v12.0.txt.gz" \
             "STRING database" \
             "Hypsibius exemplaris" \
             "2072580" \
             "v12.0"

echo "STRING PPI data download completed."
echo