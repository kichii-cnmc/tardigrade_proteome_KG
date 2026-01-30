#!/bin/bash

# This script is used to download the STRING PPI data for Ramazzottius varieornatus (missing for Hypsibius exemplaris)
echo "Downloading STRING PPI data..."

# Function to log download metadata to JSON file
log_download() {
    local file="$1"
    local url="$2"
    local source="$3"
    local organism="$4"
    local taxonomy_id="$5"
    local version="$6"
    local script_name="$(basename "$0")"
    local timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    
    mkdir -p logs
    local log_file="logs/download_metadata.json"
    
    # Create new entry
    local new_entry="{
      \"timestamp\": \"$timestamp\",
      \"script\": \"$script_name\",
      \"source\": \"$source\",
      \"file\": \"$file\",
      \"url\": \"$url\",
      \"organism\": \"$organism\",
      \"taxonomy_id\": \"$taxonomy_id\",
      \"version\": \"$version\"
    }"
    
    # Initialize file if it doesn't exist
    if [ ! -f "$log_file" ]; then
        echo '{"downloads": []}' > "$log_file"
    fi
    
    # Add new entry (simple append - assumes valid JSON structure)
    if grep -q '"downloads": \[\]' "$log_file"; then
        # Empty array - add first entry
        sed -i '' 's/"downloads": \[\]/"downloads": [\
    '"$new_entry"'\
  ]/' "$log_file"
    else
        # Has entries - add to array
        sed -i '' 's/\(.*\)\]/\1,\
    '"$new_entry"'\
  ]/' "$log_file"
    fi
    
    echo "Download logged to $log_file"
}

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
# WIP: Currently not available on STRING database, so skipping this step

echo "STRING PPI data download completed."
echo