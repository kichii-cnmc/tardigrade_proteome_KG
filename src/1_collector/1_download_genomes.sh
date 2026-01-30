#!/bin/bash

# This scripts downloads genome FASTA files for target organisms specified in config.yaml 
echo "Downloading genome FASTA files for target organisms..."

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

# Parse command line arguments for testing flag
TEST_MODE=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -t)
            TEST_MODE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [-t]"
            echo "  -t: Enable test mode (stop early)"
            exit 1
            ;;
    esac
done

# if none, make data/1_raw/rv_proteome directory
mkdir -p data/1_raw/rv_proteome

# Ramazzottius varieornatus UniProt proteome
echo "Downloading Ramazzottius varieornatus UniProt proteome..."
curl -o data/rv_uniprot.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=taxonomy_id:947166"
gunzip data/rv_uniprot.fasta.gz
mv data/rv_uniprot.fasta data/1_raw/rv_proteome/rv_uniprot_proteome.fasta

# Log the UniProt download
log_download "data/1_raw/rv_proteome/rv_uniprot_proteome.fasta" \
             "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=taxonomy_id:947166" \
             "UniProt" \
             "Ramazzottius varieornatus" \
             "947166" \
             "current"

# Ramazzottius varieornatus NCBI genome
echo "Downloading Ramazzottius varieornatus NCBI proteome..."
curl -L "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_001949185.1/download?include_annotation_type=PROT_FASTA" -o data/rv_ncbi.zip
unzip data/rv_ncbi.zip -d data/rv_ncbi
mv data/rv_ncbi/ncbi_dataset/data/GCA_001949185.1/protein.faa data/1_raw/rv_proteome/rv_ncbi_proteome.faa
rm -r data/rv_ncbi data/rv_ncbi.zip

# Log the NCBI download
log_download "data/1_raw/rv_proteome/rv_ncbi_proteome.faa" \
             "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_001949185.1/download?include_annotation_type=PROT_FASTA" \
             "NCBI Datasets" \
             "Ramazzottius varieornatus" \
             "947166" \
             "GCA_001949185.1"

# Stop here if testing flag is set
if [ "$TEST_MODE" = true ]; then
    echo "Test Mode - Stopping Downloads Early"
    exit 0
fi

# Hypsibius exemplaris genome directory
mkdir -p data/1_raw/he_proteome

# Hypsibius exemplaris UniProt proteome
echo "Downloading Hypsibius exemplaris UniProt proteome..."
curl -o data/he_uniprot.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=taxonomy_id:2072580"
gunzip data/he_uniprot.fasta.gz
mv data/he_uniprot.fasta data/1_raw/he_proteome/he_uniprot_proteome.fasta

# Log the UniProt download
log_download "data/1_raw/he_proteome/he_uniprot_proteome.fasta" \
             "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=taxonomy_id:2072580" \
             "UniProt" \
             "Hypsibius exemplaris" \
             "2072580" \
             "current"

# Hypsibius exemplaris NCBI genome
echo "Downloading Hypsibius exemplaris NCBI proteome..."
curl -L "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_002082055.1/download?include_annotation_type=PROT_FASTA" -o data/he_ncbi.zip
unzip data/he_ncbi.zip -d data/he_ncbi
mv data/he_ncbi/ncbi_dataset/data/GCA_002082055.1/protein.faa data/1_raw/he_proteome/he_ncbi_proteome.faa
rm -r data/he_ncbi data/he_ncbi.zip

# Log the NCBI download
log_download "data/1_raw/he_proteome/he_ncbi_proteome.faa" \
             "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_002082055.1/download?include_annotation_type=PROT_FASTA" \
             "NCBI Datasets" \
             "Hypsibius exemplaris" \
             "2072580" \
             "GCA_002082055.1"

echo "Genome downloads completed."