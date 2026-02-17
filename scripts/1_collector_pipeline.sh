#!/bin/bash

# This script is used to run the entire collector pipeline for downloading and processing genome data
# it also includes a test mode to stop early for testing purposes

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

# clear logging file
LOG_FILE = "logs/download_metadata.json"
if [ -f "$LOG_FILE" ]; then
    rm "$LOG_FILE"
    echo "Cleared existing download log file."
fi

# clear 3_organized directory
rm -rf data/3_organized/*.tsv
echo "Cleared existing organized data files."
echo

# run the genome download script
bash src/1_collector/1_download_genomes.sh $(if [ "$TEST_MODE" = true ]; then echo "-t"; fi)
echo "Genome download completed."
echo

# run the merge script on each species folders - only if proteomes were downloaded
if [[ -f "data/1_raw/rv_proteome/rv_ncbi_proteome.faa" && -f "data/1_raw/rv_proteome/rv_uniprot_proteome.fasta" ]]; then
    python3 src/1_collector/2_cdhit_merge_fastas.py data/1_raw/rv_proteome/rv_ncbi_proteome.faa data/1_raw/rv_proteome/rv_uniprot_proteome.fasta data/2_merged/merged_RV.tsv --identity 1.0
fi

if [[ -f "data/1_raw/he_proteome/he_ncbi_proteome.faa" && -f "data/1_raw/he_proteome/he_uniprot_proteome.fasta" ]]; then
    python3 src/1_collector/2_cdhit_merge_fastas.py data/1_raw/he_proteome/he_ncbi_proteome.faa data/1_raw/he_proteome/he_uniprot_proteome.fasta data/2_merged/merged_HE.tsv --identity 1.0
fi
echo "Merging of FASTA files completed."
echo

# if in test mode, make the merged files smaller for testing
if [ "$TEST_MODE" = true ]; then
    echo "Test mode enabled: Reducing merged files to first 50 entries for testing."
    if [[ -f "data/2_merged/merged_RV.tsv" ]]; then
        head -n 51 data/2_merged/merged_RV.tsv > data/2_merged/test_RV.tsv
    fi
    if [[ -f "data/2_merged/merged_HE.tsv" ]]; then
        head -n 51 data/2_merged/merged_HE.tsv > data/2_merged/test_HE.tsv
    fi
fi

# run the protein info pull script on each test file
if [ "$TEST_MODE" = true ]; then
    echo "Test mode enabled: Pulling protein info for test files."
    if [[ -f "data/2_merged/test_RV.tsv" ]]; then
    python3 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/test_RV.tsv data/3_organized/test_UniProt_RV.tsv
    fi
    if [[ -f "data/2_merged/test_HE.tsv" ]]; then
        python3 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/test_HE.tsv data/3_organized/test_UniProt_HE.tsv
    fi
fi

# run the protein info pull script on each full file if not in test mode
if [ "$TEST_MODE" = false ]; then
    if [[ -f "data/2_merged/merged_RV.tsv" ]]; then
        python3 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/merged_RV.tsv data/3_organized/UniProt_RV.tsv
    fi  
    if [[ -f "data/2_merged/merged_HE.tsv" ]]; then
        python3 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/merged_HE.tsv data/3_organized/UniProt_HE.tsv
    fi
fi

# collect STRING PPI data
bash src/1_collector/4_download_string_ppi.sh

# process STRING PPI data
if [[ -f "data/1_raw/string_ppi/rv_ppi.txt" ]]; then
    python3 src/1_collector/5_string_ppi_processor.py data/1_raw/string_ppi/rv_ppi.txt data/3_organized/STRING_RV_PPI.tsv
fi
if [[ -f "data/1_raw/string_ppi/he_ppi.txt" ]]; then
    python3 src/1_collector/5_string_ppi_processor.py data/1_raw/string_ppi/he_ppi.txt data/3_organized/STRING_HE_PPI.tsv
fi

echo "Protein information retrieval completed."