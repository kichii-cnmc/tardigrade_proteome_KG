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

# run the genome download script
bash src/1_collector/1_download_genomes.sh $(if [ "$TEST_MODE" = true ]; then echo "-t"; fi)
echo "Genome download completed."
echo

# run the merge script on each species folders - only if proteomes were downloaded
if [[ -f "data/1_raw/rv_proteome/rv_ncbi_proteome.faa" && -f "data/1_raw/rv_proteome/rv_uniprot_proteome.fasta" ]]; then
    python3.11 src/1_collector/2_cdhit_merge_fastas.py data/1_raw/rv_proteome/rv_ncbi_proteome.faa data/1_raw/rv_proteome/rv_uniprot_proteome.fasta data/2_merged/merged_RV.tsv --identity 1.0
fi

if [[ -f "data/1_raw/he_proteome/he_ncbi_proteome.faa" && -f "data/1_raw/he_proteome/he_uniprot_proteome.fasta" ]]; then
    python3.11 src/1_collector/2_cdhit_merge_fastas.py data/1_raw/he_proteome/he_ncbi_proteome.faa data/1_raw/he_proteome/he_uniprot_proteome.fasta data/2_merged/merged_HE.tsv --identity 1.0
fi
echo "Merging of FASTA files completed."
echo

# run the protein info pull script on each merged file
if [[ -f "data/2_merged/merged_RV.tsv" ]]; then
    python3.11 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/merged_RV.tsv data/3_organized/protein_info_RV.tsv --method batch
fi

if [[ -f "data/2_merged/merged_HE.tsv" ]]; then
    python3.11 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/merged_HE.tsv data/3_organized/protein_info_HE.tsv --method batch
fi