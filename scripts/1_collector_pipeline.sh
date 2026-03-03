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
LOG_FILE="logs/download_metadata.json"
if [ -f "$LOG_FILE" ]; then
    rm "$LOG_FILE"
    echo "Cleared existing download log file."
fi

# clear 3_organized directory
# rm -rf data/3_organized/*.tsv
# echo "Cleared existing organized data files."
# echo

# announce if test mode is enabled
if [ "$TEST_MODE" = true ]; then
    echo "Test mode enabled: The genome download and processing will run on a smaller subset of data for testing purposes."
    echo
fi

# run the genome download script
bash src/1_collector/1_download_genomes.sh
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

# if in test mode, make the merged files smaller and remove original merged files to save space
if [ "$TEST_MODE" = true ]; then
    echo "Test mode enabled: Reducing merged files to first 100 entries for testing."
    if [[ -f "data/2_merged/merged_RV.tsv" ]]; then
        head -n 101 data/2_merged/merged_RV.tsv > data/2_merged/temp_RV.tsv
        mv data/2_merged/temp_RV.tsv data/2_merged/merged_RV.tsv
        rm -f data/2_merged/temp_RV.tsv
    fi
    if [[ -f "data/2_merged/merged_HE.tsv" ]]; then
        head -n 101 data/2_merged/merged_HE.tsv > data/2_merged/temp_HE.tsv
        mv data/2_merged/temp_HE.tsv data/2_merged/merged_HE.tsv
        rm -f data/2_merged/temp_HE.tsv
    fi
fi

# run the protein info pull script on each file
python3 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/merged_RV.tsv data/3_organized/RV.tsv
python3 src/1_collector/3_pull_uniprot_protein_info.py data/2_merged/merged_HE.tsv data/3_organized/HE.tsv

# collect STRING PPI data
# delete existing string PPI files
rm -f data/1_raw/string_ppi/rv_ppi.txt
rm -f data/1_raw/string_ppi/he_ppi.txt
bash src/1_collector/4_download_string_ppi.sh

# process STRING PPI data
# delete existing processed ppi files if they exist to avoid confusion
python3 src/1_collector/5_string_ppi_processor.py data/1_raw/string_ppi/rv_ppi.txt data/3_organized/RV_STRING_ppi.tsv --upid_folder data/3_organized
python3 src/1_collector/5_string_ppi_processor.py data/1_raw/string_ppi/he_ppi.txt data/3_organized/HE_STRING_ppi.tsv --upid_folder data/3_organized
echo "Protein information retrieval completed."

# collect DeepGO annotations
mkdir -p data/1_raw/deepgo_annotations
python3 src/1_collector/6_collect_deepgo.py data/3_organized/ --output_file data/1_raw/deepgo_annotations/deepgo_annotations.tsv $(if [ "$TEST_MODE" = true ]; then echo "--test_mode"; fi)
echo "DeepGO annotation collection completed."

# process DeepGO annotations to create edge list for graph building
python3 src/1_collector/7_process_deepgo.py data/1_raw/deepgo_annotations/deepgo_annotations.tsv data/3_organized/DeepGO.tsv --upid_folder data/3_organized
echo "DeepGO annotation processing completed."

# end of script
echo "Collector pipeline completed."
