#!/bin/bash

# This scripts downloads genome FASTA files for target organisms specified in config.yaml 

# if testing flag is set, only download Ramazzottius varieornatus data


# Ramazzottius varieornatus UniProt proteome
curl -o data/rv_uniprot.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=taxonomy_id:947166"
gunzip data/rv_uniprot.fasta.gz

# Ramazzottius varieornatus NCBI genome
curl -L "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_001949185.1/download?include_annotation_type=PROT_FASTA" -o data/rv_ncbi.zip
unzip data/rv_ncbi.zip -d data/rv_ncbi
mv data/rv_ncbi/ncbi_dataset/data/GCA_001949185.1/protein.faa data/rv_ncbi.faa
rm -r data/rv_ncbi data/rv_ncbi.zip

# stop here if testing
if [ "$1" == "test" ]; then
    echo "Test mode: only downloaded Ramazzottius varieornatus data."

