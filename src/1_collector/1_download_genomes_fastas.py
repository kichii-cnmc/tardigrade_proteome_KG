# used to download the initial genome FASTA files from NCBI and UniProt
import os
import argparse
import requests
from tqdm import tqdm
import yaml

def extract_genome_targets(input_yaml):
    '''Extracts the names, sources, and IDs of the genome targets from the config YAML file.'''
    if not os.path.exists(input_yaml):
        raise FileNotFoundError(f"Config file {input_yaml} not found")
    
    try:
        with open(input_yaml, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML file {input_yaml}: {e}")
    
    if config is None:
        raise ValueError(f"Config file {input_yaml} is empty or invalid")
    
    genome_targets = config.get('genome_targets', [])
    return genome_targets

def pull_uniprot_fasta(tax_id, output_fasta):
    '''Downloads the entire proteome of the given taxonomic ID from UniProt as a single FASTA file.'''
    url = f"https://www.uniprot.org/uniprot/?query=organism:{tax_id}&format=fasta"
    response = requests.get(url, stream=True)
    if response.status_code != 200:
        raise Exception(f"Failed to download FASTA from UniProt for tax ID {tax_id}. Status code: {response.status_code}")
    with open(output_fasta, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

def pull_ncbi_fasta(assembly_id, output_fasta):
    pass

if __name__ == "__main__":
    # used for testing purposes, should not be run directly
    print("This script is not intended to be run directly and is only for testing. Please use the main pipeline script.")

    # testing
    genome_targets = extract_genome_targets("config.yaml")
    genome_targets = genome_targets[:1]
    for target in genome_targets:
        organism = target['organism']
        source = target['source']
        id_ = target['id']
        output_fasta = f"data/{organism.replace(' ', '_')}_{source}_test.fasta"
        if source == "UniProt":
            pull_uniprot_fasta(id_, output_fasta)
        elif source == "NCBI":
            pull_ncbi_fasta(id_, output_fasta)
        else:
            print(f"Unknown source {source} for organism {organism}")

    print("Testing completed.")