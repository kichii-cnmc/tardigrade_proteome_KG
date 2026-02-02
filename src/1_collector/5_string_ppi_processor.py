# This script is used to process STRING PPI data txt files into TSV files for KG building.

import pandas as pd

def isolate_protein_id(protein_string):
    '''Isolates the protein ID from a STRING protein string.'''
    return protein_string.split('.')[1] if '.' in protein_string else protein_string

def process_string_ppi_file(input_file, output_file):
    '''Processes a STRING PPI text file and saves it as a TSV file.'''
    df = pd.read_csv(input_file, sep=' ')

    # isolate UniProt IDs from the STRING protein names
    df['protein1'] = df['protein1'].apply(isolate_protein_id)
    df['protein2'] = df['protein2'].apply(isolate_protein_id)

    # removes scores less than a threshold (current: 400)
    df = df[df['combined_score'] >= 400]

    # divide score by 1000 to normalize between 0 and 1
    df['combined_score'] = df['combined_score'] / 1000.0

    # rename columns to match KG building expectations
    df = df.rename(columns={'protein1': 'PPI_source', 'protein2': 'PPI_target', 'combined_score': 'score'})

    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process STRING PPI data file into TSV format.')
    parser.add_argument('input_file', type=str, help='Path to the input STRING PPI text file.')
    parser.add_argument('output_file', type=str, help='Path to the output TSV file.')
    args = parser.parse_args()

    print("Processing STRING PPI data...")
    process_string_ppi_file(args.input_file, args.output_file)
    print(f"Processed STRING PPI data saved to {args.output_file}.")