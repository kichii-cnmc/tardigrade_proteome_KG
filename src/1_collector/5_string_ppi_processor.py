# This script is used to process STRING PPI data txt files into TSV files for KG building.

import pandas as pd
import os

def find_uniprot_tsvs(folder_path):
    '''Finds all TSV files in the specified folder with UniProt IDs (as second column) and returns a list of their paths.'''
    tsv_files = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.tsv'):
            file_path = os.path.join(folder_path, filename)
            try:
                df = pd.read_csv(file_path, sep='\t')
                if len(df.columns) > 1 and df.columns[1].startswith('UniProt_ID'):
                    tsv_files.append(file_path)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
    return tsv_files

def isolate_protein_id(protein_string):
    '''Isolates the protein ID from a STRING protein string.'''
    return protein_string.split('.')[1] if '.' in protein_string else protein_string

def process_string_ppi_file(input_file, output_file, upid_folder=None, merged_ids_file=None):
    '''Processes a STRING PPI text file and saves it as a TSV file.'''
    df = pd.read_csv(input_file, sep=' ')

    # isolate UniProt IDs from the STRING protein names
    df['protein1'] = df['protein1'].apply(isolate_protein_id)
    df['protein2'] = df['protein2'].apply(isolate_protein_id)

    # removes scores less than a threshold (current: 400)
    df = df[df['combined_score'] >= 400]

    # divide score by 1000 to normalize between 0 and 1
    df['combined_score'] = df['combined_score'] / 1000.0

    # normalize all scores to be between 0 and 1 by dividing by the maximum score in the dataset
    max_score = df['combined_score'].max()
    if max_score > 0:
        df['combined_score'] = df['combined_score'] / max_score

    # rename columns to match KG building expectations
    df = df.rename(columns={'protein1': 'PPI_source', 'protein2': 'PPI_target', 'combined_score': 'score'})

    # apply any conversion of merged IDs
    if merged_ids_file:
        conversion_dict = build_merged_ids_conversion_dict(merged_ids_file)
        df = convert_merged_ids_in_ppi(df, conversion_dict)

    if upid_folder:
        upid_files = find_uniprot_tsvs(upid_folder)
        if upid_files:
            print(f"Found {len(upid_files)} UniProt ID TSV files for sequence filtering.")
            valid_uniprot_ids = set()
            for upid_file in upid_files:
                # build a set of valid UniProt IDs from the files
                upid_df = pd.read_csv(upid_file, sep='\t')
                if len(upid_df.columns) > 1 and upid_df.columns[1].startswith('UniProt_ID'):
                    valid_uniprot_ids.update(upid_df.iloc[:, 1].dropna().astype(str).tolist())
            print(f"Total unique valid UniProt IDs collected: {len(valid_uniprot_ids)}")
            # filter the PPI DataFrame to keep only rows where both source and target are in the valid UniProt ID set
            df = df[df['PPI_source'].isin(valid_uniprot_ids) & df['PPI_target'].isin(valid_uniprot_ids)]
            print(f"PPI data filtered to {len(df)} interactions after applying UniProt ID filtering.")
        else:
            print(f"No UniProt ID TSV files found in {upid_folder}. Skipping UniProt ID filtering.")

    df.to_csv(output_file, sep='\t', index=False)

def build_merged_ids_conversion_dict(input_file):
    '''Builds a dictionary to convert merged UniProt IDs to their current IDs.'''
    conversion_dict = {}
    merged_df = pd.read_csv(input_file, sep='\t')
    for index, row in merged_df.iterrows():
        current_id = row['UniProt_ID']
        merged_ids = str(row['Other_IDs']).split(';') if pd.notna(row['Other_IDs']) else []
        for merged_id in merged_ids:
            merged_id = merged_id.strip()
            if merged_id:
                conversion_dict[merged_id] = current_id
    return conversion_dict

def convert_merged_ids_in_ppi(df, conversion_dict):
    '''Converts merged UniProt IDs in the PPI DataFrame to their current IDs using the provided conversion dictionary.'''
    df['PPI_source'] = df['PPI_source'].apply(lambda x: conversion_dict.get(x, x))
    df['PPI_target'] = df['PPI_target'].apply(lambda x: conversion_dict.get(x, x))
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process STRING PPI data file into TSV format.')
    parser.add_argument('input_file', type=str, help='Path to the input STRING PPI text file.')
    parser.add_argument('output_file', type=str, help='Path to the output TSV file.')
    parser.add_argument('--upid_folder', type=str, default=None, help='Path to the folder containing UniProt IDs TSV files for ID filtering.')
    parser.add_argument('--merge', type=str, default=None, help='Path to the merged UniProt IDs TSV file.')
    args = parser.parse_args()

    print("Processing STRING PPI data...")
    process_string_ppi_file(args.input_file, args.output_file, args.upid_folder, args.merge)
    print(f"Processed STRING PPI data saved to {args.output_file}.")