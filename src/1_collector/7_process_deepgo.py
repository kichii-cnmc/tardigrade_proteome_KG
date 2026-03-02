# script converts raw deepgo output into a edge list format for graph building, integrating scores.

import os

import pandas as pd
import argparse

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

def process_deepgo_annotations(input_file, output_file, upid_folder=None):
    '''Processes raw DeepGO annotations into an edge list format for graph building.'''
    # add to df
    df = pd.read_csv(input_file, sep='\t')

    # confirm that columns are correct
    expected_columns = ['SwissProt ID', 'Function Type', 'GO Term', 'Score']
    if not all(col in df.columns for col in expected_columns):
        raise ValueError(f"Input file must contain the following columns: {expected_columns}")
    
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
            # filter the PPI DataFrame to keep only rows where the SwissProt ID is in the valid UniProt ID set
            df = df[df['SwissProt ID'].isin(valid_uniprot_ids)]
            print(f"DeepGO data filtered to {len(df)} annotations after applying UniProt ID filtering.")
        else:
            print(f"No UniProt ID TSV files found in {upid_folder}. Skipping UniProt ID filtering.")

    # separate by function type (BP, MF, CC)
    function_types = df['Function Type'].unique()
    for function_type in function_types:
        subset = df[df['Function Type'] == function_type]

        # remove the function type column
        subset = subset.drop(columns=['Function Type'])

        function_shorthand_dict = {'Cellular Component': 'CC', 'Molecular Function': 'MF', 'Biological Process': 'BP'}
        shorthand = function_shorthand_dict.get(function_type, function_type)

        # change column names to UniProt_ID, DeepGO_Term, Score
        subset = subset.rename(columns={'SwissProt ID': 'UniProt_ID', 'GO Term': f"DeepGO_{shorthand}", 'Score': 'Score'}) 

        # save as separate tsvs with columns: ProteinID, GO_Term, Score
        output_name = output_file.replace(".tsv", f"_{shorthand}.tsv")
        subset.to_csv(output_name, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process raw DeepGO annotations into edge list format for graph building.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file containing raw DeepGO annotations.')
    parser.add_argument('output_file', type=str, help='Base path to save the processed edge list TSV files (without extension).')
    parser.add_argument('--upid_folder', type=str, default=None, help='Path to the folder containing UniProt IDs TSV files for ID filtering.')
    args = parser.parse_args()

    process_deepgo_annotations(args.input_file, args.output_file, args.upid_folder)
    print(f"Processed DeepGO annotations saved to {args.output_file.replace('.tsv', '')}_<FunctionType>.tsv for each function type.")