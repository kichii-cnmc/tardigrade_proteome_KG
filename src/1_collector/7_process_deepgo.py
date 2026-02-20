# script converts raw deepgo output into a edge list format for graph building, integrating scores.

import pandas as pd
import argparse

def process_deepgo_annotations(input_file, output_file):
    '''Processes raw DeepGO annotations into an edge list format for graph building.'''
    # add to df
    df = pd.read_csv(input_file, sep='\t')

    # confirm that columns are correct
    expected_columns = ['SwissProt ID', 'Function Type', 'GO Term', 'Score']
    if not all(col in df.columns for col in expected_columns):
        raise ValueError(f"Input file must contain the following columns: {expected_columns}")

    # separate by function type (BP, MF, CC)
    function_types = df['Function Type'].unique()
    for function_type in function_types:
        subset = df[df['Function Type'] == function_type]

        # remove the function type column
        subset = subset.drop(columns=['Function Type'])

        # change column names to UniProt_ID, DeepGO_Term, Score
        subset = subset.rename(columns={'SwissProt ID': 'UniProt_ID', 'GO Term': 'DeepGO_Term'}) 

    # save as separate tsvs with columns: ProteinID, GO_Term, Score
        function_shorthand = {'Cellular Component': 'CC', 'Molecular Function': 'MF', 'Biological Process': 'BP'}
        output_name = output_file.replace(".tsv", f"_{function_shorthand.get(function_type, function_type)}.tsv")
        subset.to_csv(output_name, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process raw DeepGO annotations into edge list format for graph building.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file containing raw DeepGO annotations.')
    parser.add_argument('output_file', type=str, help='Base path to save the processed edge list TSV files (without extension).')
    args = parser.parse_args()

    process_deepgo_annotations(args.input_file, args.output_file)
    print(f"Processed DeepGO annotations saved to {args.output_file.replace('.tsv', '')}_<FunctionType>.tsv for each function type.")