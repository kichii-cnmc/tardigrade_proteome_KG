# script to access the PROSITE database and pull domain information for the proteins in our dataset, then save the results to a TSV file

import requests
import pandas as pd
import argparse
import time
import json

def parse_input_tsv(input_filepath):
    '''Parses the input TSV file and returns a list of Protein IDs and column name of the 2nd column.'''
    df = pd.read_csv(input_filepath, sep='\t')
    protein_ids = df.iloc[:, 0].tolist()  # get the first column as a list of protein IDs
    return protein_ids , df.columns[1]  # return the column name of the 2nd column for later use

def query_prosite_for_protein(protein_id):
    '''Queries the PROSITE database for domain information for a given protein ID.'''
    url = f"https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq={protein_id}&output=json&lowscore=1"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # raise an error for bad status codes
        return response.text
    except requests.RequestException as e:
        print(f"Error querying PROSITE for {protein_id}: {e}")
        return None
    
def process_prosite_response_json(protein_id, response_text):
    '''Processes the PROSITE response JSON and extracts relevant domain information.'''
    try:
        data = json.loads(response_text)
        domain_info_list = []
        for match in data.get('matchset', []):
            # collect score and level OR level_tag if score is not available
            domain_info = {
                'ProteinID': protein_id,
                'PROSITE_AC': match.get('signature_ac', ''),
                'Description': match.get('signature_id', ''),
                'Score': match.get('score', ''),
                'Level': match.get('level', ''),
                'LevelTag': match.get('level_tag', '')
            }
            domain_info_list.append(domain_info)
        return domain_info_list
    except ValueError as e:
        print(f"Error processing PROSITE response for {protein_id}: {e}")
        return []

def save_domain_info_to_tsv(domain_info_list, output_filepath):
    '''Saves the extracted domain information to a TSV file.'''
    if not domain_info_list:
        print("No domain information to save.")
        return
    
    df = pd.DataFrame(domain_info_list)
    df.to_csv(output_filepath, mode='a', sep='\t', index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Collect PROSITE domain information for a list of proteins.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file containing protein IDs.')
    parser.add_argument('output_file', type=str, help='Path to save the extracted domain information as a TSV file.')
    parser.add_argument('--test_mode', action='store_true', help='If set, processes only a small subset of data for testing.')

    args = parser.parse_args()

    print(f"\nParsing input file {args.input_file}...")
    protein_ids, second_column_name = parse_input_tsv(args.input_file)
    print(f"Found {len(protein_ids)} protein IDs. Second column name: {second_column_name}")

    if args.test_mode:
        protein_ids = protein_ids[:50]  # process only the first 50 protein IDs for testing
        print("Test mode enabled. Processing only the first 50 protein IDs.")

    # create output file and write header
    with open(args.output_file, 'w') as f:
        f.write('ProteinID\tPROSITE_AC\tDescription\tScore\tLevel\tLevelTag\n')

    start_time = time.time()
    for protein_id in protein_ids:
        print(f"Querying PROSITE for {protein_id}...")
        response_text = query_prosite_for_protein(protein_id)
        if response_text is None:
            print(f"Skipping {protein_id} due to query error.")
            continue
        # print(response_text)
        domain_info_list = process_prosite_response_json(protein_id, response_text)
        if domain_info_list:
            save_domain_info_to_tsv(domain_info_list, args.output_file)
        time.sleep(0.1)  # add a small delay to avoid overwhelming the server
    end_time = time.time()
    print(f"PROSITE domain information collection completed in {end_time - start_time:.2f} seconds for {len(protein_ids)} proteins. Results saved to {args.output_file}.")