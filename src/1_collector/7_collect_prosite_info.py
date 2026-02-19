# script to access the PROSITE database and pull domain information for the proteins in our dataset, then save the results to a TSV file

import requests
import pandas as pd
import argparse
import time
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from dataclasses import dataclass
import asyncio
import aiohttp

# Global lock for thread-safe file operations
file_lock = Lock()

@dataclass
class PrositeConfig:
    concurrent_requests: int = 3  # Much lower for PROSITE
    rate_limit_delay: float = 1.0  # Longer delay between requests
    timeout: int = 30
    max_retries: int = 3
    retry_delay: float = 5.0
    batch_size: int = 50

def parse_input_tsv(input_filepath):
    '''Parses the input TSV file and returns a list of Protein IDs and column name of the 2nd column.'''
    df = pd.read_csv(input_filepath, sep='\t')
    protein_ids = df.iloc[:, 0].tolist()  # get the first column as a list of protein IDs
    return protein_ids , df.columns[1]  # return the column name of the 2nd column for later use

def query_prosite_for_protein(protein_id, session=None):
    '''Queries the PROSITE database for domain information for a given protein ID.'''
    if session is None:
        session = requests.Session()
    
    url = f"https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq={protein_id}&output=json&lowscore=1"
    try:
        response = session.get(url, timeout=15)
        response.raise_for_status()  # raise an error for bad status codes
        return response.text
    except requests.RequestException as e:
        print(f"Error querying PROSITE for {protein_id}: {e}")
        return None
    
def query_prosite_for_protein_with_retry(protein_id, session, config):
    '''Queries PROSITE with retry logic.'''
    for attempt in range(config.max_retries):
        try:
            response_text = query_prosite_for_protein(protein_id, session)
            if response_text:
                return response_text
        except Exception as e:
            if attempt < config.max_retries - 1:
                print(f"Retry {attempt + 1} for {protein_id} after error: {e}")
                time.sleep(config.retry_delay * (attempt + 1))
            else:
                print(f"Failed after {config.max_retries} attempts for {protein_id}")
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

def process_single_protein(protein_id, session):
    '''Process a single protein and return domain information.'''
    print(f"Querying PROSITE for {protein_id}...")
    response_text = query_prosite_for_protein_with_retry(protein_id, session, PrositeConfig())
    if response_text is None:
        print(f"Skipping {protein_id} due to query error.")
        return []
    
    domain_info_list = process_prosite_response_json(protein_id, response_text)
    return domain_info_list

def save_domain_info_to_tsv(domain_info_list, output_filepath):
    '''Saves the extracted domain information to a TSV file.'''
    if not domain_info_list:
        return
    
    with file_lock:
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
        protein_ids = protein_ids[:500]  # process only the first 50 protein IDs for testing
        print("Test mode enabled. Processing only the first 50 protein IDs.")

    # create output file and write header
    with open(args.output_file, 'w') as f:
        f.write('ProteinID\tPROSITE_AC\tDescription\tScore\tLevel\tLevelTag\n')

    start_time = time.time()
    
    # Use ThreadPoolExecutor for concurrent requests
    max_workers = 5  # Be respectful to the server
    batch_size = 50  # Process results in batches
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Create a session for each worker
        sessions = [requests.Session() for _ in range(max_workers)]
        
        # Submit all jobs
        future_to_protein = {}
        for i, protein_id in enumerate(protein_ids):
            session = sessions[i % max_workers]
            future = executor.submit(process_single_protein, protein_id, session)
            future_to_protein[future] = protein_id
        
        # Collect results in batches
        completed_count = 0
        all_domain_info = []
        
        for future in as_completed(future_to_protein):
            protein_id = future_to_protein[future]
            try:
                domain_info_list = future.result()
                all_domain_info.extend(domain_info_list)
                completed_count += 1
                
                # Save in batches to reduce I/O
                if len(all_domain_info) >= batch_size or completed_count == len(protein_ids):
                    if all_domain_info:
                        save_domain_info_to_tsv(all_domain_info, args.output_file)
                        all_domain_info = []
                
                if completed_count % 10 == 0:
                    print(f"Processed {completed_count}/{len(protein_ids)} proteins...")
                    
            except Exception as e:
                print(f"Error processing {protein_id}: {e}")
        
        # Close all sessions
        for session in sessions:
            session.close()
    
    end_time = time.time()
    print(f"PROSITE domain information collection completed in {end_time - start_time:.2f} seconds for {len(protein_ids)} proteins. Results saved to {args.output_file}.")