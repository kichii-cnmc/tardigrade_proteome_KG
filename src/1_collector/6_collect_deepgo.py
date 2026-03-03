# script for REST API calls to DeepGOPlus for scored GO term predictions

import argparse
import csv
import os
import requests
import time
import json
import pandas as pd

def identify_sequence_files(input_dir):
    '''Identifies sequence TSV files in the input directory'''
    # lists all sequence files in the input directory, checks their 2nd column title looking for "Sequence"
    sequence_files = []
    for filename in os.listdir(input_dir):
        if filename.endswith('.tsv'):
            file_path = os.path.join(input_dir, filename)
            try:
                df = pd.read_csv(file_path, sep='\t')
                if 'Sequence' in df.columns:
                    sequence_files.append(file_path)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
    return sequence_files

def pull_sequences(file_path):
    '''Pulls protein_ID and sequences from a TSV file'''
    df = pd.read_csv(file_path, sep='\t')
    protein_ids = df['UniProt_ID'].tolist()
    sequences = df['Sequence'].tolist()
    id_sequence_set = set(zip(protein_ids, sequences))
    return id_sequence_set

def form_fasta_string(id_sequence_set):
    '''Forms a FASTA formatted string from the set of (protein_ID, sequence) tuples'''
    fasta_string = ""
    for protein_id, sequence in id_sequence_set:
        if protein_id == "0" or pd.isna(protein_id) or pd.isna(sequence):
            continue  # Skip invalid entries
        fasta_string += f">{protein_id}\n{sequence}\n"
    return fasta_string

def process_deepgo_to_tsv(json_response, output_file="deepgo_results.tsv"):
    """
    Parses DeepGO JSON response and writes a TSV file with columns:
    SwissProt ID, Function Type, GO ID, GO Name, Score.
    """
    
    # Define the headers for our TSV file
    headers = ["SwissProt ID", "Function Type", "GO Term", "Score"]
    
    file_exists = os.path.isfile(output_file)
    
    try:
        with open(output_file, mode='a', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, delimiter='\t')
            
            # Write the header row
            if not file_exists:
                writer.writerow(headers)
            
            # Loop 1: Iterate through each protein prediction
            for prediction in json_response.get("predictions", []):
                # 'protein_info' usually contains the ID/Accession
                swissprot_id = prediction.get("protein_info", "Unknown")
                
                # Loop 2: Iterate through the 3 categories (Cellular, Molecular, Biological)
                for func_category in prediction.get("functions", []):
                    function_type = func_category.get("name")
                    
                    # Loop 3: Iterate through specific GO terms in that category
                    # The API returns a list like: ["GO:0110165", "cellular anatomical structure", 0.35]
                    for go_item in func_category.get("functions", []):
                        go_term = go_item[0] + " - " + go_item[1]  # GO ID and Name
                        score = go_item[2]
                        
                        # Write the flattened row
                        writer.writerow([swissprot_id, function_type, go_term, score])
                        
        print(f"Successfully processed {len(json_response['predictions'])} proteins into '{output_file}'")

    except IOError as e:
        print(f"File error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

def query_deepgoplus(id_sequence_set, output_file="deepgoplus_predictions.tsv", version_number="1.0.27"):
    '''Queries DeepGOPlus API for GO term predictions in batches, continuously appends to TSV file as results are received.'''
    
    url_create = "https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create"
    batch_size = 100
    id_sequence_list = list(id_sequence_set)

    # failure tracking
    consective_failures = 0
    max_consecutive_failures = 5
    
    for i in range(0, len(id_sequence_list), batch_size):
        sequences = id_sequence_list[i:i+batch_size]
        
        # Use the SAME payload structure as your working test function
        payload = {
            'version': version_number,
            'data_format': 'fasta',  # Make sure this matches test
            'data': form_fasta_string(sequences),
            'threshold': '0.1'  # Use string like in test function
        }
        
        try:
            print(f"Submitting batch {i//batch_size + 1}...")
            # print(f"Payload: {payload}")  # Debug print
            
            response = requests.post(url_create, json=payload, timeout=60)
            
            #print(f"Response status: {response.status_code}")
            #print(f"Response text: {response.text[:200]}...")  # Debug print
            
            response.raise_for_status()
            
            job_data = response.json()
            process_deepgo_to_tsv(job_data, output_file=output_file)  # Process initial response for any immediate results
            consective_failures = 0  # Reset on success
                
        except requests.exceptions.RequestException as e:
            print(f"Request failed for batch {i//batch_size + 1}: {e}")
            consective_failures += 1
            
            # Retry with exponential backoff
            max_retries = 3
            for retry_count in range(1, max_retries + 1):
                wait_time = 10 * (2 ** (retry_count - 1))  # 10, 20, 40 seconds
                print(f"Retrying batch {i//batch_size + 1} in {wait_time} seconds (attempt {retry_count}/{max_retries})...")
                time.sleep(wait_time)
                
                try:
                    response = requests.post(url_create, json=payload, timeout=60)
                    response.raise_for_status()
                    job_data = response.json()
                    process_deepgo_to_tsv(job_data, output_file=output_file)
                    print(f"Batch {i//batch_size + 1} succeeded on retry attempt {retry_count}")
                    consective_failures = 0  # Reset on success
                    break  # Success, exit retry loop
                except requests.exceptions.RequestException as retry_e:
                    print(f"Retry attempt {retry_count} failed: {retry_e}")
                    if retry_count == max_retries:
                        # All retries exhausted, log and continue
                        failed_proteins = [protein_id for protein_id, _ in sequences]
                        os.makedirs("logs", exist_ok=True)
                        with open("logs/failed_deepgo_batches.log", "a") as log_file:
                            log_file.write(f"Batch {i//batch_size + 1} failed after {max_retries} retries. "
                                f"Proteins: {', '.join(failed_proteins)}. Error: {e}\n")
                        print(f"Batch {i//batch_size + 1} permanently failed after {max_retries} retries")
                except json.JSONDecodeError as retry_e:
                    print(f"Invalid JSON response on retry {retry_count}: {retry_e}")
                    if retry_count == max_retries:
                        failed_proteins = [protein_id for protein_id, _ in sequences]
                        os.makedirs("logs", exist_ok=True)
                        with open("logs/failed_deepgo_batches.log", "a") as log_file:
                            log_file.write(f"Batch {i//batch_size + 1} failed with JSON decode error after {max_retries} retries. "
                                f"Proteins: {', '.join(failed_proteins)}\n")
                except Exception as retry_e:
                    print(f"Unexpected error on retry {retry_count}: {retry_e}")
                    if retry_count == max_retries:
                        failed_proteins = [protein_id for protein_id, _ in sequences]
                        os.makedirs("logs", exist_ok=True)
                        with open("logs/failed_deepgo_batches.log", "a") as log_file:
                            log_file.write(f"Batch {i//batch_size + 1} failed with unexpected error after {max_retries} retries. "
                                f"Proteins: {', '.join(failed_proteins)}. Error: {retry_e}\n")
        except json.JSONDecodeError as e:
            print(f"Invalid JSON response for batch {i//batch_size + 1}: {e}")
            continue
        
        # Small delay between batches to be nice to the API
        time.sleep(2)

def normalize_scores(file_path):
    '''Normalizes the scores in a TSV file to be between 0 and 1 by dividing by the maximum score.'''
    df = pd.read_csv(file_path, sep='\t')
    if len(df.columns) > 3 and df.columns[3] == 'Score':
        max_score = df['Score'].max()
        if max_score > 0:
            df['Score'] = df['Score'] / max_score
            df.to_csv(file_path, sep='\t', index=False)
            print(f"Scores normalized in {file_path}")
        else:
            print(f"No valid scores to normalize in {file_path}")
    else:
        print(f"'Score' column not found in {file_path}")

def test_api_endpoint():
    '''Test the API endpoint to see what it expects'''
    try:
        # Test with minimal data
        test_payload = {
            'version': '1.0.27',
            'data_format': 'fasta',
            'data': '>test\nMKVLWALLLVWLLLLVWGTVAKTTIGA\n>test2\nMKVLLLWALLLVWLLLLVWGTVAKTTIGA\n',
            'threshold': '0.1'
        }
        
        response = requests.post(
            "https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create",
            json=test_payload,
            timeout=30
        )
        
        print(f"Test response status: {response.status_code}")
        print(f"Test response: {response.text}")
        
    except Exception as e:
        print(f"Test failed: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect GO term predictions from DeepGOPlus for protein sequences in TSV files.")
    parser.add_argument("input_dir", type=str, help="Directory containing TSV files with protein sequences.")
    parser.add_argument("--output_file", type=str, default="deepgoplus_predictions.tsv", help="Output TSV file for GO term predictions.")
    parser.add_argument("--version_number", type=str, default="1.0.27", help="DeepGOPlus version number to use.")
    parser.add_argument("--test_mode", action="store_true", help="If set, runs in test mode with limited data.")
    args = parser.parse_args()

    # print("Testing API endpoint...")
    # test_api_endpoint()
    # exit(0)

    sequence_files = identify_sequence_files(args.input_dir)
    all_id_sequence_set = set()
    for file_path in sequence_files:
        id_sequence_set = pull_sequences(file_path)
        all_id_sequence_set.update(id_sequence_set)

    if args.test_mode:
        # Sort by protein ID for consistent, deterministic selection
        all_id_sequence_list = sorted(list(all_id_sequence_set), key=lambda x: x[0])  # Sort by protein_id (first element)
        test_sequences = all_id_sequence_list[:100] or []  # Take first 100 or all if less than 100
        all_id_sequence_set = set(test_sequences)
        print(f"Test mode: Using first 100 sequences (sorted by protein ID)")

    start_time = time.time()
    query_deepgoplus(all_id_sequence_set, output_file=args.output_file, version_number=args.version_number)
    end_time = time.time()
    normalize_scores(args.output_file)
    print(f"GO term prediction of {len(all_id_sequence_set)} sequences completed in {end_time - start_time:.2f} seconds. Results saved to {args.output_file}.")

