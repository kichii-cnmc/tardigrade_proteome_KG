# check sequence identity using CD-HIT, merging sequences above a given threshold 100%
# merges were checked with Clustal Omega MSA to confirm identical sequences
import os
import argparse
import subprocess
import pandas as pd

def combine_fastas(fasta_input_1, fasta_input_2, combined_fasta):
    """
    Combines two FASTA files into one.

    Parameters:
    fasta_input_1 (str): Path to the first input FASTA file.
    fasta_input_2 (str): Path to the second input FASTA file.
    combined_fasta (str): Path to the output combined FASTA file.
    """
    with open(combined_fasta, 'w') as outfile:
        for fname in [fasta_input_1, fasta_input_2]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def identify_accession_id(id_string):
    ''' Identifies whether the ID string is from, returns 0 for UniProt, 1 for others'''
    id_string = id_string.lstrip('>')
    if id_string.find('|') != -1:
        parts = id_string.split('|')
        if len(parts) >= 3:
            db = parts[0]
            accession = parts[1]
            if db in ['sp', 'tr']:
                return 0, accession
            else:
                return 1, accession
    else:
        return 1, id_string

def extract_accession_list(fasta_file):
    '''Extracts a list of accession IDs from a FASTA file.'''
    accession_list = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                id_string = line[1:].strip().split()[0]
                accession_list.append(id_string)
    return accession_list

def check_all_accessions_processed(fasta_input1, fasta_input2, accession_list):
    '''Checks that all accesssions across every step were processed.'''
    accession_set = set(accession_list)
    input_accessions_1 = set(extract_accession_list(fasta_input1))
    input_accessions_2 = set(extract_accession_list(fasta_input2))
    all_input_accessions = input_accessions_1.union(input_accessions_2)
    missing_accessions = all_input_accessions - accession_set
    if missing_accessions:
        print("Warning: The following accessions were not processed:")
        for acc in missing_accessions:
            print(acc)
    else:
        print("All accessions were processed successfully.")
        print(f"Input Accessions: {len(input_accessions_1)} and {len(input_accessions_2)}, Processed: {len(accession_set)}")

def parse_cdhit_cluster_file(cluster_file):
    """
    Parses CD-HIT cluster file (.clstr) and returns a DataFrame similar to pycdhit output.
    
    Parameters:
    cluster_file (str): Path to the .clstr file
    
    Returns:
    pd.DataFrame: DataFrame with columns ['cluster', 'identifier', 'size', 'identity']
    """
    clusters = []
    current_cluster = -1
    
    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[-1])
            elif line and not line.startswith('>'):
                # Parse sequence entry: 0	426aa, >WP_011923077.1... at 100.00%
                parts = line.split('\t')[1].split(',')
                length_part = parts[0].strip()
                size = int(length_part.replace('aa', '').replace('nt', ''))
                
                # Extract identifier
                id_part = parts[1].strip()
                if id_part.startswith('>'):
                    identifier = id_part[1:].split('...')[0]
                else:
                    identifier = id_part.split('...')[0]
                
                # Extract identity (if present)
                identity = 100.0  # Default for representative sequence
                if 'at' in line and '%' in line:
                    identity_str = line.split('at ')[-1].replace('%', '')
                    try:
                        identity = float(identity_str)
                    except ValueError:
                        identity = 100.0
                
                clusters.append({
                    'cluster': current_cluster,
                    'identifier': identifier,
                    'size': size,
                    'identity': identity
                })
    
    return pd.DataFrame(clusters)

def cluster_sequences(input_fasta, output_tsv, identity_threshold=1.0):
    """
    Clusters sequences in the input FASTA file using CD-HIT and writes the clustered sequences to the output FASTA file.

    Parameters:
    input_fasta (str): Path to the input FASTA file containing sequences to be clustered.
    output_tsv (str): Path to the output TSV file where cluster information will be saved.
    identity_threshold (float): Sequence identity threshold for clustering (default is 1.0).
    """
    # Run CD-HIT using subprocess
    temp_output = "temp.fasta"
    temp_cluster = temp_output + ".clstr"
    
    cmd = [
        "cd-hit",
        "-i", input_fasta,
        "-o", temp_output,
        "-c", str(identity_threshold),
        "-d", "0",  # sequence name length (0 = full)
        "-sc", "1"  # sort clusters by size
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("CD-HIT completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"CD-HIT failed with error: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise
    except FileNotFoundError:
        print("CD-HIT not found. Please install it using: micromamba install cd-hit")
        raise
    
    # Parse cluster file
    df_clstr = parse_cdhit_cluster_file(temp_cluster)
    
    # write tsv output, for each cluster value finds the longest sequence for each accession source
    accession_list = []
    with open(output_tsv, 'w') as out_f:
        out_f.write("UniProt_ID\tUniProt_Length\tNCBI_ID\tNCBI_Length\tOther_IDs\n")
        for cluster_id, group in df_clstr.groupby('cluster'):
            uniprot_id = ""
            uniprot_length = 0
            ncbi_id = ""
            ncbi_length = 0
            other_ids = []
            for idx, row in group.iterrows():
                id_type, accession = identify_accession_id(row['identifier'])
                if accession == "A0A1W0W9Q4":
                    print(id_type, row['identifier'], accession)
                accession_list.append(row['identifier'])
                if id_type == 0:  # UniProt
                    if row['size'] > uniprot_length:
                        if uniprot_id != "" and uniprot_id != accession:
                            other_ids.append(uniprot_id)
                        uniprot_id = accession
                        uniprot_length = row['size']
                    else:
                        other_ids.append(accession)
                elif id_type == 1:  # NCBI or others
                    if row['size'] > ncbi_length:
                        if ncbi_id != "" and ncbi_id != accession:
                            other_ids.append(ncbi_id)
                        ncbi_id = accession
                        ncbi_length = row['size']
                    else:
                        other_ids.append(accession)
                else:
                    other_ids.append(accession)
            if uniprot_id == "": # if no uniprot id, remove from list
                continue
            out_f.write(f"{uniprot_id}\t{uniprot_length}\t{ncbi_id}\t{ncbi_length}\t{';'.join(other_ids)}\n")
    
    return accession_list
        
if __name__ == "__main__":
    # use CLI arguments
    parser = argparse.ArgumentParser(description="Combine two FASTA files and cluster sequences using CD-HIT.")
    parser.add_argument("fasta1", help="Path to the first input FASTA file.")
    parser.add_argument("fasta2", help="Path to the second input FASTA file.")
    parser.add_argument("output", help="Path to the output clustered FASTA file.")
    parser.add_argument("--identity", type=float, default=1.0, help="Sequence identity threshold for clustering (default: 0.9).")
    args = parser.parse_args()  

    combined_fasta = "combined_temp.fasta"
    print("Combining FASTA files...")
    combine_fastas(args.fasta1, args.fasta2, combined_fasta)
    print("Clustering sequences with CD-HIT...")
    accession_list = cluster_sequences(combined_fasta, args.output, identity_threshold=args.identity)  
    print("Checking all accessions were processed...")
    check_all_accessions_processed(args.fasta1, args.fasta2, accession_list)
    print("Cleaning up temporary files...")
    if os.path.exists("temp.fasta"):
        os.remove("temp.fasta") 
    if os.path.exists("temp.fasta.clstr"):
        os.remove("temp.fasta.clstr")
    if os.path.exists(combined_fasta):
        os.remove(combined_fasta)
    print("Process completed.")