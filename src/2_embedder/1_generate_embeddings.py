# Script used to generate ESM-2 embeddings for protein sequences.

import torch
import esm
import time # measure time taken to genereate embeddings per protein
import os
import pandas as pd
import transformers
import numpy as np
import argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer
import matplotlib.pyplot as plt

def generate_df_list_of_sequences(folder_filepath, column_name = 'Sequence'):
    '''Generates a list of DataFrames containing protein sequences from TSV files in the specified folder.'''
    df_list = []
    for filename in os.listdir(folder_filepath):
        if filename.endswith('.tsv'):
            file_path = os.path.join(folder_filepath, filename)
            df = pd.read_csv(file_path, sep='\t')
            if len(df.columns) > 1 and column_name == df.columns[1]:  # ensure at least two columns before accessing the second
                df_list.append(df)
    return df_list

def initialize_esm_model(model_name='esm2_t33_650M_UR50D'):
    '''Initializes and returns the specified ESM model and its batch converter.'''
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    # move to GPU if available
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    print(f"Using device: {device}")

    return model, batch_converter, device

def get_long_sequence_embedding(sequence, model, batch_converter, device):
    '''Generates an embedding for a long protein sequence using a sliding window approach.'''
    # Hard limit for ESM-2
    max_len = 1022 # 1024 minus <cls> and <eos>
    overlap = 256
    
    # if short process normally
    if len(sequence) <= max_len:
        data = [("temp", sequence)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        
        with torch.no_grad():
            results = model(batch_tokens.to(device), repr_layers=[33])
        token_representations = results["representations"][33]
        # Remove special tokens and take mean
        return token_representations[0, 1:-1, :].mean(dim=0).cpu().numpy()

    # if long, use a sliding window
    residue_embeddings = []
    
    # step through the sequence with an overlap
    for i in range(0, len(sequence), max_len - overlap):
        chunk = sequence[i : i + max_len]
        if len(chunk) < 50: # skip tiny fragments at the very end
            continue
            
        data = [("temp", chunk)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        
        with torch.no_grad():
            results = model(batch_tokens.to(device), repr_layers=[33])
        
        token_representations = results["representations"][33]
        # remove special tokens and store residue embeddings
        chunk_emb = token_representations[0, 1:-1, :] # remove special tokens
        residue_embeddings.append(chunk_emb.cpu())

    # combine all chunks and take the global mean
    if not residue_embeddings:
        # Fall back to processing the full sequence
        data = [("temp", sequence)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        
        with torch.no_grad():
            results = model(batch_tokens.to(device), repr_layers=[33])
        token_representations = results["representations"][33]
        return token_representations[0, 1:-1, :].mean(dim=0).cpu().numpy()
    
    combined = torch.cat(residue_embeddings, dim=0)
    return combined.mean(dim=0).numpy()

def build_embedding_dict(df_list, model, batch_converter, device, sequence_column='Sequence', id_column='UniProt_ID'):
    '''Builds a dictionary of protein embeddings from a list of DataFrames.'''
    embedding_dict = {}
    for df in df_list:
        for index, row in df.iterrows():
            protein_id = row[id_column]
            sequence = row[sequence_column]
            if pd.isna(sequence) or not isinstance(sequence, str) or len(sequence) == 0:
                continue
            
            print(f"Processing {protein_id}...")
            # Generate embedding
            embedding = get_long_sequence_embedding(sequence, model, batch_converter, device)
            embedding_dict[protein_id] = embedding
    return embedding_dict

def save_embeddings_to_file(embedding_dict, output_filepath):
    '''Saves the embedding dictionary to a file in NumPy .npz format.'''
    np.savez_compressed(output_filepath, **embedding_dict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate ESM-2 embeddings for protein sequences.")
    parser.add_argument("input_folder", help="Path to the folder containing TSV files with protein sequences.")
    parser.add_argument("output_file", default = "embeddings.npz", help="Path to the output .npz file to save embeddings.")
    parser.add_argument("--test_mode", action='store_true', help="If set, processes only a small subset of data for testing.")
    args = parser.parse_args()

    print("\nLoading protein sequences from TSV files...")
    df_list = generate_df_list_of_sequences(args.input_folder)
    print(f"Found {len(df_list)} sequence TSV files.")

    print(f"Initializing ESM model...")
    model, batch_converter, device = initialize_esm_model()

    if args.test_mode:
        print("Test mode enabled: limiting to first 100 sequences per file. Total 100 sequences.")
        df_list = [df.head(100) for df in df_list]
    
    start_time = time.time()
    print("Generating embeddings for protein sequences...")
    embedding_dict = build_embedding_dict(df_list, model, batch_converter, device=device)
    end_time = time.time()
    print(f"Generated embeddings for {len(embedding_dict)} proteins in {end_time - start_time:.2f} seconds.")

    print(f"Saving embeddings to {args.output_file}...")
    # if no directory present, create it
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    save_embeddings_to_file(embedding_dict, args.output_file)

    print("Embedding generation completed.")