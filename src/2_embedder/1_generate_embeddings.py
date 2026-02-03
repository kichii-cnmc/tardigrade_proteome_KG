# Script used to generate ESM-2 embeddings for protein sequences.

import torch
import esm
import time # measure time taken to genereate embeddings per protein
import os
import pandas as pd
import transformers
import numpy as np

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
    return model, batch_converter

def get_long_sequence_embedding(sequence, model, batch_converter):
    '''Generates an embedding for a long protein sequence using a sliding window approach.'''
    # Hard limit for ESM-2
    max_len = 1022 # 1024 minus <cls> and <eos>
    overlap = 256
    
    # if short process normally
    if len(sequence) <= max_len:
        data = [("temp", sequence)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33])
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
            results = model(batch_tokens, repr_layers=[33])
        
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
            results = model(batch_tokens, repr_layers=[33])
        token_representations = results["representations"][33]
        return token_representations[0, 1:-1, :].mean(dim=0).cpu().numpy()
    
    combined = torch.cat(residue_embeddings, dim=0)
    return combined.mean(dim=0).numpy()

def build_embedding_dict(df_list, model, batch_converter, sequence_column='Sequence', id_column='UniProt_ID'):
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
            embedding = get_long_sequence_embedding(sequence, model, batch_converter)
            embedding_dict[protein_id] = embedding
    return embedding_dict

def save_embeddings_to_file(embedding_dict, output_filepath):
    '''Saves the embedding dictionary to a file in NumPy .npz format.'''
    np.savez_compressed(output_filepath, **embedding_dict)

if __name__ == "__main__":
    # Load ESM-2 model
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    # Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
    data = [
        ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
        ("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE")
    ]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Measure time taken to generate embeddings
    start_time = time.time()

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    sequence_representations = []
    for i, tokens_len in enumerate(batch_lens):
        sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken to generate embeddings for {len(data)} proteins: {elapsed_time:.2f} seconds")

    print("Generated ESM-2 embeddings for the following proteins:")
    for i, label in enumerate(batch_labels):
        print(f"{label}: {sequence_representations[i].shape}")

    # test build_embedding_dict function to see if it gives the same reuslts as above
    df_test = pd.DataFrame(data, columns=['UniProt_ID', 'Sequence'])
    df_list = [df_test]
    embedding_dict = build_embedding_dict(df_list, model, batch_converter)
    for label in batch_labels:
        embedding = embedding_dict[label]
        print(f"{label}: {embedding.shape}")
        if label in embedding_dict:
            emb = embedding_dict[label]
            print(f"{label}: {emb.shape}")
        else:
            print(f"{label} not found in embedding dictionary.")
    
    # check if the embeddings match
    for i, label in enumerate(batch_labels):
        emb1 = sequence_representations[i].cpu().numpy()
        emb2 = embedding_dict[label]
        if np.allclose(emb1, emb2):
            print(f"Embeddings match for {label}.")
        else:
            print(f"Embeddings do NOT match for {label}.")