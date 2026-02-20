# script to directly compare embeddings to each other using cosine similarity and output the similarity matrix as a .csv file

from time import time
import numpy as np
import pandas as pd
import argparse
from sklearn.preprocessing import Normalizer
from sklearn.metrics.pairwise import cosine_similarity

def extract_embeddings_from_npz(npz_filepath):
    '''Extracts embeddings from a .npz file and returns them as a dictionary.'''
    data = np.load(npz_filepath, allow_pickle=True)
    embedding_dict = {key: data[key] for key in data.files}
    return embedding_dict

def calculate_cosine_similarity_matrix(embedding_dict, normalize=True):
    '''Calculates the cosine similarity matrix for the given embedding dictionary.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))
    
    # sklearn's cosine_similarity handles normalization internally
    similarity_matrix = cosine_similarity(embeddings)
    
    return protein_ids, similarity_matrix

def flatten_matrix(matrix):
    '''Flattens a 2D NxN matrix into an edge list format (protein1, protein2, similarity), excluding self-similarity.'''
    # Create indices for all pairs
    i_indices, j_indices = np.meshgrid(range(len(matrix)), range(len(matrix)), indexing='ij')
    
    # Create mask to exclude diagonal (self-similarity)
    mask = i_indices != j_indices
    
    # Extract indices and values using the mask
    i_flat = i_indices[mask]
    j_flat = j_indices[mask]
    values_flat = matrix[mask]
    
    # Stack into edge list format
    edge_list = np.column_stack((i_flat, j_flat, values_flat))
    return edge_list

def save_similarity_matrix_to_tsv(protein_ids, similarity_matrix, output_filepath):
    '''Saves the cosine similarity matrix to a TSV file with protein IDs as headers.'''
    df = pd.DataFrame(similarity_matrix, index=protein_ids, columns=protein_ids)
    df.to_csv(output_filepath, sep='\t', float_format='%.4f')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate cosine similarity matrix from embeddings.')
    parser.add_argument('input_file', type=str, help='Path to the input .npz file containing embeddings.')
    parser.add_argument('output_dir', type=str, help='Path to save the cosine similarity matrix as a .tsv file and edge list format.')
    parser.add_argument('--normalize', action='store_true', help='Whether to normalize embeddings before calculating similarity (default: False).')
    parser.add_argument('--test_mode', action='store_true', help='If set, processes only a small subset of data for testing.')

    args = parser.parse_args()
    output_csv_path = f"{args.output_dir}/embedding_similarity_matrix.tsv"
    output_list_path = f"{args.output_dir}/embedding_similarity_edge_list.tsv"

    print(f"\nLoading embeddings from {args.input_file}...")
    embedding_dict = extract_embeddings_from_npz(args.input_file)
    print(f"Loaded {len(embedding_dict)} embeddings at {len(embedding_dict[next(iter(embedding_dict))])} dimensions.")
    
    if args.test_mode:
        # Process only a small subset of data for testing
        test_size = min(100, len(embedding_dict))
        embedding_dict = {k: embedding_dict[k] for k in list(embedding_dict.keys())[:test_size]}
        print(f"Test mode enabled: processing only {test_size} embeddings.")

    start_time = time()
    print("Calculating cosine similarity matrix...")
    protein_ids, similarity_matrix = calculate_cosine_similarity_matrix(embedding_dict, normalize=args.normalize)
    end_time = time()

    print(f"Calculated cosine similarity matrix for {len(protein_ids)} proteins in {end_time - start_time:.2f} seconds.")

    print(f"Cosine similarity matrix shape: {similarity_matrix.shape}")

    # show quick histogram of similarity values to verify distribution
    import matplotlib.pyplot as plt
    plt.hist(similarity_matrix.flatten(), bins=50, range=(-1, 1))
    plt.title('Histogram of Cosine Similarity Values')
    plt.xlabel('Cosine Similarity')
    plt.ylabel('Frequency')
    # plt.show()

    print(f"Saving cosine similarity matrix to {output_csv_path}...")
    # Use numpy directly for faster I/O, especially for large matrices
    np.savetxt(output_csv_path, similarity_matrix.astype(np.float32), delimiter='\t', 
               header='\t'.join(protein_ids), comments='', fmt='%.4f')
    
    print("Cosine similarity calculation completed.")

    # save the matrix in edge list format as well for easier loading in graph databases
    time_start = time()
    print("Saving edge list format of similarity matrix...")
    
    # get entire matrix (bidirectional) excluding self-similarity & negative similarities
    indices = np.where((~np.eye(similarity_matrix.shape[0], dtype=bool)) & (similarity_matrix > 0))  # get indices of non-diagonal elements with positive similarity
    
    # Create edge list directly from indices and similarity values
    edge_list_df = pd.DataFrame({
        'protein1': [protein_ids[i] for i in indices[0]],
        'protein2': [protein_ids[j] for j in indices[1]], 
        'embedding_similarity': similarity_matrix[indices]
    })
    
    edge_list_df.to_csv(output_list_path, sep='\t', index=False, float_format='%.4f')
    time_end = time()
    print(f"Saved edge list format in {time_end - time_start:.2f} seconds.")