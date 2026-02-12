# script is for applying PCA to reduce dimensionality of embeddings before clustering

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer
import matplotlib.pyplot as plt
import argparse
import time


def apply_pca_reduction(embedding_dict, n_components = 100):
    '''Applies PCA to reduce the dimensionality of the embeddings in a dict.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))

    # normalize embeddings before PCA
    normalizer = Normalizer()
    normalized_embeddings = normalizer.fit_transform(embeddings)

    pca = PCA(n_components=n_components)  # use MLE to automatically determine the number of components to retain 90% variance
    reduced_embeddings = pca.fit_transform(normalized_embeddings)
    print(f"PCA reduced embeddings from {embeddings.shape[1]} to {reduced_embeddings.shape[1]} dimensions.")
    # change float precision to save space
    reduced_embeddings = reduced_embeddings.astype(np.float16)
    return dict(zip(protein_ids, reduced_embeddings))

def find_optimal_number_of_components(embedding_dict, variance_threshold = 0.90):
    '''Determines the optimal number of PCA components by generating a scree plot and finding the elbow point.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))
    pca = PCA().fit(embeddings)
    plt.figure()
    plt.plot(np.cumsum(pca.explained_variance_ratio_), marker='o')
    plt.title('Cumulative Explained Variance by PCA Components')
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.grid()
    plt.savefig('pca_scree_plot.png')
    optimal_components = np.argmax(np.cumsum(pca.explained_variance_ratio_) >= variance_threshold) + 1  # +1 because index starts at 0
    print(f"Optimal number of PCA components to retain {variance_threshold*100}% variance: {optimal_components}")
    return optimal_components

def extract_embeddings_from_npz(npz_filepath):
    '''Extracts embeddings from a .npz file and returns them as a dictionary.'''
    data = np.load(npz_filepath, allow_pickle=True)
    embedding_dict = {key: data[key] for key in data.files}
    # print(embedding_dict.keys())  # print keys to verify contents
    return embedding_dict

def save_embeddings_to_file(embedding_dict, output_filepath):
    '''Saves the embedding dictionary to a file in NumPy .npz format.'''
    np.savez_compressed(output_filepath, **embedding_dict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Apply PCA to reduce dimensionality of embeddings.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input .npz file containing embeddings.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save the PCA-reduced embeddings as a .npz file.')
    parser.add_argument('--variance_threshold', type=float, default=0.90, help='Variance threshold for determining optimal number of PCA components (default: 0.90).')
    args = parser.parse_args()

    print(f"Loading embeddings from {args.input_file}...")
    embedding_dict = extract_embeddings_from_npz(args.input_file)
    print(f"Loaded {len(embedding_dict)} embeddings at {len(embedding_dict[next(iter(embedding_dict))])} dimensions.")

    print("Applying PCA for dimensionality reduction...")
    optimal_components = find_optimal_number_of_components(embedding_dict, variance_threshold=args.variance_threshold)  # optional: find optimal components and show scree plot
    reduced_embedding_dict = apply_pca_reduction(embedding_dict, n_components=optimal_components)

    print(f"Saving PCA-reduced embeddings to {args.output_file}...")
    save_embeddings_to_file(reduced_embedding_dict, args.output_file)