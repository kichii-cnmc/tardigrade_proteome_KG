# Script used to take the compressed ESM-2 embeddings and calculate clusters using PCA and K-means.

# PCA is applied, cluster by kmeans cosine similarity, visualize clusters with UMAP, save cluster assignments and centroids for later use in the KG construction step.
# Note: Embeddings were normalized before PCA.

import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn import metrics
from scipy.spatial.distance import cdist

def extract_embeddings_from_npz(npz_filepath):
    '''Extracts embeddings from a .npz file and returns them as a dictionary.'''
    data = np.load(npz_filepath, allow_pickle=True)
    embedding_dict = {key: data[key] for key in data.files}
    # print(embedding_dict.keys())  # print keys to verify contents
    return embedding_dict

def calculate_clusters(embedding_dict, n_clusters = 10):
    '''Calculates clusters from the PCA-reduced embeddings using K-means.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(embeddings)
    return dict(zip(protein_ids, cluster_labels)), kmeans.cluster_centers_

def save_cluster_assignments(cluster_labels, output_filepath):
    '''Saves the cluster assignments to a CSV file.'''
    cluster_df = pd.DataFrame(list(cluster_labels.items()), columns=['Protein_ID', 'Cluster_Label'])
    cluster_df.to_csv(output_filepath, index=False)

def calculate_silhouette_score(embedding_dict, cluster_labels):
    '''Calculates the silhouette score for the clustering.'''
    embeddings = np.array(list(embedding_dict.values()))
    labels = np.array([cluster_labels[pid] for pid in embedding_dict.keys()])
    score = silhouette_score(embeddings, labels)
    return score

def determine_optimal_clusters_ss(embedding_dict, min_clusters = 2,max_clusters=20, increment = 1):
    '''Determines the optimal number of clusters using silhouette scores.'''
    embeddings = np.array(list(embedding_dict.values()))
    silhouette_scores = []
    cluster_range = range(min_clusters, max_clusters + 1, increment)
    for n_clusters in cluster_range:
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        labels = kmeans.fit_predict(embeddings)
        score = silhouette_score(embeddings, labels)
        silhouette_scores.append(score)
    
    # Plotting the silhouette scores
    plt.figure()
    plt.plot(cluster_range, silhouette_scores, marker='o')
    plt.title('Silhouette Score vs Number of Clusters')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.grid()
    #plt.show()

    optimal_clusters = cluster_range[np.argmax(silhouette_scores)]

    # save the distribution of silhouette scores to a file
    silhouette_df = pd.DataFrame({
        'n_clusters': list(cluster_range),
        'silhouette_score': silhouette_scores
    })
    silhouette_df.to_csv('silhouette_scores.csv', index=False)

    return optimal_clusters

def determine_optimal_clusters_em(embedding_dict, min_clusters = 2,max_clusters=20, increment = 1):
    '''Determines the optimal number of clusters using the elbow method.'''
    embeddings = np.array(list(embedding_dict.values()))
    distortions = []
    inertia_scores = []
    cluster_range = range(min_clusters, max_clusters + 1, increment)
    for n_clusters in cluster_range:
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        kmeans.fit(embeddings)
        distortions.append(sum(np.min(cdist(embeddings, kmeans.cluster_centers_, 'euclidean'), axis=1)) / embeddings.shape[0])
        inertia_scores.append(kmeans.inertia_)
    # Plotting the elbow method results
    plt.figure()
    plt.plot(cluster_range, distortions, marker='o')
    plt.title('Elbow Method: Distortion vs Number of Clusters')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Distortion')
    plt.grid()
    plt.savefig('elbow_method_distortion.png')
    plt.figure()
    plt.plot(cluster_range, inertia_scores, marker='o')
    plt.title('Elbow Method: Inertia vs Number of Clusters')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Inertia')
    plt.grid()
    plt.savefig('elbow_method_inertia.png')
    # Determine optimal clusters using the elbow method (looking for the "elbow" point)
    optimal_clusters = cluster_range[np.argmin(np.diff(distortions))]
    return optimal_clusters

if __name__ == "__main__":
    # input file, cluster range (optional).
    parser = argparse.ArgumentParser(description='Calculate clusters from PCA-reduced embeddings.')
    parser.add_argument('input_file', type=str, help='Path to the input .npz file containing PCA-reduced embeddings.')
    parser.add_argument('--min_clusters', type=int, default=2, help='Minimum number of clusters to evaluate (default: 2).')
    parser.add_argument('--max_clusters', type=int, default=20, help='Maximum number of clusters to evaluate (default: 20).')
    parser.add_argument('--increment', type=int, default=1, help='Increment for the number of clusters to evaluate (default: 1).')
    args = parser.parse_args()

    print(f"Loading embeddings from {args.input_file}...")
    embedding_dict = extract_embeddings_from_npz(args.input_file)
    print(f"Loaded {len(embedding_dict)} embeddings at {len(embedding_dict[next(iter(embedding_dict))])} dimensions.")
    optimal_cluster_ct_em = determine_optimal_clusters_ss(embedding_dict, args.min_clusters, args.max_clusters, args.increment)
    print(f"Optimal number of clusters determined by elbow method: {optimal_cluster_ct_em}")

    # save_cluster_assignments(cluster_labels, 'cluster_assignments.csv')
    # print("Cluster assignments saved to cluster_assignments.csv")




