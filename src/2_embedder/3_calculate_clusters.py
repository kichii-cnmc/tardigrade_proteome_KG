# Script used to take the compressed ESM-2 embeddings and calculate clusters using PCA and K-means.

# PCA is applied, cluster by kmeans cosine similarity, visualize clusters with UMAP, save cluster assignments and centroids for later use in the KG construction step.
# Note: Embeddings were normalized before PCA.

import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score
from sklearn import metrics
from scipy.spatial.distance import cdist
from sklearn.preprocessing import Normalizer
import hdbscan

def extract_embeddings_from_npz(npz_filepath):
    '''Extracts embeddings from a .npz file and returns them as a dictionary.'''
    data = np.load(npz_filepath, allow_pickle=True)
    embedding_dict = {key: data[key] for key in data.files}
    # print(embedding_dict.keys())  # print keys to verify contents
    return embedding_dict

def normalize_embeddings(embedding_dict):
    '''Normalizes the embeddings in the dictionary using L2 normalization.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))
    normalizer = Normalizer()
    normalized_embeddings = normalizer.fit_transform(embeddings)
    return dict(zip(protein_ids, normalized_embeddings))

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

def x_means_clustering(embedding_dict, initial_clusters = 2, max_clusters = 20):
    '''Performs X-means clustering to automatically determine the optimal number of clusters.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))
    # Start with initial K-means clustering
    kmeans = KMeans(n_clusters=initial_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(embeddings)
    current_clusters = initial_clusters
    print(f"Starting X-means clustering with min {initial_clusters} clusters and max {max_clusters} clusters.")
    while current_clusters < max_clusters:
        # For each cluster, try splitting it and evaluate with AIC/BIC
        new_labels = cluster_labels.copy()
        for cluster in range(current_clusters):
            cluster_indices = np.where(cluster_labels == cluster)[0]
            if len(cluster_indices) <= 1:
                continue  # skip small clusters
            kmeans_split = KMeans(n_clusters=2, random_state=42)
            split_labels = kmeans_split.fit_predict(embeddings[cluster_indices])
            # Evaluate split with AIC/BIC)
            original_inertia = kmeans.inertia_
            split_inertia = kmeans_split.inertia_
            if split_inertia < original_inertia:  # simple criterion for accepting split
                new_labels[cluster_indices[split_labels == 0]] = cluster
                new_labels[cluster_indices[split_labels == 1]] = current_clusters
                current_clusters += 1
                print(f"Cluster {cluster} split into 2 clusters. Total clusters: {current_clusters}")
                if current_clusters >= max_clusters:
                    break
        # If no splits were accepted, break the loop
        if np.array_equal(cluster_labels, new_labels):
            break
        cluster_labels = new_labels
        current_clusters += 1
    print(f"X-means determined optimal number of clusters: {current_clusters}")
    return dict(zip(protein_ids, cluster_labels)), kmeans.cluster_centers_

def dbscan_clustering(embedding_dict, eps=0.25, min_samples=5):
    '''Performs DBSCAN clustering to automatically determine clusters based on density.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    cluster_labels = dbscan.fit_predict(embeddings)
    print(max(cluster_labels) + 1, "clusters found (including noise if present).")
    return dict(zip(protein_ids, cluster_labels)), None  # DBSCAN does not have cluster centers

def hdbscan_clustering(embedding_dict, min_cluster_size=5):
    '''Performs HDBSCAN clustering to automatically determine clusters based on density.'''
    protein_ids = list(embedding_dict.keys())
    embeddings = np.array(list(embedding_dict.values()))
    hdbscan_clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
    cluster_labels = hdbscan_clusterer.fit_predict(embeddings)
    print(max(cluster_labels) + 1, "clusters found (including noise if present).")
    return dict(zip(protein_ids, cluster_labels)), None  # HDBSCAN does not have cluster centers

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
    print("Normalizing embeddings...")
    embedding_dict = normalize_embeddings(embedding_dict)
    print("Calculating clusters using DBSCAN...")
    # for i in range(1, 20):
    #      print(f"DBSCAN with eps={0.1 * i}...")
    #      cluster_labels, _ = dbscan_clustering(embedding_dict, eps=0.1 * i)
    #      cluster_count = max(cluster_labels.values()) + 1  # +1 because cluster labels start at 0, and -1 is noise
    #      print(f"Calculated clusters for {cluster_count} clusters (including noise if present).")
    #      if cluster_count > 1:  # silhouette score is only valid if there are at least 2 clusters 
    #         silhouette_score_value = calculate_silhouette_score(embedding_dict, cluster_labels)
    #         print(f"Silhouette Score for DBSCAN with eps={0.1 * i}: {silhouette_score_value:.4f}")

    # try HDBSCAN 
    print("Calculating clusters using HDBSCAN...")
    cluster_labels, _ = hdbscan_clustering(embedding_dict, min_cluster_size=5)
    silhouette_score_value = calculate_silhouette_score(embedding_dict, cluster_labels)
    print(f"Silhouette Score for HDBSCAN with min_cluster_size=5: {silhouette_score_value:.4f}")

    # optimal_cluster_ct_em = determine_optimal_clusters_ss(embedding_dict, args.min_clusters, args.max_clusters, args.increment)
    # print(f"Optimal number of clusters determined by elbow method: {optimal_cluster_ct_em}")

    # save_cluster_assignments(cluster_labels, 'cluster_assignments.csv')
    # print("Cluster assignments saved to cluster_assignments.csv")




