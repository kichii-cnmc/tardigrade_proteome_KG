# Used to build knowledge graph from protein info TSV files or folder.

import pandas as pd
import networkx as nx
import argparse
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from collections import defaultdict, Counter
import seaborn as sns

def load_tsv_into_df(file_path, selection_percentage=100):
    '''Load protein info from TSV file into a pandas DataFrame.'''
    df = pd.read_csv(file_path, sep="\t")
    # keep only the top X% of rows if selection_percentage < 100
    if selection_percentage < 100:
        num_rows = int(len(df) * (selection_percentage / 100))
        df = df.iloc[:num_rows].reset_index(drop=True)
    return df

def load_tsvs_from_folder(folder_path, selection_percentage=100):
    '''Load all TSV files from a folder into a list of pandas DataFrames.'''
    all_files = glob.glob(os.path.join(folder_path, "*.tsv"))
    df_list = []
    for file in all_files:
        df = pd.read_csv(file, sep="\t")
        # keep only the top X% of rows if selection_percentage < 100
        if selection_percentage < 100:
            num_rows = int(len(df) * (selection_percentage / 100))
            df = df.iloc[:num_rows].reset_index(drop=True)
        df_list.append(df)
    return df_list

def collect_header_names(df_list, n = 0):
    '''Collect a list of the nth header names from a set of DataFrames.'''
    header_names = []
    for df in df_list:
        if n < len(df.columns):
            header_names.append(df.columns[n])
    return header_names

def build_kg_target_triples_df_list(df_list, list_of_targets):
    '''List version of build_kg_target_triples_df, used to process multiple targets.'''
    all_triples_df = []
    failed_targets_ct = len(list_of_targets)
    target_names = [target[0] for target in list_of_targets]
    df_names = collect_header_names(df_list, 1)
    for i, df in enumerate(df_list):
        if df_names[i] in target_names:
            target_index = target_names.index(df_names[i])
            target, edge_type, weight = list_of_targets[target_index]
            triples_df = build_kg_target_triples_df(df, target, edge_type, weight)
            all_triples_df.append(triples_df)
            failed_targets_ct -= 1
    if failed_targets_ct > 0:
        print(f"Warning: {failed_targets_ct} items in list_of_targets were not found in the df_list")
    return all_triples_df

def build_kg_node_attr_df_list(df_list, list_of_node_attr):
    '''List version of build_kg_node_attr_df, used to process multiple node attributes.'''
    all_nodes_df = []
    failed_targets_ct = len(list_of_node_attr)
    node_attr_names = [node_attr[0] for node_attr in list_of_node_attr]
    df_names = collect_header_names(df_list, 1)
    for i, df in enumerate(df_list):
        if df_names[i] in node_attr_names:
            node_attr_index = node_attr_names.index(df_names[i])
            id_column, node_attr_value = list_of_node_attr[node_attr_index]
            nodes_df = build_kg_node_attr_df(df, id_column, node_attr_value)
            all_nodes_df.append(nodes_df)
            failed_targets_ct -= 1 
    if failed_targets_ct > 0:
        print(f"Warning: {failed_targets_ct} items in list_of_node_attr were not found in the df_list")
    return all_nodes_df

def build_kg_node_alias_attr_df_list(df_list, list_of_node_alias_attr):
    '''List version of build_kg_node_alias_attr_df, used to process multiple alias attributes.'''
    all_alias_df = []
    failed_targets_ct = len(list_of_node_alias_attr)
    alias_attr_names = [alias_attr_item[1] for alias_attr_item in list_of_node_alias_attr]
    df_names = collect_header_names(df_list, 1)
    for i, df in enumerate(df_list):
        if df_names[i] in alias_attr_names:
            alias_attr_index = alias_attr_names.index(df_names[i])
            id_column, alias_attr_column = list_of_node_alias_attr[alias_attr_index]
            alias_df = build_kg_node_alias_attr_df(df, id_column, alias_attr_column)
            all_alias_df.append(alias_df)
            failed_targets_ct -= 1
            print(f"Built alias attribute DataFrame for {df_names[i]}")
    if failed_targets_ct > 0:
        print(f"Warning: {failed_targets_ct} items in list_of_node_alias_attr were not found in the df_list")
    return all_alias_df

def build_kg_target_triples_df(df, target, edge_type = None, weight = None):
    '''Builds triples from the DataFrame based on a target column, edge type, and weight.'''
    triples_list = []
    for index, row in df.iterrows():
        source_id = row.iloc[0]  # assuming the first column is the source ID
        target_values = str(row[target]).split(';') if pd.notna(row[target]) else []
        for target_value in target_values:
            target_value = target_value.strip()
            if weight is None:
                weight = row.iloc[2] if len(row) > 2 else 1  # default weight from third column or 1
            if target_value:
                triples_list.append((source_id, target_value, edge_type, weight))
    triples_df = pd.DataFrame(triples_list, columns=['source', 'target', 'edge_type', 'weight'])
    return triples_df

def build_kg_node_attr_df(df, id_column, node_attr_value):
    '''Builds a DataFrame for assigning the same node_attr value to all ids under the id_column.'''
    nodes_set = set()
    for index, row in df.iterrows():
        id_values = str(row[id_column]).split(';') if pd.notna(row[id_column]) else []
        for id_value in id_values:
            id_value = id_value.strip()
            if id_value:
                nodes_set.add((id_value, node_attr_value))
    nodes_df = pd.DataFrame(list(nodes_set), columns=['id', 'node_attr'])
    return nodes_df

def build_kg_node_alias_attr_df(df, id_column, alias_attr_column):
    '''Builds a DataFrame for assigning an alias attribute to nodes of a specific type.'''
    alias_set = set()
    for index, row in df.iterrows():
        id_values = str(row[id_column]).split(';') if pd.notna(row[id_column]) else []
        for id_value in id_values:
            id_value = id_value.strip()
            if id_value:
                alias_set.add((id_value, row[alias_attr_column]))
    alias_df = pd.DataFrame(list(alias_set), columns=['id', 'alias_attr'])
    return alias_df

def build_knowledge_graph(triples_df_list, node_attr_df_list, node_alias_attr_df_list):
    '''Builds a NetworkX knowledge graph from triples and node attributes.'''
    G = nx.MultiDiGraph()
    
    # First, add all nodes from node attributes to ensure they exist
    for node_attr_df in node_attr_df_list:
        for index, row in node_attr_df.iterrows():
            if row['id'] not in G.nodes:
                G.add_node(row['id'])
            if 'node_attr' not in G.nodes[row['id']]:
                G.nodes[row['id']]['node_attr'] = []
            G.nodes[row['id']]['node_attr'].append(row['node_attr'])
    
    # Add alias attributes (also ensures nodes exist)
    for alias_attr_df in node_alias_attr_df_list:
        for index, row in alias_attr_df.iterrows():
            if row['id'] not in G.nodes:
                G.add_node(row['id'])
            if 'alias_attr' not in G.nodes[row['id']]:
                G.nodes[row['id']]['alias_attr'] = []
            G.nodes[row['id']]['alias_attr'].append(row['alias_attr'])
    
    # Then add triples as edges
    for triples_df in triples_df_list:
        for index, row in triples_df.iterrows():
            G.add_edge(row['source'], row['target'], edge_type=row['edge_type'], weight=row['weight'])
    
    return G

def get_node_colors_and_sizes(G, subgraph_nodes, center_node=None):
    """Generate colors and sizes for nodes based on their attributes and importance."""
    node_colors = []
    node_sizes = []
    
    # Color mapping for different node types
    color_map = {
        'Protein': '#FF6B6B',  # Red
        'GO': '#4ECDC4',       # Teal
        'KEGG': '#45B7D1',     # Blue  
        'Pfam': '#96CEB4',     # Green
        'PDB': '#FECA57',      # Yellow
        'STRING': '#FF9FF3',   # Pink
        'UniProt': '#54A0FF',  # Light blue
        'NCBI': '#5F27CD',     # Purple
        'DeepGO': '#00D2D3',   # Cyan
        'Embedding': '#FF5722', # Orange
        'AlphaFold': '#795548'  # Brown
    }
    
    for node in subgraph_nodes:
        # Determine node type and color
        node_attr = G.nodes[node].get('node_attr', 'Unknown')
        color = '#95A5A6'  # Default gray
        
        for node_type, type_color in color_map.items():
            if node_type in node_attr:
                color = type_color
                break
                
        node_colors.append(color)
        
        # Determine node size based on degree and centrality
        degree = G.degree(node)
        base_size = 300
        size = base_size + (degree * 50)  # Scale with degree
        
        # Highlight center node
        if center_node and node == center_node:
            size *= 2
            
        node_sizes.append(min(size, 2000))  # Cap maximum size
        
    return node_colors, node_sizes

def get_edge_colors_and_widths(G, subgraph):
    """Generate colors and widths for edges based on their types and weights."""
    edge_colors = []
    edge_widths = []
    
    edge_color_map = {
        'GO_annotation': '#3498DB',
        'KEGG_pathway': '#E74C3C', 
        'Pfam_domain': '#2ECC71',
        'PDB_structure': '#F39C12',
        'STRING_interaction': '#9B59B6',
        'sequence_similarity': '#1ABC9C',
        'embedding_similarity': '#E67E22',
        'DeepGO_prediction': '#34495E'
    }
    
    for edge in subgraph.edges(data=True):
        edge_type = edge[2].get('edge_type', 'unknown')
        weight = edge[2].get('weight', 1.0)
        
        # Color based on edge type
        color = edge_color_map.get(edge_type, '#7F8C8D')
        edge_colors.append(color)
        
        # Width based on weight
        width = max(0.5, min(5.0, float(weight) * 2))  # Scale and cap width
        edge_widths.append(width)
        
    return edge_colors, edge_widths

def create_legend(G, subgraph):
    """Create legend for node types and edge types."""
    # Node type legend
    node_types = set()
    for node in subgraph.nodes():
        node_attr = G.nodes[node].get('node_attr', 'Unknown')
        node_types.add(node_attr)
    
    # Edge type legend  
    edge_types = set()
    for edge in subgraph.edges(data=True):
        edge_type = edge[2].get('edge_type', 'unknown')
        edge_types.add(edge_type)
    
    return list(node_types), list(edge_types)

def visualize_subgraph_advanced(G, center_node=None, method='degree', radius=2, 
                               max_nodes=50, layout='spring', figsize=(16, 12),
                               show_labels=True, show_edge_labels=True, 
                               save_path=None, node_filter=None, edge_filter=None):
    """
    Advanced subgraph visualization with multiple customization options.
    
    Parameters:
    -----------
    G : networkx.Graph
        The full knowledge graph
    center_node : str, optional
        Central node to build subgraph around. If None, picks random node.
    method : str, default='degree'
        Method for subgraph selection: 'degree', 'proximity', 'community', 'centrality'
    radius : int, default=2
        For degree method: number of hops from center node
    max_nodes : int, default=50
        Maximum number of nodes to include in subgraph
    layout : str, default='spring'
        Layout algorithm: 'spring', 'circular', 'kamada_kawai', 'planar', 'spectral'
    figsize : tuple, default=(16, 12)
        Figure size
    show_labels : bool, default=True
        Whether to show node labels
    show_edge_labels : bool, default=True
        Whether to show edge labels
    save_path : str, optional
        Path to save the visualization
    node_filter : callable, optional
        Function to filter nodes (takes node, attributes as input)
    edge_filter : callable, optional
        Function to filter edges (takes edge data as input)
    """
    
    if G.number_of_nodes() == 0:
        print("The graph is empty. No nodes to visualize.")
        return None
    
    # Select center node if not provided
    if center_node is None:
        center_node = max(G.nodes(), key=lambda x: G.degree(x))  # Highest degree node
    elif center_node not in G.nodes():
        print(f"Node '{center_node}' not found in graph. Using random node.")
        center_node = max(G.nodes(), key=lambda x: G.degree(x))
    
    print(f"Building subgraph around node: {center_node}")
    
    # Build subgraph based on method
    if method == 'degree':
        # Include nodes within radius hops
        nodes_to_include = set([center_node])
        current_layer = {center_node}
        
        for _ in range(radius):
            next_layer = set()
            for node in current_layer:
                neighbors = set(G.neighbors(node))
                if G.is_directed():
                    neighbors.update(G.predecessors(node))
                next_layer.update(neighbors)
            nodes_to_include.update(next_layer)
            current_layer = next_layer
            
            if len(nodes_to_include) > max_nodes:
                break
                
    elif method == 'proximity':
        # Include closest nodes by shortest path
        distances = nx.single_source_shortest_path_length(G, center_node, cutoff=radius)
        nodes_to_include = set(list(distances.keys())[:max_nodes])
        
    elif method == 'centrality':
        # Include nodes with highest centrality scores
        centrality = nx.degree_centrality(G)
        sorted_nodes = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
        nodes_to_include = set([center_node] + [node for node, _ in sorted_nodes[:max_nodes-1]])
        
    elif method == 'community':
        # Include nodes from same community
        try:
            communities = nx.community.greedy_modularity_communities(G)
            center_community = None
            for community in communities:
                if center_node in community:
                    center_community = community
                    break
            nodes_to_include = set(list(center_community)[:max_nodes]) if center_community else {center_node}
        except:
            # Fallback to degree method
            print("Community detection failed. Using degree method.")
            method = 'degree'
            return visualize_subgraph_advanced(G, center_node, 'degree', radius, max_nodes, layout, figsize, show_labels, show_edge_labels, save_path, node_filter, edge_filter)
    
    # Apply node filter if provided
    if node_filter:
        nodes_to_include = {node for node in nodes_to_include 
                          if node_filter(node, G.nodes[node])}
    
    # Limit to max_nodes
    if len(nodes_to_include) > max_nodes:
        # Prioritize nodes closer to center
        distances = nx.single_source_shortest_path_length(G, center_node)
        nodes_to_include = set(sorted(nodes_to_include, 
                                    key=lambda x: distances.get(x, float('inf')))[:max_nodes])
    
    # Create subgraph
    subgraph = G.subgraph(nodes_to_include)
    
    # Apply edge filter if provided
    if edge_filter:
        edges_to_remove = []
        for edge in subgraph.edges(data=True):
            if not edge_filter(edge[2]):
                edges_to_remove.append((edge[0], edge[1]))
        subgraph = subgraph.copy()
        subgraph.remove_edges_from(edges_to_remove)
    
    print(f"Subgraph created with {subgraph.number_of_nodes()} nodes and {subgraph.number_of_edges()} edges")
    
    if subgraph.number_of_nodes() == 0:
        print("No nodes in subgraph after filtering.")
        return None
    
    # Create visualization
    fig, ax = plt.subplots(figsize=figsize)
    
    # Choose layout
    if layout == 'spring':
        pos = nx.spring_layout(subgraph, k=3, iterations=50, seed=42)
    elif layout == 'circular':
        pos = nx.circular_layout(subgraph)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(subgraph)
    elif layout == 'spectral':
        pos = nx.spectral_layout(subgraph)
    else:
        pos = nx.spring_layout(subgraph, seed=42)
    
    # Get node colors and sizes
    node_colors, node_sizes = get_node_colors_and_sizes(G, subgraph.nodes(), center_node)
    
    # Get edge colors and widths
    edge_colors, edge_widths = get_edge_colors_and_widths(G, subgraph)
    
    # Draw nodes
    nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors, 
                          node_size=node_sizes, alpha=0.8, ax=ax)
    
    # Draw edges
    nx.draw_networkx_edges(subgraph, pos, edge_color=edge_colors, 
                          width=edge_widths, alpha=0.6, 
                          arrows=G.is_directed(), arrowsize=20, ax=ax)
    
    # Draw labels
    if show_labels:
        # Truncate long labels
        labels = {node: node[:20] + '...' if len(str(node)) > 20 else str(node) 
                 for node in subgraph.nodes()}
        nx.draw_networkx_labels(subgraph, pos, labels=labels, font_size=8, 
                               font_weight='bold', ax=ax)
    
    # Draw edge labels
    if show_edge_labels and subgraph.number_of_edges() < 50:  # Only show if not too crowded
        edge_labels = {(u, v): data.get('edge_type', '')[:10] 
                      for u, v, data in subgraph.edges(data=True)}
        nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels, 
                                    font_size=6, font_color='red', ax=ax)
    
    # Add title and stats
    stats_text = (f"Method: {method} | Nodes: {subgraph.number_of_nodes()} | "
                 f"Edges: {subgraph.number_of_edges()} | Center: {center_node}")
    plt.title(f"Knowledge Graph Subgraph\n{stats_text}", fontsize=14, fontweight='bold')
    
    # Add legend
    node_types, edge_types = create_legend(G, subgraph)
    
    # Node type legend
    color_map = {
        'Protein': '#FF6B6B', 'GO': '#4ECDC4', 'KEGG': '#45B7D1', 
        'Pfam': '#96CEB4', 'PDB': '#FECA57', 'STRING': '#FF9FF3',
        'UniProt': '#54A0FF', 'NCBI': '#5F27CD', 'DeepGO': '#00D2D3',
        'Embedding': '#FF5722', 'AlphaFold': '#795548'
    }
    
    legend_elements = []
    for node_type in sorted(set(node_types)):
        color = '#95A5A6'  # Default
        for key, val in color_map.items():
            if key in node_type:
                color = val
                break
        legend_elements.append(mpatches.Patch(color=color, label=node_type[:15]))
    
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), 
                 title="Node Types", fontsize=8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {save_path}")
    
    plt.show()
    
    return subgraph

def visualize_knowledge_graph_degrees(G, node_n=0, n_degree=2):
    '''Legacy function - calls advanced visualization with degree method.'''
    if G.number_of_nodes() == 0:
        print("The graph is empty. No nodes to visualize.")
        return
    
    center_node = list(G.nodes())[node_n % G.number_of_nodes()]
    return visualize_subgraph_advanced(G, center_node=center_node, method='degree', 
                                     radius=n_degree, figsize=(12, 12))

def visualize_knowledge_graph_proximity(G, node_n=0, max_protein_nodes=5):
    '''Legacy function - calls advanced visualization with proximity method.'''
    if G.number_of_nodes() == 0:
        print("The graph is empty. No nodes to visualize.")
        return
        
    center_node = list(G.nodes())[node_n % G.number_of_nodes()]
    return visualize_subgraph_advanced(G, center_node=center_node, method='proximity', 
                                     max_nodes=max_protein_nodes*10, figsize=(12, 12))

def analyze_subgraph_statistics(G, subgraph):
    """Print comprehensive statistics about the subgraph."""
    print("\n=== SUBGRAPH ANALYSIS ===")
    print(f"Nodes: {subgraph.number_of_nodes()}")
    print(f"Edges: {subgraph.number_of_edges()}")
    print(f"Density: {nx.density(subgraph):.4f}")
    
    if subgraph.number_of_nodes() > 1:
        print(f"Average clustering coefficient: {nx.average_clustering(subgraph):.4f}")
        
        if nx.is_connected(subgraph):
            print(f"Average shortest path length: {nx.average_shortest_path_length(subgraph):.4f}")
            print(f"Diameter: {nx.diameter(subgraph)}")
        else:
            print("Graph is not connected")
            print(f"Number of connected components: {nx.number_connected_components(subgraph)}")
    
    # Node type distribution
    node_types = Counter()
    for node in subgraph.nodes():
        node_attr = G.nodes[node].get('node_attr', 'Unknown')
        node_types[node_attr] += 1
    
    print("\nNode type distribution:")
    for node_type, count in node_types.most_common():
        print(f"  {node_type}: {count}")
    
    # Edge type distribution
    edge_types = Counter()
    for _, _, data in subgraph.edges(data=True):
        edge_type = data.get('edge_type', 'unknown')
        edge_types[edge_type] += 1
    
    print("\nEdge type distribution:")
    for edge_type, count in edge_types.most_common():
        print(f"  {edge_type}: {count}")
    
    # Top degree nodes
    degrees = dict(subgraph.degree())
    top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:5]
    print("\nTop 5 nodes by degree:")
    for node, degree in top_nodes:
        node_type = G.nodes[node].get('node_attr', 'Unknown')
        print(f"  {node} ({node_type}): degree {degree}")
    
    print("=" * 25)

def demonstrate_visualization_modes(G, sample_node=None):
    """
    Demonstrate different visualization modes available in the advanced subgraph function.
    """
    if G.number_of_nodes() == 0:
        print("Graph is empty - cannot demonstrate visualizations.")
        return
    
    if sample_node is None:
        sample_node = max(G.nodes(), key=lambda x: G.degree(x))
    
    print(f"\n=== DEMONSTRATION: Different Visualization Modes ===")
    print(f"Using center node: {sample_node}")
    
    # 1. Degree-based visualization
    print("\n1. Degree-based subgraph (2-hop neighborhood):")
    subgraph1 = visualize_subgraph_advanced(
        G, center_node=sample_node, method='degree', radius=2, 
        max_nodes=30, layout='spring', figsize=(14, 10)
    )
    
    # 2. Proximity-based visualization  
    print("\n2. Proximity-based subgraph (closest nodes):")
    subgraph2 = visualize_subgraph_advanced(
        G, center_node=sample_node, method='proximity', radius=3,
        max_nodes=25, layout='kamada_kawai', figsize=(14, 10)
    )
    
    # 3. Centrality-based visualization
    print("\n3. Centrality-based subgraph (high centrality nodes):")
    subgraph3 = visualize_subgraph_advanced(
        G, center_node=sample_node, method='centrality', 
        max_nodes=20, layout='circular', figsize=(12, 12)
    )
    
    # 4. Filtered visualization (protein nodes only)
    print("\n4. Protein-focused subgraph (proteins and direct connections):")
    def protein_filter(node, attrs):
        return 'Protein' in attrs.get('node_attr', '') or any('Protein' in neighbor_attr.get('node_attr', '') 
                                                               for neighbor_attr in [G.nodes[n] for n in G.neighbors(node)])
    
    subgraph4 = visualize_subgraph_advanced(
        G, center_node=sample_node, method='degree', radius=1,
        max_nodes=15, layout='spring', figsize=(12, 10),
        node_filter=protein_filter, show_edge_labels=False
    )
    
    # Analyze the last subgraph
    if subgraph4:
        analyze_subgraph_statistics(G, subgraph4)  



def save_knowledge_graph(G, output_file):
    '''Saves the knowledge graph to a file in GraphML format.'''
    nx.write_graphml(G, output_file)
    print(f"Knowledge graph saved to {output_file}")

if __name__ == "__main__":
    # CLI
    parser = argparse.ArgumentParser(description="Build a knowledge graph from protein info TSV file or folder.")
    parser.add_argument("file_path", type=str, help="Path to the protein info TSV file or folder containing it.")
    parser.add_argument("--graph_input", type=str, help="Path to a GraphML file to load an existing graph instead of building a new one.")
    parser.add_argument("--selection_percentage", type=int, default=10, help="Percentage of proteins to select randomly.")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed for selection.")
    parser.add_argument('--demo_visualizations', action='store_true', 
                       help="Run demonstration of different visualization modes.")
    parser.add_argument('--custom_visualization', action='store_true',
                       help="Run custom visualization with user-specified parameters.")
    parser.add_argument('--center_node', help="Specify center node for visualization.")
    parser.add_argument('--vis_method', choices=['degree', 'proximity', 'centrality', 'community'],
                       default='degree', help="Visualization method to use.")
    parser.add_argument('--vis_radius', type=int, default=2, help="Radius for degree/proximity methods.")
    parser.add_argument('--vis_max_nodes', type=int, default=30, help="Maximum nodes in visualization.")
    parser.add_argument('--vis_layout', choices=['spring', 'circular', 'kamada_kawai', 'spectral'],
                       default='spring', help="Layout algorithm for visualization.")
    parser.add_argument('--save_viz', help="Path to save visualization image.")
    args = parser.parse_args()

    # Example usage
    file_path = args.file_path
    graph_input = args.graph_input
    selection_percentage = args.selection_percentage
    random_seed = args.random_seed

    print(f"Selection Percentage set at: {selection_percentage}%")
    
    list_of_targets = [
        ('DeepGO_MF', 'has_molecular_function', None),
        ('DeepGO_CC', 'located_in_cellular_component', None),
        ('DeepGO_BP', 'involved_in_biological_process', None),
        ('Pfam_domains', 'has_pfam_domain', 1),
        ('KEGG_pathways', 'in_kegg_pathway', 1),
        ('PROSITE_annotations', 'has_prosite_annotation', 1),
        ('PPI_target', 'interacts_with', None),
        ('embedding_similarity', 'embeddings_similar_to', None)
    ]
    list_of_node_attr = [
        ('UniProt_ID', 'Protein'),
        ('DeepGO_MF', 'Molecular_Function'),
        ('DeepGO_CC', 'Cellular_Component'),
        ('DeepGO_BP', 'Biological_Process'),
        ('Pfam_domains', 'Pfam_Domain'),
        ('KEGG_pathways', 'KEGG_Pathway'),
        ('PROSITE_annotations', 'PROSITE_Annotation')
    ]
    list_of_node_alias_attr = [
        ('UniProt_ID', 'NCBI_ID')
        # ('UniProt_ID', 'Organism'),
        # ('UniProt_ID', 'Gene_Name')
    ]

    # Load existing graph or build new graph with error handling
    if graph_input and os.path.exists(graph_input):
        try:
            # Check if file is not empty
            if os.path.getsize(graph_input) > 0:
                print(f"Loading existing graph from {graph_input}...")
                G = nx.read_graphml(graph_input)
            else:
                print(f"Graph file {graph_input} is empty. Creating new graph...")
                G = nx.Graph()
        except Exception as e:
            print(f"Error loading graph file {graph_input}: {e}")
            print("Creating new graph...")
            # load all tsv files from folder
            df_list = load_tsvs_from_folder(file_path, selection_percentage=selection_percentage) if os.path.isdir(file_path) else [load_tsv_into_df(file_path, selection_percentage=selection_percentage)]
            print(f"Loaded {len(df_list)} TSV files from {file_path}.")

            # go through the list of tsvs/dataframes and build triples and node attribute dataframes
            triples_df_list = build_kg_target_triples_df_list(df_list, list_of_targets)
            node_attr_df_list = build_kg_node_attr_df_list(df_list, list_of_node_attr)
            node_alias_attr_df_list = build_kg_node_alias_attr_df_list(df_list, list_of_node_alias_attr)

            # build knowledge graph
            G = build_knowledge_graph(triples_df_list, node_attr_df_list, node_alias_attr_df_list)
    else:
        # load all tsv files from folder
        df_list = load_tsvs_from_folder(file_path, selection_percentage=selection_percentage) if os.path.isdir(file_path) else [load_tsv_into_df(file_path, selection_percentage=selection_percentage)]
        print(f"Loaded {len(df_list)} TSV files from {file_path}.")

        # go through the list of tsvs/dataframes and build triples and node attribute dataframes
        triples_df_list = build_kg_target_triples_df_list(df_list, list_of_targets)
        node_attr_df_list = build_kg_node_attr_df_list(df_list, list_of_node_attr)
        node_alias_attr_df_list = build_kg_node_alias_attr_df_list(df_list, list_of_node_alias_attr)

        # build knowledge graph
        G = build_knowledge_graph(triples_df_list, node_attr_df_list, node_alias_attr_df_list)

    print(f"Knowledge graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    
    # Handle different visualization modes
    if args.demo_visualizations:
        demonstrate_visualization_modes(G)
    elif args.custom_visualization:
        # Custom visualization with user parameters
        subgraph = visualize_subgraph_advanced(
            G, center_node=args.center_node, method=args.vis_method, 
            radius=args.vis_radius, max_nodes=args.vis_max_nodes,
            layout=args.vis_layout, save_path=args.save_viz
        )
        if subgraph:
            analyze_subgraph_statistics(G, subgraph)
    else:
        # Default: simple degree visualization
        visualize_knowledge_graph_degrees(G, node_n=0, n_degree=2)

        # print a specific node's attributes and connections for verification
        if G.number_of_nodes() > 0:
            sample_node = list(G.nodes())[0]
            print(f"Sample node: {sample_node}")
            print(f"Attributes: {G.nodes[sample_node]}")
            print(f"Connections: {list(G.edges(sample_node, data=True))}")

    # save knowledge graph
    output_file = "knowledge_graph.graphml"
    save_knowledge_graph(G, output_file)
    



    
    