from data_processor import (
    load_protein_info,
    randomize_protein_info_selection,
    build_set_of_kg_triples,
    build_set_of_node_types
)
import pandas as pd
import networkx as nx
import argparse
import os
import glob
import matplotlib.pyplot as plt

def build_triples_and_nodes_from_protein_info(file_path, 
                              list_of_targets, 
                              list_of_node_types, 
                              selection_percentage=100, 
                              random_seed=42):
    """Builds knowledge graph triples and node type DataFrames from protein info TSV file or all files in a folder."""

    # Load protein info
    if os.path.isdir(file_path):
        all_files = glob.glob(os.path.join(file_path, "*.tsv"))
        df_list = [load_protein_info(f) for f in all_files]
        df = pd.concat(df_list, ignore_index=True)
    else:
        df = load_protein_info(file_path)
    
    # Randomize selection if needed
    df = randomize_protein_info_selection(df, percentage=selection_percentage, seed=random_seed)
    
    # Build triples DataFrames
    triples_dfs = build_set_of_kg_triples(df, list_of_targets)
    
    # Build node type DataFrames
    nodes_dfs = build_set_of_node_types(df, list_of_node_types)
    
    return triples_dfs, nodes_dfs

def build_knowledge_graph(triples_dfs, nodes_dfs):
    """Builds a NetworkX knowledge graph from triples and node type DataFrames."""
    G = nx.DiGraph()
    
    # Add nodes with types
    for nodes_df in nodes_dfs:
        for _, row in nodes_df.iterrows():
            G.add_node(row['id'], node_type=row['node_type'])
    
    # Add edges from triples
    for triples_df in triples_dfs:
        for _, row in triples_df.iterrows():
            G.add_edge(row['source'], row['target'], edge_type=row['edge_type'], weight=row['weight'])
    
    return G

def visualize_knowledge_graph(G, n_degree=2):
    '''Visualizes a subgraph of KG to the extent of n_degree from a random node.'''
    if G.number_of_nodes() == 0:
        print("The graph is empty. No nodes to visualize.")
        return
    random_node = list(G.nodes())[0]
    nodes_to_include = set([random_node])
    for _ in range(n_degree):
        neighbors = set()
        for node in nodes_to_include:
            neighbors.update(G.neighbors(node))
            neighbors.update(G.predecessors(node))
        nodes_to_include.update(neighbors)
    subgraph = G.subgraph(nodes_to_include)
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(subgraph, seed=42)
    nx.draw(subgraph, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=8, font_weight='bold', arrows=True)
    edge_labels = nx.get_edge_attributes(subgraph, 'edge_type')
    nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels, font_color='red', font_size=6)
    plt.title(f"Subgraph of Knowledge Graph (n_degree={n_degree})")
    plt.show()

def save_knowledge_graph(G, output_path):
    '''Saves the knowledge graph to a GraphML file.'''
    nx.write_graphml(G, output_path)
    print(f"Knowledge graph saved to {output_path}")

if __name__ == "__main__":
    # CLI
    parser = argparse.ArgumentParser(description="Build a knowledge graph from protein info TSV file or folder.")
    parser.add_argument("file_path", type=str, help="Path to the protein info TSV file or folder containing it.")
    parser.add_argument("--selection_percentage", type=int, default=10, help="Percentage of proteins to select randomly.")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed for selection.")
    args = parser.parse_args()

    # Example usage
    file_path = args.file_path
    selection_percentage = args.selection_percentage
    random_seed = args.random_seed
    
    list_of_targets = [
        ('GO_mf', 'has_molecular_function', 1),
        ('GO_cc', 'located_in_cellular_component', 1),
        ('GO_bp', 'involved_in_biological_process', 1),
        ('Pfam_domains', 'has_pfam_domain', 1),
        ('KEGG_pathways', 'in_kegg_pathway', 1),
        ('PROSITE_annotations', 'has_prosite_annotation', 1)
    ]
    list_of_node_types = [
        ('UniProt_ID', 'Protein'),
        ('GO_mf', 'Molecular_Function'),
        ('GO_cc', 'Cellular_Component'),
        ('GO_bp', 'Biological_Process'),
        ('Pfam_domains', 'Pfam_Domain'),
        ('KEGG_pathways', 'KEGG_Pathway'),
        ('PROSITE_annotations', 'PROSITE_Annotation')
    ]
    
    triples_dfs, nodes_dfs = build_triples_and_nodes_from_protein_info(
        file_path, 
        list_of_targets, 
        list_of_node_types, 
        selection_percentage=selection_percentage, 
        random_seed=random_seed
    )
    
    kg = build_knowledge_graph(triples_dfs, nodes_dfs)
    
    print(f"Knowledge graph has {kg.number_of_nodes()} nodes and {kg.number_of_edges()} edges.")

    # print random sample of nodes and edges (n = 5)
    print("\nSample nodes:")
    for node in list(kg.nodes(data=True))[:5]:
        print(node)
    print("\nSample edges:")
    for edge in list(kg.edges(data=True))[:5]:
        print(edge)

    # Visualize a small subgraph
    visualize_knowledge_graph(kg, n_degree=1)

    # Save the knowledge graph
    save_knowledge_graph(kg, "knowledge_graph.graphml")