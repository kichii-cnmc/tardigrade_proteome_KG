from data_processor import (
    load_protein_info,
    randomize_protein_info_selection,
    build_set_of_kg_triples,
    build_set_of_node_types
)
import pandas as pd
import networkx as nx

def build_triples_and_nodes_from_protein_info(file_path, 
                              list_of_targets, 
                              list_of_node_types, 
                              selection_percentage=100, 
                              random_seed=42):
    """Builds knowledge graph triples and node type DataFrames from protein info TSV file."""
    # Load protein info
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

if __name__ == "__main__":
    # Example usage
    file_path = "data/3_organized/protein_info_RV.tsv"
    
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
        selection_percentage=10, 
        random_seed=42
    )
    
    kg = build_knowledge_graph(triples_dfs, nodes_dfs)
    
    print(f"Knowledge graph has {kg.number_of_nodes()} nodes and {kg.number_of_edges()} edges.")