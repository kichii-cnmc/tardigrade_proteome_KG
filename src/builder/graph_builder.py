# Used to build knowledge graph from protein info TSV files or folder.

import pandas as pd
import networkx as nx
import argparse
import os
import glob
import matplotlib.pyplot as plt

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
    target_names = [target[0] for target in list_of_targets]
    df_names = collect_header_names(df_list, 1)
    for i, df in enumerate(df_list):
        if df_names[i] in target_names:
            target_index = target_names.index(df_names[i])
            target, edge_type, weight = list_of_targets[target_index]
            triples_df = build_kg_target_triples_df(df, target, edge_type, weight)
            all_triples_df.append(triples_df)
    return all_triples_df

def build_kg_node_attr_df_list(df_list, list_of_node_attr):
    '''List version of build_kg_node_attr_df, used to process multiple node attributes.'''
    all_nodes_df = []
    node_attr_names = [node_attr[0] for node_attr in list_of_node_attr]
    df_names = collect_header_names(df_list, 1)
    for i, df in enumerate(df_list):
        if df_names[i] in node_attr_names:
            node_attr_index = node_attr_names.index(df_names[i])
            id_column, node_attr_value = list_of_node_attr[node_attr_index]
            nodes_df = build_kg_node_attr_df(df, id_column, node_attr_value)
            all_nodes_df.append(nodes_df)
    return all_nodes_df

def build_kg_node_alias_attr_df_list(df_list, list_of_node_alias_attr):
    '''List version of build_kg_node_alias_attr_df, used to process multiple alias attributes.'''
    all_alias_df = []
    alias_attr_names = [alias_attr[0] for alias_attr in list_of_node_alias_attr]
    df_names = collect_header_names(df_list, 1)
    for i, df in enumerate(df_list):
        if df_names[i] in alias_attr_names:
            alias_attr_index = alias_attr_names.index(df_names[i])
            id_column, alias_attr_column = list_of_node_alias_attr[alias_attr_index]
            alias_df = build_kg_node_alias_attr_df(df, id_column, alias_attr_column)
            all_alias_df.append(alias_df)
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
    
    # Add triples as edges
    for triples_df in triples_df_list:
        for index, row in triples_df.iterrows():
            G.add_edge(row['source'], row['target'], edge_type=row['edge_type'], weight=row['weight'])
    
    # Add node attributes
    for node_attr_df in node_attr_df_list:
        for index, row in node_attr_df.iterrows():
            if 'node_attr' not in G.nodes[row['id']]:
                G.nodes[row['id']]['node_attr'] = []
            G.nodes[row['id']]['node_attr'].append(row['node_attr'])
    
    # Add alias attributes
    for alias_attr_df in node_alias_attr_df_list:
        for index, row in alias_attr_df.iterrows():
            if 'alias_attr' not in G.nodes[row['id']]:
                G.nodes[row['id']]['alias_attr'] = []
            G.nodes[row['id']]['alias_attr'].append(row['alias_attr'])
    
    return G

def visualize_knowledge_graph(G, node_n = 0, n_degree=2):
    '''Visualizes a subgraph of KG to the extent of n_degree from a random node.'''
    if G.number_of_nodes() == 0:
        print("The graph is empty. No nodes to visualize.")
        return
    random_node = list(G.nodes())[node_n % G.number_of_nodes()]
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

    print(f"Selection Percentage set at: {selection_percentage}%")
    
    list_of_targets = [
        ('GO_mf', 'has_molecular_function', 1),
        ('GO_cc', 'located_in_cellular_component', 1),
        ('GO_bp', 'involved_in_biological_process', 1),
        ('Pfam_domains', 'has_pfam_domain', 1),
        ('KEGG_pathways', 'in_kegg_pathway', 1),
        ('PROSITE_annotations', 'has_prosite_annotation', 1),
        ('PPI_target', 'interacts_with', None)
    ]
    list_of_node_attr = [
        ('UniProt_ID', 'Protein'),
        ('GO_mf', 'Molecular_Function'),
        ('GO_cc', 'Cellular_Component'),
        ('GO_bp', 'Biological_Process'),
        ('Pfam_domains', 'Pfam_Domain'),
        ('KEGG_pathways', 'KEGG_Pathway'),
        ('PROSITE_annotations', 'PROSITE_Annotation')
    ]
    list_of_node_alias_attr = [
        ('UniProt_ID', 'NCBI_ID')
        # ('UniProt_ID', 'Organism'),
        # ('UniProt_ID', 'Gene_Name')
    ]

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
    # visualize knowledge graph 1
    visualize_knowledge_graph(G, node_n = 0, n_degree=2)
    # visualize starting from another node
    visualize_knowledge_graph(G, node_n = 1, n_degree=2)

    # print a specific node's attributes and connections for verification: A0A1D1W4Z0
    sample_node = list(G.nodes())[0]
    print(f"Sample node: {sample_node}")
    print(f"Attributes: {G.nodes[sample_node]}")
    print(f"Connections: {list(G.edges(sample_node, data=True))}")
    



    
    