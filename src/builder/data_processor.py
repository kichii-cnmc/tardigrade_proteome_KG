# helper functions used to build knowledge graph from protein info tsv files
import pandas as pd

def load_protein_info(file_path):
    """Load protein info from TSV file into a pandas DataFrame."""
    df = pd.read_csv(file_path, sep="\t")
    return df

def randomize_protein_info_selection(df, percentage = 100, seed = 42, orphan_removal = False):
    '''Randomly shuffle the DataFrame rows and select a percentage of rows.'''
    if percentage < 100 and not orphan_removal:
        df = df.sample(frac=percentage/100, random_state=seed).reset_index(drop=True)
    if orphan_removal:
        # remove rows where all target columns are NaN or empty
        original_length = len(df)
        target_columns = ['GO_mf', 'GO_cc', 'GO_bp', 'Pfam_domains', 'KEGG_pathways', 'PROSITE_annotations']
        df = df.dropna(subset=target_columns, how='all')
        for col in target_columns:
            df = df[~df[col].astype(str).str.strip().eq('')]
        df = df.reset_index(drop=True)
        if percentage < 100:
            percentage/((len(df)/original_length)*100)
            df = df.sample(frac=percentage/((len(df)/original_length)*100), random_state=seed).reset_index(drop=True)
    return df

def build_kg_triple_dataframe(df, target, edge_type = None, weight = 1):
    '''Selects a target column from the input DataFrame and builds a triples DataFrame.'''
    triples_list = []
    for index, row in df.iterrows():
        protein_id = row.iloc[0] # uses first column for ID
        target_values = str(row[target]).split(';') if pd.notna(row[target]) else []
        for target_value in target_values:
            target_value = target_value.strip()
            if target_value:
                triples_list.append({
                    'source': protein_id,
                    'target': target_value,
                    'edge_type': edge_type,
                    'weight': weight
                })
    triples_df = pd.DataFrame(triples_list)
    return triples_df

def build_set_of_kg_triples(df, list_of_targets):
    '''Builds a list of triples DataFrames based on a list of tuples holding target columns, edges, and weights.'''
    all_triples_df = []
    for target, edge_type, weight in list_of_targets:
        triples_df = build_kg_triple_dataframe(df, target, edge_type, weight)
        all_triples_df.append(triples_df)
    return all_triples_df

def build_node_type_dataframe(df, id_column, node_type):
    '''Builds a DataFrame for nodes of a specific type, separating by semicolons if needed.'''
    nodes_set = set()
    for index, row in df.iterrows():
        id_values = str(row[id_column]).split(';') if pd.notna(row[id_column]) else []
        for id_value in id_values:
            id_value = id_value.strip()
            if id_value:
                nodes_set.add((id_value, node_type))
    nodes_df = pd.DataFrame(list(nodes_set), columns=['id', 'node_type'])
    return nodes_df

def build_set_of_node_types(df, list_of_node_types):
    '''Builds a list of node type DataFrames based on a list of tuples holding id columns and node types.'''
    all_nodes_df = []
    for id_column, node_type in list_of_node_types:
        nodes_df = build_node_type_dataframe(df, id_column, node_type)
        all_nodes_df.append(nodes_df)
    return all_nodes_df

if __name__ == "__main__":
    print("Warning: This module is intended to be imported and used in other scripts.\n")

    # Example usage
    file_path = "data/3_organized/protein_info_RV.tsv"
    df = load_protein_info(file_path)
    df = randomize_protein_info_selection(df, percentage=10, seed=42, orphan_removal=True)
    print(f"Selected {len(df)} proteins after randomization and orphan removal.")

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

    triples_df = build_set_of_kg_triples(df, list_of_targets)
    for i, triples in enumerate(triples_df):
        print(f"\nTriples DataFrame {i+1}:\n", triples.head())
        print(f"Number of triples: {len(triples)}")

    nodes_df = build_set_of_node_types(df, list_of_node_types)
    for i, nodes in enumerate(nodes_df):
        print(f"\nNodes DataFrame {i+1}:\n", nodes.head())
        print(f"Number of nodes: {len(nodes)}")