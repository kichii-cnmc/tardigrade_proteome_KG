# script that uses igraph to build and manage the graph of the project
import igraph as ig
import os
import json
import pandas as pd
import glob
import argparse
import matplotlib
matplotlib.use('Agg')  # non-interactive backend — must be set before pyplot import
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

class AliasManager:
    def __init__(self):
        self.alias_to_name = {}
        self.names_to_aliases = {}
    
    def add_alias(self, alias, node_name):
        # add alias to node_name mapping
        if alias in self.alias_to_name:
            if self.alias_to_name[alias] != node_name:
                print(f"Warning: Alias '{alias}' already mapped to '{self.alias_to_name[alias]}', cannot remap to '{node_name}'.")
        else:
            self.alias_to_name[alias] = node_name
            self.names_to_aliases.setdefault(node_name, set()).add(alias)
    
    def get_node_name(self, alias):
        return self.alias_to_name.get(alias, None)
    
    def get_aliases(self, node_name):
        return self.names_to_aliases.get(node_name, set())

class IGraphBuilder:
    DIRECTED_GRAPH = True
    GRAPH_EDGE_LABELS = { # column name: (edge label, weight)
        'DeepGO_MF': ('has_molecular_function', None),
        'DeepGO_CC': ('located_in_cellular_component', None),
        'DeepGO_BP': ('involved_in_biological_process', None),
        'Pfam_domains': ('has_pfam_domain', 1),
        'KEGG_pathways': ('in_kegg_pathway', 1),
        'PROSITE_annotations': ('has_prosite_annotation', 1),
        'PPI_target': ('interacts_with', None),
        'protein2': ('embeddings_similar_to', None)
    }
    GRAPH_NODE_LABELS = { # column name: node label
        'UniProt_ID': 'Protein',
        'DeepGO_MF': 'Molecular_Function',
        'DeepGO_CC': 'Cellular_Component',
        'DeepGO_BP': 'Biological_Process',
        'Pfam_domains': 'Pfam_Domain',
        'KEGG_pathways': 'KEGG_Pathway',
        'PROSITE_annotations': 'PROSITE_Annotation',
        'PPI_target': 'Protein',
        'protein2': 'Protein'
    }
    GRAPH_NODE_ALIAS_TYPES = ['NCBI_ID', 'Gene_Name', 'AF_structures'] # columns that can be added as node attributes
    ALIAS_MAPPER = AliasManager()

    def __init__(self):
        self.graph = ig.Graph(directed=self.DIRECTED_GRAPH)

    def add_node(self, node_id, **attributes):
        if not self.graph.vs.select(name=node_id):
            self.graph.add_vertex(name=node_id, **attributes)

    def add_edge(self, source_id, target_id, **attributes):
        self.add_node(source_id)
        self.add_node(target_id)
        if not self.graph.es.select(_source=self.graph.vs.find(name=source_id).index, _target=self.graph.vs.find(name=target_id).index):
            self.graph.add_edge(source_id, target_id, **attributes)

    def map_names_to_graph_ids(self, name_to_id):
        '''Maps node names to their corresponding graph IDs.'''
        for vertex in self.graph.vs:
            name_to_id[vertex['name']] = vertex.index
    
    def add_nodes_from_list(self, nodes_list, node_type=None):
        '''Add nodes from a list of node names. New nodes are inserted with node_type;
        existing nodes have their node_type updated if node_type is provided.'''
        existing_names = set(self.graph.vs['name']) if self.graph.vcount() > 0 else set()
        seen = set()
        new_nodes, update_nodes = [], []
        for node in nodes_list:
            if node not in seen:
                seen.add(node)
                (update_nodes if node in existing_names else new_nodes).append(node)
        if new_nodes:
            self.graph.add_vertices(
                len(new_nodes),
                attributes={'name': new_nodes, 'node_type': [node_type] * len(new_nodes)}
            )
        if node_type is not None and update_nodes:
            update_set = set(update_nodes)
            vs = self.graph.vs.select(lambda v: v['name'] in update_set)
            vs['node_type'] = [node_type] * len(vs)

    def add_edges_from_list(self, edges_list, edge_type=None):
        '''Add edges from a list of (source, target, weight) tuples.'''
        if not edges_list:
            return

        # Ensure all nodes exist in one bulk call
        all_nodes = [node for source, target, _ in edges_list for node in (source, target)]
        self.add_nodes_from_list(all_nodes)

        # Build name -> vertex index mapping in one pass
        name_to_idx = dict(zip(self.graph.vs['name'], self.graph.vs.indices))

        # Get existing edges as a set of (src_idx, tgt_idx) for O(1) lookup
        existing_edges = set(self.graph.get_edgelist())

        # Filter to new unique edges only
        new_edges = []
        weights = []
        for source, target, weight in edges_list:
            key = (name_to_idx[source], name_to_idx[target])
            if key not in existing_edges:
                new_edges.append(key)
                weights.append(weight if weight is not None else 1)
                existing_edges.add(key)  # prevent duplicates within edges_list

        if new_edges:
            self.graph.add_edges(
                new_edges,
                attributes={'edge_type': [edge_type] * len(new_edges), 'weight': weights}
            )

    def add_alias_attributes_from_list(self, alias_list):
        '''Add alias list to mapping object. Alias list = (alias_value, node_name)'''
        for alias_value, node_name in alias_list:
            self.ALIAS_MAPPER.add_alias(alias_value, node_name)

    def build_go_term_alias_mapping(self, node_names_list):
        '''Build mapping of GO term IDs to their names for better interpretability in the graph.'''
        alias_list = []
        for name in node_names_list:
            if name.startswith("GO:"):
                go_id = name[:10]
                alias_list.append((go_id, name))
        return alias_list

    def add_to_graph_tsv(self, tsv_file_path):
        '''Add nodes and edges from a TSV file, determines the correct format based on the 2nd column name.'''
        if not os.path.exists(tsv_file_path):
            print(f"File {tsv_file_path} does not exist.")
            return
        df = pd.read_csv(tsv_file_path, sep='\t')
        if df.empty:
            print(f"File {tsv_file_path} is empty.")
            return
        second_col = df.columns[1]
        # Handle pandas automatic column renaming for duplicates (e.g., 'UniProt_ID.1')
        second_col_clean = second_col.split('.')[0] if '.' in second_col else second_col
        # if second column matches an edge label, add edges
        if second_col_clean in self.GRAPH_EDGE_LABELS:
            edge_label, weight_col = self.GRAPH_EDGE_LABELS[second_col]
            edges_list = []
            for _, row in df.iterrows():
                source = row.iloc[0]
                target = row.iloc[1]
                weight = round(row.iloc[2], 4) if weight_col is None else 1
                edges_list.append((source, target, weight))
            print(f"Adding edges of type '{edge_label}' from {tsv_file_path}...")
            # check column name for node type mapping, default to 'UnknownType' if not found
            node_type = self.GRAPH_NODE_LABELS.get(second_col_clean, "UnknownType")
            if node_type == "UnknownType":
                print(f"    Warning: No node type mapping found for column '{second_col_clean}', defaulting to 'UnknownType'.")
            self.add_nodes_from_list([target for _, target, _ in edges_list], node_type=node_type)
            if node_type in ["Molecular_Function", "Cellular_Component", "Biological_Process"]:
                # inputs a non-tuple list of the GO terms added to the graph as nodes to build the alias mapping for better interpretability
                alias_list = self.build_go_term_alias_mapping(list(set(target for _, target, _ in edges_list)))
                self.add_alias_attributes_from_list(alias_list)
            self.add_edges_from_list(edges_list, edge_type=edge_label)
        # if second column matches a node label, add nodes
        elif second_col_clean in self.GRAPH_NODE_LABELS:
            node_type = self.GRAPH_NODE_LABELS[second_col_clean]
            print(f"Adding nodes of type '{node_type}' from {tsv_file_path}...")
            self.add_nodes_from_list(df.iloc[:, 1].tolist(), node_type=node_type)
        # if second column matches an alias type, add to alias mapping
        elif second_col_clean in self.GRAPH_NODE_ALIAS_TYPES:
            print(f"Adding alias attributes from {tsv_file_path}...")
            alias_list = list(zip(df.iloc[:, 1].tolist(), df.iloc[:, 0].tolist())) # swap to (alias_value, node_name) format
            self.add_alias_attributes_from_list(alias_list)
        else:
            print(f"Unrecognized format in {tsv_file_path}, skipping.")

    def save_graph(self, file_path):
        with open(file_path, 'w') as f:
            json.dump(self.graph.to_dict(), f)

    def load_graph(self, file_path):
        with open(file_path, 'r') as f:
            graph_dict = json.load(f)
            self.graph = ig.Graph.from_dict(graph_dict)

    def evaluate_graph(self):
        '''Evaluate the graph by printing basic statistics.'''
        print(f"Number of nodes: {self.graph.vcount()}")
        print(f"Number of edges: {self.graph.ecount()}")
        print(f"Graph density: {self.graph.density():.4f}")
        print(f"Average degree: {sum(self.graph.degree()) / self.graph.vcount():.2f}")
        print(f"Connected components: {len(self.graph.components())}")
        print()
        print("Example Graph Triples:")
        for edge in self.graph.es[:5]:
            source = self.graph.vs[edge.source]["name"]
            target = self.graph.vs[edge.target]["name"]
            weight = edge['weight']
            print(f"({source}, {edge['edge_type']}, {target}, weight={weight})")
        print("\nExample Nodes:")
        for v in self.graph.vs[:5]:
            print(f"({v['name']}, {v['node_type']})")

    def get_graph(self):
        return self.graph
    
    def sample_representative_subgraph(self, k=6, max_edges=600):
        '''Return a representative subgraph of the top-k highest-degree nodes per node_type.
        Node selection is based on inter-type degree (edges to nodes of a different type)
        so the chosen nodes are structurally meaningful across the KG.
        Subgraph construction uses neighbor traversal to avoid scanning all 3M edges.'''
        import random
        import numpy as np

        # Compute inter-type degree via numpy (one vectorised pass over all edges)
        edge_arr = np.array(self.graph.get_edgelist())  # (E, 2) int array
        node_types = np.array(self.graph.vs['node_type'], dtype=object)
        if edge_arr.size > 0:
            inter_mask = node_types[edge_arr[:, 0]] != node_types[edge_arr[:, 1]]
            inter_degrees = np.bincount(
                edge_arr[inter_mask].ravel(), minlength=self.graph.vcount()
            ).tolist()
        else:
            inter_degrees = [0] * self.graph.vcount()

        # Group vertices by type and select top-k by inter-type degree
        type_to_indices = {}
        for v in self.graph.vs:
            type_to_indices.setdefault(v['node_type'], []).append(v.index)
        selected = []
        for indices in type_to_indices.values():
            top_k = sorted(indices, key=lambda i: inter_degrees[i], reverse=True)[:k]
            selected.extend(top_k)

        selected_set = set(selected)
        old_to_new = {old: new for new, old in enumerate(selected)}

        # Build fresh subgraph via neighbor traversal — never copies the full edge list
        sub = ig.Graph(directed=self.DIRECTED_GRAPH)
        sub.add_vertices(
            len(selected),
            attributes={
                'name':      [self.graph.vs[i]['name']      for i in selected],
                'node_type': [self.graph.vs[i]['node_type'] for i in selected],
            }
        )

        new_edges, etypes, weights = [], [], []
        seen = set()
        for old_src in selected:
            for old_tgt in self.graph.successors(old_src):
                if old_tgt in selected_set:
                    key = (old_to_new[old_src], old_to_new[old_tgt])
                    if key not in seen:
                        seen.add(key)
                        eid = self.graph.get_eid(old_src, old_tgt)
                        new_edges.append(key)
                        etypes.append(self.graph.es[eid]['edge_type'])
                        weights.append(self.graph.es[eid]['weight'])

        if len(new_edges) > max_edges:
            keep = sorted(random.sample(range(len(new_edges)), max_edges))
            new_edges = [new_edges[i] for i in keep]
            etypes    = [etypes[i]    for i in keep]
            weights   = [weights[i]   for i in keep]

        if new_edges:
            sub.add_edges(new_edges, attributes={'edge_type': etypes, 'weight': weights})
        return sub

    def visualize_graph(self, g, output_path):
        NODE_TYPE_COLORS = {
            'Protein':               '#4C72B0',
            'Molecular_Function':    '#DD8452',
            'Cellular_Component':    '#55A868',
            'Biological_Process':    '#C44E52',
            'Pfam_Domain':           '#8172B2',
            'KEGG_Pathway':          '#937860',
            'PROSITE_Annotation':    '#DA8BC3',
            None:                    '#CCCCCC',
        }
        EDGE_TYPE_COLORS = {
            'has_molecular_function': '#DD8452',
            'located_in_cellular_component': '#55A868',
            'involved_in_biological_process': '#C44E52',
            'has_pfam_domain': '#8172B2',
            'in_kegg_pathway': '#937860',
            'has_prosite_annotation': '#DA8BC3',
            'interacts_with': '#7B4173',
            'embeddings_similar_to': "#2D00E1",
            None: '#CCCCCC',
        }
        # Find the indices of vertices with a degree of 0
        isolated_vertices_indices = [v.index for v in g.vs if v.degree() == 0]
        # Delete the identified vertices
        g.delete_vertices(isolated_vertices_indices)
        
        vertex_colors = [NODE_TYPE_COLORS.get(v['node_type'], '#CCCCCC') for v in g.vs]
        degrees = g.degree()
        max_deg = max(degrees) if degrees else 1
        vertex_sizes = [20 + 60 * (d / max_deg) for d in degrees]
        labels = [str(v['name'])[:12] for v in g.vs]
        # edge_labels = [f"{e['edge_type']} ({e['weight']:.2f})" for e in g.es]
        edge_labels = [f"{e['weight']:.2f}" for e in g.es]
        node_connectivity = [g.degree(v.index) for v in g.vs]

        layout = g.layout("fr")
        ig.plot(
            g,
            output_path,
            layout=layout,
            bbox=(4800, 4800),
            margin=120,
            vertex_color=vertex_colors,
            vertex_size=vertex_sizes,
            vertex_label=[f"{label}\n({connectivity})" for label, connectivity in zip(labels, node_connectivity)],
            vertex_label_size=20,
            vertex_label_color='#111111',
            edge_width=1.6,
            edge_arrow_size=1.2,
            edge_color=[EDGE_TYPE_COLORS.get(e['edge_type'], '#CCCCCC') for e in g.es],
                edge_label=edge_labels,
                edge_label_size=16,
                edge_label_color='#333333'
        )

        present_types = {v['node_type'] for v in g.vs}
        patches = [mpatches.Patch(color=NODE_TYPE_COLORS.get(t, '#CCCCCC'), label=t or 'Unknown')
                   for t in present_types]
        fig, ax = plt.subplots(figsize=(3, len(patches) * 0.4 + 0.5))
        ax.legend(handles=patches, loc='center', frameon=False)
        ax.axis('off')
        legend_path = output_path.rsplit('.', 1)[0] + '_legend.png'
        fig.savefig(legend_path, bbox_inches='tight', dpi=150)
        plt.close(fig)

    def sample_node_subgraph(self, node_name):
        node_id = self.identify_node(node_name)
        if node_id is None:
            print(f"Cannot visualize graph for '{node_name}' because the node does not exist.")
            return

        neighbors = self.graph.neighbors(node_id, mode="all")
        subgraph = self.graph.subgraph([node_id] + neighbors)
        return subgraph

    def output_triples(self, output_path):
        with open(output_path, 'w') as f:
            for edge in self.graph.es:
                source = self.graph.vs[edge.source]["name"]
                target = self.graph.vs[edge.target]["name"]
                f.write(f"{source}\t{edge['edge_type']}\t{target}\n")

    def filter_edges_by_weight(self, min_weight):
        '''Removes all edges with weight below the specified threshold.'''
        edges_to_remove = [e.index for e in self.graph.es if e['weight'] < min_weight]
        self.graph.delete_edges(edges_to_remove)
        print(f"Removed {len(edges_to_remove)} edges with weight below {min_weight}. Remaining edges: {self.graph.ecount()}")

    def filter_nodes_by_degree(self, min_degree):
        '''Removes all nodes with degree below the specified threshold.'''
        nodes_to_remove = [v.index for v in self.graph.vs if self.graph.degree(v.index) < min_degree]
        self.graph.delete_vertices(nodes_to_remove)
        print(f"Removed {len(nodes_to_remove)} nodes with degree below {min_degree}. Remaining nodes: {self.graph.vcount()}")

    def output_nodes(self, output_path):
        with open(output_path, 'w') as f:
            for vertex in self.graph.vs:
                name = vertex["name"]
                node_type = vertex["node_type"]
                f.write(f"{name}\t{node_type}\n")

    def identify_node(self, search_name):
        '''Returns the graph ID of the node with the given name, or None if not found.'''
        search_name = self.ALIAS_MAPPER.get_node_name(search_name) or search_name
        node_id = self.graph.vs.find(name=search_name).index if self.graph.vs.select(name=search_name) else None
        if node_id is None:
            print(f"Node '{search_name}' not found in graph.")
        return node_id
        
    def get_node_relations(self, search_name):
        '''Returns node relations for the node with the given name.'''
        node_id = self.identify_node(search_name)
        if node_id is not None:
            node_name = self.ALIAS_MAPPER.get_node_name(search_name) or search_name
            neighbors = self.graph.neighbors(node_id, mode="all")
            relations = []
            for neighbor in neighbors:
                edge = self.graph.es.select(_source=node_id, _target=neighbor) or self.graph.es.select(_source=neighbor, _target=node_id)
                if edge:
                    edge_type = edge[0]['edge_type']
                    neighbor_name = self.graph.vs[neighbor]['name']
                    weight = edge[0]['weight']
                    sentence = f"{node_name} {edge_type} {neighbor_name} ({weight})"
                    relations.append(sentence)
            aliases = self.ALIAS_MAPPER.get_aliases(node_name)
            if aliases:
                relations.append(f"Aliases for {node_name}: {', '.join(aliases)}")
            return relations
        else:
            print(f"No relations found for '{search_name}' because the node does not exist.")
            return []

    def save_graphml(self, file_path):
        self.graph.write_graphml(file_path)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Build a graph from TSV files using igraph.")
    argparser.add_argument("tsv_folder", type=str, help="Path to the folder containing TSV files.")
    argparser.add_argument("output_graph", type=str, help="Path to save the output graph JSON file.")
    argparser.add_argument("--selection_percentage", type=float, default=100, help="Percentage of rows to select from each TSV file (default: 100).")
    args = argparser.parse_args()

    graph_builder = IGraphBuilder()

    files_in_folder = glob.glob(os.path.join(args.tsv_folder, "*.tsv"))
    for file in files_in_folder:
        print(f"Processing file: {file}")
        graph_builder.add_to_graph_tsv(file)
    graph_builder.filter_edges_by_weight(min_weight=0.3)  # Example threshold, adjust as needed
    graph_builder.filter_nodes_by_degree(min_degree=1)  # Example threshold, adjust as needed
    graph_builder.evaluate_graph()
    print()
    node_relations = graph_builder.get_node_relations("GO:0016773")  # Example node name, adjust as needed
    for relation in node_relations:
        print(relation)
    print()
    graph_builder.visualize_graph(graph_builder.sample_representative_subgraph(k = 4), "graph_visualization.png")
    graph_builder.visualize_graph(graph_builder.sample_node_subgraph("GAV06484.1"), "GAV06484.1_subgraph.png")
    graph_builder.output_triples("graph_triples.tsv")
    graph_builder.output_nodes("graph_nodes.tsv")
    graph_builder.save_graphml("graph.graphml")