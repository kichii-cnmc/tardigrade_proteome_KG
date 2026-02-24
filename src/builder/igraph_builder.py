# script that uses igraph to build and manage the graph of the project
import igraph as ig
import os
import json
import pandas as pd
import glob
import argparse

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
        'UniProt_ID': 'Protein',
        'PPI_target': 'Protein'
    }
    GRAPH_NODE_ALIAS_TYPES = ['NCBI_ID', 'Organism', 'Gene_Name', 'AF'] # columns that can be added as node attributes

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
        if second_col in self.GRAPH_EDGE_LABELS:
            edge_label, weight_col = self.GRAPH_EDGE_LABELS[second_col]
            edges_list = []
            for _, row in df.iterrows():
                source = row.iloc[0]
                target = row.iloc[1]
                weight = row[weight_col] if weight_col and weight_col in df.columns else 1
                edges_list.append((source, target, weight))
            print(f"Adding edges of type '{edge_label}' from {tsv_file_path}...")
            # add nodes first to ensure all vertices exist before adding edges, then add edges in bulk
            node_type = self.GRAPH_NODE_LABELS.get(second_col, "UnknownType")
            if node_type == "UnknownType":
                print(f"    Warning: No node type mapping found for column '{second_col}', defaulting to 'UnknownType'.")
            self.add_nodes_from_list([target for _, target, _ in edges_list], node_type=node_type)
            self.add_edges_from_list(edges_list, edge_type=edge_label)
        elif second_col in self.GRAPH_NODE_LABELS:
            node_type = self.GRAPH_NODE_LABELS[second_col]
            print(f"Adding nodes of type '{node_type}' from {tsv_file_path}...")
            self.add_nodes_from_list(df.iloc[:, 1].tolist(), node_type=node_type)
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
            print(f"({source}, {edge['edge_type']}, {target})")
        print("Example Nodes:")
        for v in self.graph.vs[:5]:
            print(f"({v['name']}, {v['node_type']})")

    def get_graph(self):
        return self.graph
    
    def _sample_representative_subgraph(self, k=6, max_edges=600):
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

    def visualize_graph(self, output_path, sample_k=15):
        import matplotlib
        matplotlib.use('Agg')  # non-interactive backend — must be set before pyplot import
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

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

        g = self._sample_representative_subgraph(k=sample_k)

        vertex_colors = [NODE_TYPE_COLORS.get(v['node_type'], '#CCCCCC') for v in g.vs]
        degrees = g.degree()
        max_deg = max(degrees) if degrees else 1
        vertex_sizes = [10 + 30 * (d / max_deg) for d in degrees]
        labels = [str(v['name'])[:12] for v in g.vs]

        layout = g.layout("fr")
        ig.plot(
            g,
            output_path,
            layout=layout,
            bbox=(2400, 2400),
            margin=80,
            vertex_color=vertex_colors,
            vertex_size=vertex_sizes,
            vertex_label=labels,
            vertex_label_size=7,
            vertex_label_color='#111111',
            edge_width=0.6,
            edge_arrow_size=0.4,
            edge_color='#AAAAAA',
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

    def output_triples(self, output_path):
        with open(output_path, 'w') as f:
            for edge in self.graph.es:
                source = self.graph.vs[edge.source]["name"]
                target = self.graph.vs[edge.target]["name"]
                f.write(f"{source}\t{edge['edge_type']}\t{target}\n")

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
    graph_builder.evaluate_graph()
    graph_builder.visualize_graph("graph_visualization.png")
    graph_builder.output_triples("graph_triples.tsv")
