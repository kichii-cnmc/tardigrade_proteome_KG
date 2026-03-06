# script used to make predictions using the KG built by the builder

import torch
import igraph as ig
import numpy as np
import argparse
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.utils import from_networkx
from torch_geometric.loader import DataLoader


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Run R-GCN link prediction on a graphml file")
    argparser.add_argument("--graphml_file", type=str, required=True, help="Path to the input graphml file")
    argparser.add_argument("--query_nodes", type=str, required=True, help="List of query node names for link prediction separated by commas")
    argparser.add_argument("--top_k", type=int, default=10, help="Number of top predictions to output")
    args = argparser.parse_args()

    # open graph from graphml file with igraph
    g = ig.Graph.Read_GraphML(args.graphml_file)
    # extract edge index, edge type, and edge norm as PyTorch tensors
    edge_index = torch.tensor(g.get_edgelist(), dtype=torch.long).t().contiguous()
    edge_types = torch.tensor(g.es["edge_type"], dtype=torch.long)
    edge_norm = torch.tensor(g.es["weight"], dtype=torch.float)
    node_types = torch.tensor(g.vs["node_type"], dtype=torch.long)

    # define relational graph convolutional network (R-GCN) model


    # add a link prediction layer to the R-GCN model (DistMult or ComplEx)


    # train the model on the graph data


    # input query node set (protein/traits) and rank correlation to the query node set


    # output top-k predicted links with their probabilities and supporting evidence from the graph (e.g., neighboring nodes and edges)
    pass