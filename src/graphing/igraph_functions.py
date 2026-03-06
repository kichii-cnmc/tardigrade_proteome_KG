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

    # open graph from graphml file with igraph

    # extract edge index, edge type, and edge norm as PyTorch tensors

    # define relational graph convolutional network (R-GCN) model


    # add a link prediction layer to the R-GCN model (DistMult or ComplEx)


    # train the model on the graph data


    # input query node set (protein/traits) and rank link probabilities to the query node set


    # output top-k predicted links with their probabilities and supporting evidence from the graph (e.g., neighboring nodes and edges)

