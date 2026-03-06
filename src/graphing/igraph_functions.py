# script used to make predictions using the KG built by the builder

import torch
import igraph as ig
import numpy as np
import argparse
import pykeen
from pykeen.triples import TriplesFactory


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Predict links in a knowledge graph using R-GCN")
    argparser.add_argument("graphml_file", type=str, help="Path to the input graphml file")
    argparser.add_argument("--query_nodes", type=str, default="Dsup", help="Comma-separated list of query node IDs (e.g., 'protein1,trait1')")
    argparser.add_argument("--top_k", type=int, default=10, help="Number of predicted top linking proteins to output")
    args = argparser.parse_args()

    # open graph from graphml file with igraph

    # extract triples (source, label, target), edge-weights, and node-types from the graph

    # use TriplesFactory.from_labeled_triples to create a TriplesFactory for pykeen
    # tf = TriplesFactory.from_labeled_triples(triples)

    # fit an R-GCN model to the data using pykeen

    # make predictions for the specified query nodes and output the top-k predicted linking proteins for the query set

