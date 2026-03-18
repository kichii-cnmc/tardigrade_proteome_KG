#!/bin/bash
# This script is used to run the entire graphing pipeline for building and analyzing the knowledge graph

# run the graph building script to build the graph from the data/3_organized directory and save as graphml file
python3 src/graphing/igraph_builder.py data/3_organized/

# run the pyg graph function script to perform prediction and ranking
python3 src/graphing/pyg_graph_functions.py graph.graphml --query_nodes A0A1D1V419