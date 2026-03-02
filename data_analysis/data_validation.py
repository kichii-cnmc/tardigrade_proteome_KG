# script used to generate a validation report and figures for data in 3_organized

import pandas as pd
import argparse

'''
Script Checks:
- Number of unique proteins in each dataset (for every dataset)
- UniProt ID and Sequence Coverage
- Consistency of IDs across datasets (for every dataset, check that all UniProt IDs are present in the protein info dataset)

Relationship-type Dataset Checks (STRING PPI, ESM-2 similarity):
- Number of unique interactions
- Score distribution (for STRING PPI and ESM-2 similarity)
- Clustering of interactions; are there distinct clusters of associations/interactions?

Assignment-type Dataset Checks (DeepGO, Pfam, PROSITE):
- Number of unique protein-function/domain assignments
- Score distribution (for DeepGO)
- Histogram of number of functions/domains assigned per protein; are there proteins with many assignments vs. few?
    - Mean, median, range of number of assignments per protein
- Histogram of number of proteins assigned per function/domain; are there functions/domains that are very common vs. rare?
    - Mean, median, range of number of proteins per function/domain
'''

