# script used to generate a validation report and figures for data in 3_organized

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

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

def load_protein_info(protein_info_file):
    '''Loads the protein info dataset and returns a DataFrame.'''
    df = pd.read_csv(protein_info_file, sep='\t')
    return df


class DataValidator:
    def __init__(self, filepath, base_id_list, dataset_name):
        self.filepath = filepath
        self.dataset_name = dataset_name
        self.data = load_protein_info(filepath)
        self.base_id_list = base_id_list # list of valid UniProt IDs to check against for consistency
        self.report = {}

    def return_name(self):
        return self.dataset_name

    def count_protein_sources(self):
        '''Checks the number of unique and non-unique protein sources in the dataset and adds it to the report.'''
        unique_proteins = set()
        non_unique_proteins_count = 0
        if self.data.columns[0] == 'UniProt_ID':
            unique_proteins = set(self.data['UniProt_ID'].dropna().astype(str).tolist())
            non_unique_proteins_count = len(self.data) - len(unique_proteins)
        elif self.data.columns[0] == 'PPI_source':
            unique_proteins = set(self.data['PPI_source'].dropna().astype(str).tolist() + self.data['PPI_target'].dropna().astype(str).tolist())
        self.report['Number of Unique Proteins'] = len(unique_proteins)
        self.report['Number of Non-Unique or Duplicate Proteins'] = non_unique_proteins_count
        return len(unique_proteins)
    
    def check_id_consistency(self):
        '''Checks that all UniProt IDs in the dataset are present in the base ID list.'''
        if 'UniProt_ID' in self.data.columns:
            dataset_ids = set(self.data['UniProt_ID'].dropna().astype(str).tolist())
        elif 'PPI_source' in self.data.columns and 'PPI_target' in self.data.columns:
            dataset_ids = set(self.data['PPI_source'].dropna().astype(str).tolist() + self.data['PPI_target'].dropna().astype(str).tolist())
        else:
            print(f"Warning: No UniProt_ID or PPI_source/PPI_target columns found in {self.dataset_name} for consistency check.")
            return None
        # ids present in dataset but not in base list
        inconsistent_ids = dataset_ids - set(self.base_id_list)
        self.report['Number of IDs Not in Base List'] = len(inconsistent_ids)
        self.report['Percentage of IDs Not in Base List'] = (len(inconsistent_ids) / len(dataset_ids) * 100) if len(dataset_ids) > 0 else 0
        # ids present in base list but not in dataset (coverage gap for dataset)
        missing_ids = set(self.base_id_list) - dataset_ids
        self.report['Number of IDs in Base List Not in Dataset'] = len(missing_ids)
        self.report['Percentage of IDs in Base List Not in Dataset'] = (len(missing_ids) / len(self.base_id_list) * 100) if len(self.base_id_list) > 0 else 0
        return inconsistent_ids, missing_ids
    
    def run_common_checks(self):
        '''Runs all common or universal checks and returns the report.'''
        self.count_protein_sources()
        self.check_id_consistency()
        return self.report
    
    def count_unique_relationships(self):
        '''Counts the number of unique relationships/interactions in the dataset (single direction) and adds it to the report.'''
        if 'PPI_source' in self.data.columns and 'PPI_target' in self.data.columns:
            unique_relationships = list(zip(self.data['PPI_source'].dropna().astype(str).tolist(), self.data['PPI_target'].dropna().astype(str).tolist()))
            set_unique_relationships = set(unique_relationships)
            if len(unique_relationships) != len(set_unique_relationships):
                print(f"Warning: Found {len(unique_relationships) - len(set_unique_relationships)} duplicate relationships in {self.dataset_name}.")
                self.report['!! Number of Duplicate Relationships'] = len(unique_relationships) - len(set_unique_relationships)
            self.report['Number of Unique Relationships'] = len(set_unique_relationships)
            return len(set_unique_relationships)
        elif 'protein1' in self.data.columns and 'protein2' in self.data.columns:
            unique_relationships = list(zip(self.data['protein1'].dropna().astype(str).tolist(), self.data['protein2'].dropna().astype(str).tolist()))
            set_unique_relationships = set(unique_relationships)
            if len(unique_relationships) != len(set_unique_relationships):
                print(f"Warning: Found {len(unique_relationships) - len(set_unique_relationships)} duplicate relationships in {self.dataset_name}.")
                self.report['!! Number of Duplicate Relationships'] = len(unique_relationships) - len(set_unique_relationships)
            self.report['Number of Unique Relationships'] = len(set_unique_relationships)
            return len(set_unique_relationships)
        else:
            print(f"Warning: No PPI_source/PPI_target or protein1/protein2 columns found in {self.dataset_name} for relationship counting.")
            return None
        
    def measure_relationship_score_distribution(self, plot_output=None):
        '''Measures the distribution of scores in the dataset and adds it to the report.'''
        if 'score' in self.data.columns:
            self.report['Score Mean'] = self.data['score'].mean()
            self.report['Score Median'] = self.data['score'].median()
            self.report['Score Range'] = self.data['score'].min(), self.data['score'].max()
            if plot_output:
                plt.figure(figsize=(10, 6))
                sns.histplot(self.data['score'], bins=50, kde=True)
                plt.title(f'Score Distribution for {self.dataset_name.replace("_", " ").replace(".tsv", "")}')
                plt.xlabel('Score')
                plt.ylabel('Frequency')
                plt.savefig(plot_output)
            return self.report
        elif 'embedding_similarity' in self.data.columns:
            self.report['Similarity Mean'] = self.data['embedding_similarity'].mean()
            self.report['Similarity Median'] = self.data['embedding_similarity'].median()
            self.report['Similarity Range'] = self.data['embedding_similarity'].min(), self.data['embedding_similarity'].max()
            if plot_output:
                plt.figure(figsize=(10, 6))
                sns.histplot(self.data['embedding_similarity'], bins=50, kde=True)
                plt.title(f'Embedding Similarity Distribution for {self.dataset_name.replace("_", " ").replace(".tsv", "")}')
                plt.xlabel('Embedding Similarity')
                plt.ylabel('Frequency')
                plt.savefig(plot_output)
            return self.report
        else:
            print(f"Warning: No score or embedding_similarity column found in {self.dataset_name} for score distribution measurement.")
            return None
        


if __name__ == "__main__":
    # Load base UniProt ID list from the protein info datasets for consistency checks
    rv_uniprot_id_list = load_protein_info('data/3_organized/RV_UniProt_ID.tsv')['UniProt_ID'].tolist()
    he_uniprot_id_list = load_protein_info('data/3_organized/HE_UniProt_ID.tsv')['UniProt_ID'].tolist()
    merged_uniprot_id_list = rv_uniprot_id_list + he_uniprot_id_list

    # Create DataValidator instances for each dataset for each organism
    # sequence (2)
    # geneid (2)
    # genename (2)
    # af structures (1)
    # deepgo (BP, MF, CC) (1)
    # pfam (2)
    # prosite (2)
    # string ppi (2)
    # esm-2 similarity (1)


