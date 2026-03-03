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
    
    def output_report(self, output_file):
        '''Outputs the validation report as a TSV file.'''
        report_df = pd.DataFrame(list(self.report.items()), columns=['Check', 'Result'])
        print(f"\nValidation Report for {self.dataset_name}:\n")
        print(report_df)
        report_df.to_csv(output_file, sep='\t', index=False)

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
        self.report['Base ID List Size'] = len(self.base_id_list)
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
        self.report['Percentage of IDs Not in Base List'] = str((len(inconsistent_ids) / len(dataset_ids) * 100) if len(dataset_ids) > 0 else 0) + "%"
        # ids present in base list but not in dataset (coverage gap for dataset)
        missing_ids = set(self.base_id_list) - dataset_ids
        self.report['Number of IDs in Base List Not in Dataset'] = len(missing_ids)
        self.report['Percentage of IDs in Base List Not in Dataset'] = str((len(missing_ids) / len(self.base_id_list) * 100) if len(self.base_id_list) > 0 else 0) + "%"
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
            self.report['Score Mean'] = round(self.data['score'].mean(), 4)
            self.report['Score Median'] = round(self.data['score'].median(), 4)
            self.report['Score Range'] = round(self.data['score'].min(), 4), round(self.data['score'].max(), 4)
            if plot_output:
                plt.figure(figsize=(10, 6))
                sns.histplot(self.data['score'], bins=50, kde=True)
                plt.title(f'Score Distribution for {self.dataset_name}')
                plt.xlabel('Score')
                plt.ylabel('Frequency')
                plt.savefig(plot_output)
            return self.report
        elif 'embedding_similarity' in self.data.columns:
            self.report['Similarity Mean'] = round(self.data['embedding_similarity'].mean(), 4)
            self.report['Similarity Median'] = round(self.data['embedding_similarity'].median(), 4)
            self.report['Similarity Range'] = round(self.data['embedding_similarity'].min(), 4), round(self.data['embedding_similarity'].max(), 4)
            if plot_output:
                plt.figure(figsize=(10, 6))
                sns.histplot(self.data['embedding_similarity'], bins=50, kde=True)
                plt.title(f'Embedding Similarity Distribution for {self.dataset_name}')
                plt.xlabel('Embedding Similarity')
                plt.ylabel('Frequency')
                plt.savefig(plot_output)
            return self.report
        else:
            print(f"Warning: No score or embedding_similarity column found in {self.dataset_name} for score distribution measurement.")
            return None
        
    def measure_relationships_per_source(self, plot_output=None):
        '''Measures the number of relationships/interactions per source protein and adds it to the report, plotting a histogram if specified.'''
        if 'PPI_source' in self.data.columns:
            relationships_per_source = self.data.groupby('PPI_source').size()
            self.report['Relationships per Source Mean'] = round(relationships_per_source.mean(), 2)
            self.report['Relationships per Source Median'] = round(relationships_per_source.median(), 2)
            self.report['Relationships per Source Range'] = round(relationships_per_source.min(), 2), round(relationships_per_source.max(), 2)
            if plot_output:
                plt.figure(figsize=(10, 6))
                sns.histplot(relationships_per_source, bins=50, kde=True)
                plt.title(f'Relationships per Source Protein for {self.dataset_name}')
                plt.xlabel('Number of Relationships')
                plt.ylabel('Frequency')
                plt.savefig(plot_output)
            return self.report
        elif 'protein1' in self.data.columns:
            relationships_per_source = self.data.groupby('protein1').size()
            self.report['Relationships per Source Mean'] = round(relationships_per_source.mean(), 2)
            self.report['Relationships per Source Median'] = round(relationships_per_source.median(), 2)
            self.report['Relationships per Source Range'] = round(relationships_per_source.min(), 2), round(relationships_per_source.max(), 2)
            if plot_output:
                plt.figure(figsize=(10, 6))
                sns.histplot(relationships_per_source, bins=50, kde=True)
                plt.title(f'Relationships per Source Protein for {self.dataset_name}')
                plt.xlabel('Number of Relationships')
                plt.ylabel('Frequency')
                plt.savefig(plot_output)
            return self.report
        else:
            print(f"Warning: No PPI_source or protein1 column found in {self.dataset_name} for relationships per source measurement.")
            return None

if __name__ == "__main__":
    # Load base UniProt ID list from the protein info datasets for consistency checks
    rv_uniprot_id_list = load_protein_info('data/3_organized/RV_UniProt_ID.tsv')['UniProt_ID'].tolist()
    he_uniprot_id_list = load_protein_info('data/3_organized/HE_UniProt_ID.tsv')['UniProt_ID'].tolist()
    merged_uniprot_id_list = rv_uniprot_id_list + he_uniprot_id_list

    # Create DataValidator instances for each dataset for each organism
    # sequence (2)
    rv_sequence_validator = DataValidator('data/3_organized/RV_Sequence.tsv', rv_uniprot_id_list, 'RV Sequence')
    he_sequence_validator = DataValidator('data/3_organized/HE_Sequence.tsv', he_uniprot_id_list, 'HE Sequence')
    rv_sequence_validator.run_common_checks()
    he_sequence_validator.run_common_checks()
    rv_sequence_validator.output_report('data/4_analysis/RV_Sequence_validation_report.tsv')
    he_sequence_validator.output_report('data/4_analysis/HE_Sequence_validation_report.tsv')
    # geneid (2)
    rv_geneid_validator = DataValidator('data/3_organized/RV_GeneID.tsv', rv_uniprot_id_list, 'RV GeneID')
    he_geneid_validator = DataValidator('data/3_organized/HE_GeneID.tsv', he_uniprot_id_list, 'HE GeneID')
    rv_geneid_validator.run_common_checks()
    he_geneid_validator.run_common_checks()
    rv_geneid_validator.output_report('data/4_analysis/RV_GeneID_validation_report.tsv')
    he_geneid_validator.output_report('data/4_analysis/HE_GeneID_validation_report.tsv')
    # genename (2)
    # rv_genename_validator = DataValidator('data/3_organized/RV_GeneName.tsv', rv_uniprot_id_list, 'RV GeneName')
    # he_genename_validator = DataValidator('data/3_organized/HE_GeneName.tsv', he_uniprot_id_list, 'HE GeneName')
    # rv_genename_validator.run_common_checks()
    # he_genename_validator.run_common_checks()
    # rv_genename_validator.output_report('data/4_analysis/RV_GeneName_validation_report.tsv')
    # he_genename_validator.output_report('data/4_analysis/HE_GeneName_validation_report.tsv')
    # af structures (2)
    rv_af_validator = DataValidator('data/3_organized/RV_AF_structures.tsv', rv_uniprot_id_list, 'RV AF Structures')
    he_af_validator = DataValidator('data/3_organized/HE_AF_structures.tsv', he_uniprot_id_list, 'HE AF Structures')
    rv_af_validator.run_common_checks()
    he_af_validator.run_common_checks()
    rv_af_validator.output_report('data/4_analysis/RV_AF_structures_validation_report.tsv')
    he_af_validator.output_report('data/4_analysis/HE_AF_structures_validation_report.tsv')
    # deepgo (BP, MF, CC) (1)
    # pfam (2)
    # prosite (2)
    # string ppi (2)
    # esm-2 similarity (1)


