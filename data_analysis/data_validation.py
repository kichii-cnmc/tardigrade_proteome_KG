# script used to generate a validation report and figures for data in 3_organized

import pandas as pd

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

if __name__ == "__main__":
    # Load base UniProt ID list from the protein info datasets for consistency checks
    rv_uniprot_id_list = load_protein_info('data/3_organized/RV_UniProt_ID.tsv')['UniProt_ID'].tolist()
    he_uniprot_id_list = load_protein_info('data/3_organized/HE_UniProt_ID.tsv')['UniProt_ID'].tolist()
    merged_uniprot_id_list = rv_uniprot_id_list + he_uniprot_id_list

    # Create DataValidator instances for each dataset for each organism
    # sequence (2)
    # geneid (2)
    # af structures (1)
    # deepgo (BP, MF, CC) (1)
    # pfam (2)
    # prosite (2)
    # string ppi (2)
    # esm-2 similarity (1)


