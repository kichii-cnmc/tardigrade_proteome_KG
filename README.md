# tardigrade_proteome_KG
Used to build knowledge graphs for the classification and identification of proteins within the tardigrade proteome.

```
my-kg-project/
├── data/                   # Local raw data samples (gitignored)
│   ├── raw/                # Holds raw fasta files downloaded
│   ├── merged/             # Holds protein IDs and sequences after removing duplicates and fragments
│   ├── organized/          # Holds final protein information for KGs.
|   └── graph_csvs/         # Holds CSV versions of data directly used to build KGs.
├── notebooks/              # For R&D and data exploration
|   ├── merged_tsv_validation.ipynb     # used to validate the merging of the UniProt and NCBI genome FASTAs into a single CSV
|   └── protein_info_validation.ipynb   # used to validate the protein info collected via UniProt API.
├── scripts/                # Utility scripts for database setup/bash tasks
|   └── 1_collector_pipeline.sh         # used to collect and process the FASTA files / genome data.
├── src/                    # Main source code
│   ├── collector/          # MODULE 1: Data Collection & Ingestion
│   │   ├── 1_download_genomes.sh       # bash script to download genomes from uniprot and ncbi
|   |   ├── 2_cdhit_merge_fastas.py     # used to merge uniprot and ncbi sequences, removing duplicates and fragments
|   |   └── 3_pull_uniprot_info.py      # references uniprot IDs to pull protein information into TSV file.
│   ├── builder/            # MODULE 2: KG Construction
│   │   ├── 1_build_KG_tables.py        # used to build the CSV files that are directly used to build KGs.
├── tests/                  # Unit and integration tests
├── environment.yaml        # Python dependencies
└── README.md               # Documentation
```
