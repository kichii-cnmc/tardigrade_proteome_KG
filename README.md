# tardigrade_proteome_KG
Used to build knowledge graphs for the classification and identification of proteins within the tardigrade proteome.

my-kg-project/
├── data/                   # Local raw data samples (gitignored)
│   ├── raw/                # Holds raw fasta files downloaded
│   ├── merged/             # Holds protein IDs and sequences after removing duplicates and fragments
│   ├── organized/          # Holds final protein information for KGs.
├── notebooks/              # For R&D and data exploration
├── scripts/                # Utility scripts for database setup/bash tasks
├── src/                    # Main source code
│   ├── collector/          # MODULE 1: Data Collection & Ingestion
│   │   ├── 1_download_genomes.sh       # bash script to download genomes from uniprot and ncbi
|   |   ├── 2_cdhit_merge_fastas.py     # used to merge uniprot and ncbi sequences, removing duplicates and fragments
|   |   └── 3_pull_uniprot_info.py      # references uniprot IDs to pull protein information into TSV file.
│   ├── builder/            # MODULE 2: KG Construction
│   │   ├── 
├── tests/                  # Unit and integration tests
├── environment.yaml        # Python dependencies
└── README.md               # Documentation