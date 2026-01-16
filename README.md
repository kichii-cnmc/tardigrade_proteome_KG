# tardigrade_proteome_KG
Used to build knowledge graphs for the classification and identification of proteins within the tardigrade proteome.

my-kg-project/
├── data/                   # Local raw data samples (gitignored)
├── notebooks/              # For R&D and data exploration
├── scripts/                # Utility scripts for database setup/bash tasks
├── src/                    # Main source code
│   ├── collector/          # MODULE 1: Data Collection & Ingestion
│   │   └── 
│   ├── builder/            # MODULE 2: KG Construction
│   │   ├── schema/         # Ontology/Schema definitions (JSON-LD, OWL)
│   │   ├── mapping/        # Logic to map raw data to graph entities
│   │   ├── loaders/        # Logic to push data to Neo4j, ArangoDB, etc.
│   │   └── build_kg.py     # Main entry point for KG creation
├── tests/                  # Unit and integration tests
├── config.yaml             # lists which genomes and database information to fetch
├── environment.yaml        # Python dependencies
└── README.md               # Documentation