#!/bin/bash

# This script submits the collector pipeline job to the cluster scheduler for SBATCH processing
#SBATCH --job-name=KG_data_collection_pipeline
#SBATCH --error=logs/collector_pipeline_%j.err
#SBATCH --output=logs/collector_pipeline_%j.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Initialize micromamba environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate tardigrade_proteome_KG_env

# Run the collector pipeline script
bash scripts/1_collector_pipeline.sh
