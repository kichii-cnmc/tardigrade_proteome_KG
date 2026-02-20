#!/bin/bash

# This script submits the collector pipeline job to the cluster scheduler for SBATCH processing
#SBATCH --job-name=KG_data_collection_pipeline
#SBATCH --error=logs/collector_pipeline_%j.err
#SBATCH --output=logs/collector_pipeline_%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

mkdir -p logs

# Initialize micromamba environment
module load micromamba
module load cuda # Load CUDA module if required for GPU support
eval "$(micromamba shell hook --shell bash)"
micromamba activate tardigrade_proteome_KG_env

# run the embedding pipeline script
bash scripts/2_embedding_pipeline.sh 
