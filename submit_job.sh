#!/bin/bash

# This script submits the collector pipeline job to the cluster scheduler for SBATCH processing
#SBATCH --job-name=KG_data_collection_pipeline
#SBATCH --error=logs/collector_pipeline_%j.err
#SBATCH --output=logs/collector_pipeline_%j.out
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu  # Use GPU partition

mkdir -p logs

# Initialize micromamba environment
module load micromamba
module load cuda # Load CUDA module if required for GPU support
eval "$(micromamba shell hook --shell bash)"
micromamba activate tardigrade_proteome_KG_env

# Run the collector pipeline script
python3 src/2_embedder/1_generate_embeddings.py data/3_organized/ embeddings_test_hpc.npz --test_mode
