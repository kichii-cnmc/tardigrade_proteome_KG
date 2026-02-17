#!/bin/bash

# This script is used to run the entire collector pipeline for downloading and processing genome data
# includes a test-mode to run a reduced version of the pipeline for testing purposes.

# Parse command line arguments for testing flag
TEST_MODE=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -t)
            TEST_MODE=true
            echo "Test mode enabled: The embedding generation and comparison will run on a smaller subset of data for testing purposes."
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [-t]"
            echo "  -t: Enable test mode (stop early)"
            exit 1
            ;;
    esac
done

# clear embedding generation log file
LOG_FILE="logs/embedding_generation.log"
if [ -f "$LOG_FILE" ]; then
    rm "$LOG_FILE"
    echo "Cleared existing embedding generation log file."
fi

# run the embedding generation script
python3 src/2_embedder/1_generate_embeddings.py data/3_organized data/1_raw/embeddings/raw_embeddings.npz $(if [ "$TEST_MODE" = true ]; then echo "--test_mode"; fi)
echo "Embedding generation completed."

# log generation time / size / test mode info
echo "Embedding generation log:" > logs/embedding_generation.log
echo "Date: $(date)" >> logs/embedding_generation.log
echo "Test mode: $TEST_MODE" >> logs/embedding_generation.log
if [ -f "data/1_raw/embeddings/raw_embeddings.npz" ]; then
    EMBEDDING_SIZE=$(du -h data/1_raw/embeddings/raw_embeddings.npz | cut -f1)
    echo "Embedding file size: $EMBEDDING_SIZE" >> logs/embedding_generation.log
    echo "Embedding generation completed successfully." >> logs/embedding_generation.log
else
    echo "Embedding file not found. Embedding generation may have failed." >> logs/embedding_generation.log
fi 

# apply pca dimensionality reduction to the embeddings
python3 src/2_embedder/2_pca_reduction.py data/1_raw/embeddings/raw_embeddings.npz data/1_raw/embeddings/reduced_embeddings.npz --variance_threshold 0.90

# run the embedding comparison script to calculate cosine similarity matrix & save edge list for graph building
python3 src/2_embedder/3_compare_embeddings.py data/1_raw/embeddings/reduced_embeddings.npz data/3_organized $(if [ "$TEST_MODE" = true ]; then echo "--test_mode"; fi)
echo "Embedding comparison completed."