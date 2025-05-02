#!/bin/bash

# === CONFIGURATION ===
SOURCE_DIR_BH="/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH"
SOURCE_DIR_SERUM="/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum"

TARGET_DIR_BH="/home/edman/genomeAnalysis/data/raw_data/illumina/RNA/seq_BH"
TARGET_DIR_SERUM="/home/edman/genomeAnalysis/data/raw_data/illumina/RNA/seq_serum"

# Create target directories if they don't exist
mkdir -p "$TARGET_DIR_BH"
mkdir -p "$TARGET_DIR_SERUM"

# === Create symlinks for BH ===
for file in "$SOURCE_DIR_BH"/trim_paired*; do
    ln -s "$file" "$TARGET_DIR_BH/"
done

# === Create symlinks for Serum ===
for file in "$SOURCE_DIR_SERUM"/trim_paired*; do
    ln -s "$file" "$TARGET_DIR_SERUM/"
done

echo "Symlinks created successfully!"

