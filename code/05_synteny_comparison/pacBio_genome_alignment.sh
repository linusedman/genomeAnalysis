#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J nucmer_mapping
#SBATCH -t 04:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools MUMmer/4.0.0rc1

# CONFIGURATION VARIABLES
export ASSEMBLY_FILE=/home/edman/genomeAnalysis/data/assembly_data/pacBio/run2/assembly.contigs.fasta.gz
export REFERENCE_FILE=/home/edman/genomeAnalysis/data/reference_data/GCA_001750885.1.fasta.gz
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/05_synteny_comparison/pacBio/run2/mapping
mkdir -p "$RESULT_DIR"

# Copy files to local tmp storage
cp "$ASSEMBLY_FILE" "$SNIC_TMP/"
cp "$REFERENCE_FILE" "$SNIC_TMP/"

# Change to tmp storage
cd "$SNIC_TMP" || exit 1

# Unzip assembly if compressed
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY_FILE")
gunzip -f "$ASSEMBLY_BASENAME"
UNCOMPRESSED_ASSEMBLY="${ASSEMBLY_BASENAME%.gz}"

# Also get reference basename
REFERENCE_BASENAME=$(basename "$REFERENCE_FILE")
gunzip -f "$REFERENCE_BASENAME"
UNCOMPRESSED_REFERENCE="${REFERENCE_BASENAME%.gz}"

# RUN NUCmer
nucmer --prefix=assembly_vs_ref "$UNCOMPRESSED_REFERENCE" "$UNCOMPRESSED_ASSEMBLY"

# Move output files (delta file etc.) back to result dir
mv assembly_vs_ref.* "$RESULT_DIR/"

