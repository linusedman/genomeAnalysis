#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J prokka_pacbio
#SBATCH -t 04:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools prokka/1.45-5b58020

# CONFIGURATION VARIABLES
export ASSEMBLY_FILE=/home/edman/genomeAnalysis/data/assembly_data/pacBio/assembly.contigs.fasta.gz
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/04_annotation/pacBio
mkdir -p "$RESULT_DIR"

# Copy the compressed assembly file to local temporary storage
cp "$ASSEMBLY_FILE" "$SNIC_TMP/"

# Change to the temporary directory
cd "$SNIC_TMP" || exit 1

# Get the basename of the assembly file (e.g., assembly.contigs.fasta.gz)
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY_FILE")

# Unzip the assembly file (force option to overwrite if necessary)
gunzip -f "$ASSEMBLY_BASENAME"

# Update variable to refer to the uncompressed file (strip off .gz)
UNCOMPRESSED_ASSEMBLY="${ASSEMBLY_BASENAME%.gz}"

# Running Prokka on the PacBio assembly
prokka --force --outdir "$RESULT_DIR" --prefix pacbio --cpus 4 "$UNCOMPRESSED_ASSEMBLY"

