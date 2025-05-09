#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J quast_hybrid_eval
#SBATCH -t 02:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools quast/5.0.2

# === CONFIGURATION VARIABLES ===
export ASSEMBLY_FILE=/home/edman/genomeAnalysis/data/assembly_data/hybrid/trim/contigs.fasta.gz
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/03_assembly_eval/hybrid/trim
mkdir -p "$RESULT_DIR"

# Copying the input file to local temporary storage
cp "$ASSEMBLY_FILE" "$SNIC_TMP/"

# Changing to local temporary directory
cd "$SNIC_TMP" || exit 1

# Getting the basename of the assembly file
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY_FILE")

# Unzipping the file
gunzip -f "$ASSEMBLY_BASENAME"

# Updating the variable to refer to the decompressed file
UNCOMPRESSED_ASSEMBLY="${ASSEMBLY_BASENAME%.gz}"

# Running QUAST on the uncompressed assembly    
quast.py "$UNCOMPRESSED_ASSEMBLY" -o "$RESULT_DIR" -t 4


