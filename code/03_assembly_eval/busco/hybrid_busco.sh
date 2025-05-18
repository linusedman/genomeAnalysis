#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 20:00:00
#SBATCH -J busco_hybrid_eval
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools BUSCO/5.7.1

# Set up Augustus config
source $AUGUSTUS_CONFIG_COPY

# Define paths
ASSEMBLY_FILE=/home/edman/genomeAnalysis/data/assembly_data/hybrid/contigs.fasta.gz
RESULT_DIR=/home/edman/genomeAnalysis/analyses/03_assembly_eval/hybrid/busco
LINEAGE=$BUSCO_LINEAGE_SETS/lactobacillales_odb10

# Prepare output and temp dirs
mkdir -p "$RESULT_DIR"
cp "$ASSEMBLY_FILE" "$SNIC_TMP/"
cd "$SNIC_TMP" || exit 1

# Unzip assembly
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY_FILE")
gunzip -f "$ASSEMBLY_BASENAME"
UNCOMPRESSED_ASSEMBLY="${ASSEMBLY_BASENAME%.gz}"

# Run BUSCO with offline mode enabled
busco -i "$UNCOMPRESSED_ASSEMBLY" \
      -o busco_output \
      -m genome \
      -l "$LINEAGE" \
      -f \
      --offline \
      --cpu 16 \
      --out_path "$RESULT_DIR"

