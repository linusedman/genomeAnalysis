#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 20:00:00
#SBATCH -J busco_pacbio_eval
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools BUSCO/5.7.1

# Setup AUGUSTUS config in SNIC_TMP
export AUGUSTUS_CONFIG_PATH="$SNIC_TMP/augustus_config"
mkdir -p "$AUGUSTUS_CONFIG_PATH"
cp -r "$AUGUSTUS_CONFIG_COPY_DIR"/* "$AUGUSTUS_CONFIG_PATH"

# Configuration
export ASSEMBLY_FILE=/home/edman/genomeAnalysis/data/assembly_data/pacBio/run2/assembly.contigs.fasta.gz
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/03_assembly_eval/pacBio/run2/busco
export LINEAGE=$BUSCO_LINEAGE_SETS/lactobacillales_odb10

mkdir -p "$RESULT_DIR"

# Copy input to SNIC_TMP
cp "$ASSEMBLY_FILE" "$SNIC_TMP/"
cd "$SNIC_TMP" || exit 1

# Unzip assembly
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY_FILE")
gunzip -f "$ASSEMBLY_BASENAME"
UNCOMPRESSED_ASSEMBLY="${ASSEMBLY_BASENAME%.gz}"

# Run BUSCO
busco -i "$UNCOMPRESSED_ASSEMBLY" \
      -o busco_output \
      -m genome \
      -l "$LINEAGE" \
      --cpu 8 \
      --out_path "$RESULT_DIR"

