#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J mummer_showcoords
#SBATCH -t 01:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools MUMmer/4.0.0rc1

# === CONFIGURATION VARIABLES ===
export DELTA_FILE=/home/edman/genomeAnalysis/analyses/05_synteny_comparison/pacBio/run2/mapping/assembly_vs_ref.delta
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/05_synteny_comparison/pacBio/run2/mapping
mkdir -p "$RESULT_DIR"

# Copy delta file to local tmp storage
cp "$DELTA_FILE" "$SNIC_TMP/"
cd "$SNIC_TMP" || exit 1

DELTA_BASENAME=$(basename "$DELTA_FILE")

# === Run show-coords ===
show-coords -rcl "$DELTA_BASENAME" > comparison.coords

# === Move output file back ===
mv comparison.coords "$RESULT_DIR/"

