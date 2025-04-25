#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J spades_hybrid_run
#SBATCH -t 04:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# === Load modules ===
module load bioinfo-tools spades/4.0.0

# === CONFIGURATION VARIABLES ===
export SRCDIR_ILLUMINA=/home/edman/genomeAnalysis/data/trimmed_data/dna
export SRCDIR_NANOPORE=/home/edman/genomeAnalysis/data/raw_data/nanopore
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/02_assembly/hybrid
export FASTA_DEST=/home/edman/genomeAnalysis/data/assembly_data/hybrid

# Copying input files to local temp storage
cp $SRCDIR_ILLUMINA/*_paired.fq.gz $SNIC_TMP/
cp $SRCDIR_NANOPORE/*.fasta.gz $SNIC_TMP/

# Move to local tmp
cd $SNIC_TMP || exit 1

# Unzip files
gunzip *.fq.gz
gunzip *.fasta.gz

# Get filenames (assuming only one set of paired reads and one nanopore file)
ILLUMINA_R1=$(ls *_1_paired.fq)
ILLUMINA_R2=$(ls *_2_paired.fq)
NANOPORE=$(ls *.fasta)

# Run SPAdes with both Illumina and Nanopore
spades.py \
  -1 "$ILLUMINA_R1" \
  -2 "$ILLUMINA_R2" \
  --nanopore "$NANOPORE" \
  -o spades_hybrid_output \
  -t 8 \
  -m 32

# === Post-processing ===

# Compress final contigs file
gzip spades_hybrid_output/contigs.fasta

# Create destination directories if they don't exist
mkdir -p "$FASTA_DEST"
mkdir -p "$RESULT_DIR"

# Move compressed contigs to assembly_data
mv spades_hybrid_output/contigs.fasta.gz "$FASTA_DEST/"

# Move the rest to analyses
shopt -s extglob
mv spades_hybrid_output/!(*.fasta.gz) "$RESULT_DIR/"
shopt -u extglob

