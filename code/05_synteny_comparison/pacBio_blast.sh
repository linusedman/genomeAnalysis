#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J blastn_comparison
#SBATCH -t 01:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools blast/2.15.0+

# CONFIGURATION VARIABLES
export QUERY_GZ=/home/edman/genomeAnalysis/data/assembly_data/pacBio/assembly.contigs.fasta.gz
export SUBJECT_GZ=/home/edman/genomeAnalysis/data/reference_data/GCA_001750885.1.fasta.gz
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/05_synteny_comparison/blast_output
mkdir -p "$RESULT_DIR"

# Copy compressed files to temp directory
cp "$QUERY_GZ" "$SNIC_TMP/"
cp "$SUBJECT_GZ" "$SNIC_TMP/"
cd "$SNIC_TMP" || exit 1

# Uncompress FASTA files
gunzip assembly.contigs.fasta.gz
gunzip GCA_001750885.1.fasta.gz

# Create BLAST database for reference genome
makeblastdb -in GCA_001750885.1.fasta -dbtype nucl -out ref_db

# Run BLAST using contigs as query
blastn -query assembly.contigs.fasta -db ref_db -outfmt 6 -out comparison.tab

# Move results back
mv comparison.tab "$RESULT_DIR/"
