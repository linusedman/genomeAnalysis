#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J quast_pacbio_eval
#SBATCH -t 02:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools quast/5.0.2

# CONFIGURATION VARIABLES
export ASSEMBLY_FILE=/home/edman/genomeAnalysis/data/assembly_data/pacBio/assembly.contigs.fasta.gz
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/03_assembly_eval/pacBio/quast
mkdir -p "$RESULT_DIR"

# Define the directory containing PacBio reads
export PACBIO_READS_DIR=/home/edman/genomeAnalysis/data/raw_data/pacBio

# Copy the assembly file to local temporary storage
cp "$ASSEMBLY_FILE" "$SNIC_TMP/"

# Copy all PacBio read files to local temporary storage
# (Assuming the read files have a *.subreads.fastq.gz naming pattern)
cp "$PACBIO_READS_DIR"/*.subreads.fastq.gz "$SNIC_TMP/"

# Change to the temporary directory
cd "$SNIC_TMP" || exit 1

# Get the basename of the assembly file (e.g., assembly.contigs.fasta.gz)
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY_FILE")

# Unzip the assembly file
gunzip -f "$ASSEMBLY_BASENAME"

# Update variable to refer to the uncompressed assembly file
UNCOMPRESSED_ASSEMBLY="${ASSEMBLY_BASENAME%.gz}"

# Merge the PacBio reads into a single file
# The resulting file will be 'pacbio_reads.fastq.gz'
cat *.subreads.fastq.gz > pacbio_reads.fastq.gz

# Running QUAST on the uncompressed assembly using the merged PacBio reads
quast.py "$UNCOMPRESSED_ASSEMBLY" --pacbio pacbio_reads.fastq.gz -o "$RESULT_DIR" -t 8

