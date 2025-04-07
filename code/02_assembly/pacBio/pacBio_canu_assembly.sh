#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J canu_pacbio_run
#SBATCH -t 06:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Loading  modules
module load bioinfo-tools
module load canu/2.2

# === CONFIGURATION VARIABLES ===
export SRCDIR=/home/edman/genomeAnalysis/data/raw_data/pacBio
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/02_assembly/pacBio
export FASTA_DEST=/home/edman/genomeAnalysis/data/assembly_data

# Copying input files (fastq.gz) to local temp storage
cp $SRCDIR/*.fastq.gz $SNIC_TMP/

# Moving to the local temp directory (exit if failed)
cd $SNIC_TMP || exit 1

# Unziping fastq.gz files
gunzip *.fastq.gz

# Running Canu
canu -p assembly -d canu_run genomeSize=3.5m -pacbio *.fastq useGrid=false maxThreads=4

# === Post-processing ===

# Compressing any uncompressed FASTA files from the Canu run
for f in canu_run/*.fasta; do
    [ -f "$f" ] && gzip "$f"
done

# Creating destination directories if they do not exist (the should but just in case)
mkdir -p "$FASTA_DEST"
mkdir -p "$RESULT_DIR"

# Moving compressed FASTA files to the assembly_data directory
mv canu_run/*.fasta.gz "$FASTA_DEST/"

# Moving all other files to the analyses directory
shopt -s extglob
mv canu_run/!(*.fasta.gz) "$RESULT_DIR/"
shopt -u extglob

