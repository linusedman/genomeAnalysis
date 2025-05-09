#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 07:00:00
#SBATCH -J htseq_serum
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools htseq/2.0.2 samtools

# === CONFIGURATION ===
BAM_DIR=/home/edman/genomeAnalysis/analyses/06_post_mapping/serum
GFF_ORIG=/home/edman/genomeAnalysis/analyses/04_annotation/pacBio/run2/pacbio.gff
OUT_DIR=/home/edman/genomeAnalysis/analyses/07_diff_exp_analysis/read_counting/serum
mkdir -p "$OUT_DIR"

# === Copy files to temp storage ===
cp "$GFF_ORIG" "$SNIC_TMP/"
cp "$BAM_DIR"/ERR*_sorted.bam "$SNIC_TMP/"
cp "$BAM_DIR"/ERR*_sorted.bam.bai "$SNIC_TMP/"

cd "$SNIC_TMP"

# === Strip FASTA section from GFF ===
GFF_CLEAN=pacbio_no_fasta.gff
awk '/^##FASTA/ {exit} {print}' "$(basename "$GFF_ORIG")" > "$GFF_CLEAN"

# === Sample list ===
SAMPLES=(ERR1797969 ERR1797970 ERR1797971)

# === Run htseq-count ===
for SAMPLE in "${SAMPLES[@]}"; do
    htseq-count \
      -f bam \
      -r pos \
      -s no \
      -t CDS \
      -i ID \
      "${SAMPLE}_sorted.bam" \
      "$GFF_CLEAN" > "$OUT_DIR/${SAMPLE}_counts.txt"
done

