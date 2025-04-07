cd /proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/genomics_data/Nanopore

mkdir -p /home/edman/genomeAnalysis/data/raw_data/nanopore/

for file in *; do
  ln -s "$(pwd)/$file" /home/edman/genomeAnalysis/data/raw_data/nanopore/
done
