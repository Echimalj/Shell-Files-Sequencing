#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --partition=general
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/STAR_index_%j.out
#SBATCH --error=logs/STAR_index_%j.err
#SBATCH --account=XXXXXX (Partition Here)

# ============================================================
# STAR Genome Index Generation (GRCh38)
# ============================================================
# Builds STAR index for downstream RNA-seq alignment
# Requires: genome FASTA + GTF annotation
# ============================================================

set -euo pipefail

module load star

# === Reference paths ===
GENOME_FASTA="/path/to/GRCh38/fasta/genome.fa"
GTF_FILE="/path/to/GRCh38/genes/genes.gtf"
GENOME_DIR="/path/to/GRCh38/STAR_index"

# Create output directory
mkdir -p "$GENOME_DIR"
rm -rf "${GENOME_DIR:?}"/*

echo "[$(date)] Starting STAR index generation..."

STAR \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --runMode genomeGenerate \
  --genomeDir "$GENOME_DIR" \
  --genomeFastaFiles "$GENOME_FASTA" \
  --sjdbGTFfile "$GTF_FILE" \
  --sjdbOverhang 100

echo "[$(date)] STAR index generation complete."
