#!/bin/bash
#SBATCH --job-name=RSEM_quant
#SBATCH --partition=general
#SBATCH --array=1-6
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=logs/rsem_%A_%a.out
#SBATCH --error=logs/rsem_%A_%a.err
#SBATCH --account=XXXXX (Partition Here)

# ============================================================
# RNA-seq Quantification (STAR + RSEM)
# ============================================================
# Performs alignment and quantification for each sample
# using a Slurm job array
# ============================================================

set -euo pipefail

module load rsem
module load star

# === Sample list ===
SAMPLES=(
  SRR19421304 SRR19421305 SRR19421306
  SRR19421307 SRR19421308 SRR19421309
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$IDX]}

# === Paths ===
FASTQ_DIR="/path/to/fastq"
OUT_DIR="${FASTQ_DIR}/rsem_output"
RSEM_REF="/path/to/GRCh38/RSEM_ref/GRCh38"

mkdir -p "${OUT_DIR}"

STAR_DIR=$(dirname "$(which STAR)")

echo "[$(date)] Starting RSEM for ${SAMPLE}"
echo "STAR path: ${STAR_DIR}"
echo "RSEM reference: ${RSEM_REF}"

rsem-calculate-expression \
  --paired-end \
  --star \
  --star-path "${STAR_DIR}" \
  --star-gzipped-read-file \
  --num-threads "${SLURM_CPUS_PER_TASK} \
  ${FASTQ_DIR}/${SAMPLE}_1.fastq.gz \
  ${FASTQ_DIR}/${SAMPLE}_2.fastq.gz \
  ${RSEM_REF} \
  ${OUT_DIR}/${SAMPLE}

echo "[$(date)] Finished RSEM for ${SAMPLE}"
