#!/bin/bash

#SBATCH --job-name=multiqc_star
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/multiqc_star_%j.out
#SBATCH --error=logs/multiqc_star_%j.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a/live
module load MultiQC/1.32-foss-2024a

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

INPUT_DIR="${PROJECT_ROOT}/04_alignment/bam_coord"
OUTPUT_DIR="${PROJECT_ROOT}/02_qc/multiqc_star"

cd "${PROJECT_ROOT}"

mkdir -p "${OUTPUT_DIR}"

echo "Job started: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Computer: $(hostname)"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"

multiqc \
    "${INPUT_DIR}" \
    --module star \
    --outdir "${OUTPUT_DIR}" \
    --filename "GSE164471_STAR_alignment_multiqc.html" \
    --title "GSE164471 STAR alignment quality control"

echo "MultiQC output:"
ls -lh "${OUTPUT_DIR}"

echo "Job finished: $(date)"
