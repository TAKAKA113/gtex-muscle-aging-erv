#!/bin/bash

#SBATCH --job-name=bam_qc
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-53%5
#SBATCH --output=logs/bam_qc_%A_%a.out
#SBATCH --error=logs/bam_qc_%A_%a.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a/live
module load SAMtools

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

RUN_LIST="${PROJECT_ROOT}/00_metadata/run_accessions.txt"
ALIGN_DIR="${PROJECT_ROOT}/04_alignment/bam_coord"
QC_DIR="${PROJECT_ROOT}/02_qc/samtools"

cd "${PROJECT_ROOT}"

RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUN_LIST}")

if [[ -z "${RUN}" ]]; then
    echo "ERROR: No run accession for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

BAM="${ALIGN_DIR}/${RUN}/${RUN}_Aligned.sortedByCoord.out.bam"
BAI="${BAM}.bai"
FLAGSTAT="${QC_DIR}/${RUN}.flagstat.txt"

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: BAM missing or empty: ${BAM}" >&2
    exit 1
fi

mkdir -p "${QC_DIR}"

echo "Job started: $(date)"
echo "Run: ${RUN}"
echo "BAM: ${BAM}"
echo "Samtools version:"
samtools --version | head -n 2

echo "Checking BAM structure..."
samtools quickcheck -v "${BAM}"

echo "Creating BAM index..."
samtools index \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${BAM}"

if [[ ! -s "${BAI}" ]]; then
    echo "ERROR: BAM index was not created: ${BAI}" >&2
    exit 1
fi

echo "Creating flagstat report..."
samtools flagstat \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${BAM}" \
    > "${FLAGSTAT}"

if [[ ! -s "${FLAGSTAT}" ]]; then
    echo "ERROR: flagstat report was not created: ${FLAGSTAT}" >&2
    exit 1
fi

echo "Outputs:"
ls -lh "${BAI}" "${FLAGSTAT}"

echo "Job finished: $(date)"
