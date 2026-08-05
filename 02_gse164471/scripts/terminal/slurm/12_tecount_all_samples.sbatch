#!/bin/bash

#SBATCH --job-name=tecount_all
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-53%5
#SBATCH --output=logs/tecount_%A_%a.out
#SBATCH --error=logs/tecount_%A_%a.err

set -euo pipefail

# ============================================================
# Modules
# ============================================================

module purge
module load bluebear
module load bear-apps/2023a/live
module load TEtranscripts/2.2.3-foss-2023a

# ============================================================
# Project paths
# ============================================================

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

RUN_LIST="${PROJECT_ROOT}/00_metadata/run_accessions.txt"

ALIGN_DIR="${PROJECT_ROOT}/04_alignment/bam_coord"

GENE_GTF="${PROJECT_ROOT}/03_reference/gencode_v26/gencode.v26.primary_assembly.annotation.gtf"

TE_GTF="${PROJECT_ROOT}/03_reference/te_annotation/GRCh38_GENCODE_rmsk_TE.gtf"

OUT_DIR="${PROJECT_ROOT}/05_counts/tecount"

cd "${PROJECT_ROOT}"

mkdir -p "${OUT_DIR}"

# ============================================================
# Get sample accession
# ============================================================

RUN=$(
    sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUN_LIST}" |
    tr -d '\r'
)

if [[ -z "${RUN}" ]]; then
    echo "ERROR: No run accession for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

BAM="${ALIGN_DIR}/${RUN}/${RUN}_Aligned.sortedByCoord.out.bam"

SAMPLE_OUT="${OUT_DIR}/${RUN}"

COMPLETE_FILE="${SAMPLE_OUT}/${RUN}.tecount.complete"

# ============================================================
# Validate inputs
# ============================================================

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: BAM missing or empty: ${BAM}" >&2
    exit 1
fi

if [[ ! -s "${GENE_GTF}" ]]; then
    echo "ERROR: Gene GTF missing or empty: ${GENE_GTF}" >&2
    exit 1
fi

if [[ ! -s "${TE_GTF}" ]]; then
    echo "ERROR: TE GTF missing or empty: ${TE_GTF}" >&2
    exit 1
fi

mkdir -p "${SAMPLE_OUT}"

# Skip completed samples if the job is resubmitted
if [[ -s "${COMPLETE_FILE}" ]]; then
    echo "Already completed: ${RUN}"
    exit 0
fi

# Remove incomplete outputs from a previous failed run
find "${SAMPLE_OUT}" \
    -maxdepth 1 \
    -type f \
    ! -name "${RUN}.tecount.complete" \
    -delete

# ============================================================
# Run TEcount
# ============================================================

echo "Job started: $(date)"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${RUN}"
echo "BAM: ${BAM}"
echo "Gene GTF: ${GENE_GTF}"
echo "TE GTF: ${TE_GTF}"
echo "TEcount executable: $(which TEcount)"

TEcount \
    --BAM "${BAM}" \
    --GTF "${GENE_GTF}" \
    --TE "${TE_GTF}" \
    --format BAM \
    --stranded no \
    --mode multi \
    --project "${RUN}" \
    --outdir "${SAMPLE_OUT}" \
    --sortByPos

# ============================================================
# Validate outputs
# ============================================================

OUTPUT_COUNT=$(
    find "${SAMPLE_OUT}" \
        -maxdepth 1 \
        -type f \
        -size +0c |
    wc -l
)

if [[ "${OUTPUT_COUNT}" -eq 0 ]]; then
    echo "ERROR: TEcount produced no non-empty output files for ${RUN}" >&2
    exit 1
fi

touch "${COMPLETE_FILE}"

echo
echo "TEcount output files:"
find "${SAMPLE_OUT}" \
    -maxdepth 1 \
    -type f \
    -printf '%f\t%s bytes\n' |
    sort

echo
echo "Job finished: $(date)"
