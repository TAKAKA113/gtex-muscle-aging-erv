#!/bin/bash

#SBATCH --job-name=star_all
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --array=1-53%3
#SBATCH --output=logs/star_all_%A_%a.out
#SBATCH --error=logs/star_all_%A_%a.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a/live
module load STAR/2.7.11b-GCC-13.3.0

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

RUN_LIST="${PROJECT_ROOT}/00_metadata/run_accessions.txt"
INDEX_DIR="${PROJECT_ROOT}/03_reference/star_index_gencode_v26_sjdb140"

cd "${PROJECT_ROOT}"

RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUN_LIST}")

if [[ -z "${RUN}" ]]; then
    echo "ERROR: No run accession found for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

FASTQ="${PROJECT_ROOT}/01_raw_fastq/${RUN}.fastq.gz"
FINAL_DIR="${PROJECT_ROOT}/04_alignment/bam_coord/${RUN}"
FINAL_BAM="${FINAL_DIR}/${RUN}_Aligned.sortedByCoord.out.bam"
DONE_MARKER="${FINAL_DIR}/.complete"

if [[ ! -s "${FASTQ}" ]]; then
    echo "ERROR: FASTQ is missing or empty: ${FASTQ}" >&2
    exit 1
fi

if [[ ! -s "${INDEX_DIR}/Genome" ]]; then
    echo "ERROR: STAR index is missing: ${INDEX_DIR}" >&2
    exit 1
fi

if [[ -f "${DONE_MARKER}" && -s "${FINAL_BAM}" ]]; then
    echo "Completed alignment already exists. Skipping:"
    ls -lh "${FINAL_BAM}"
    exit 0
fi

mkdir -p "${FINAL_DIR}"

WORKDIR=$(mktemp -d \
    "/scratch/${USER}_star_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX")

cleanup() {
    rm -rf "${WORKDIR}"
}

trap cleanup EXIT

SCRATCH_OUT="${WORKDIR}/out"
PREFIX="${SCRATCH_OUT}/${RUN}_"

mkdir -p "${SCRATCH_OUT}"

echo "Job started: $(date)"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Computer: $(hostname)"
echo "STAR version: $(STAR --version)"
echo "Run accession: ${RUN}"
echo "Input FASTQ: ${FASTQ}"
echo "STAR index: ${INDEX_DIR}"
echo "Temporary directory: ${WORKDIR}"
echo "Final directory: ${FINAL_DIR}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

STAR \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${INDEX_DIR}" \
    --readFilesIn "${FASTQ}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${PREFIX}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM MD \
    --outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 200 \
    --limitBAMsortRAM 30000000000

SCRATCH_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
SCRATCH_LOG="${PREFIX}Log.final.out"

if [[ ! -s "${SCRATCH_BAM}" ]]; then
    echo "ERROR: STAR BAM was not created: ${SCRATCH_BAM}" >&2
    exit 1
fi

if [[ ! -s "${SCRATCH_LOG}" ]]; then
    echo "ERROR: STAR final log was not created: ${SCRATCH_LOG}" >&2
    exit 1
fi

rm -f "${FINAL_DIR}/${RUN}_"*
cp -a "${SCRATCH_OUT}/." "${FINAL_DIR}/"

touch "${DONE_MARKER}"

echo "Final alignment outputs:"
ls -lh "${FINAL_DIR}"

echo "Job finished: $(date)"
