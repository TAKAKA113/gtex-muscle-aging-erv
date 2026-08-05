#!/bin/bash

#SBATCH --job-name=fastqc_all
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=1-53%10
#SBATCH --output=logs/fastqc_all_%A_%a.out
#SBATCH --error=logs/fastqc_all_%A_%a.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2019b/live
module load FastQC/0.11.9-Java-11

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
RUN_LIST="${PROJECT_ROOT}/00_metadata/run_accessions.txt"

cd "${PROJECT_ROOT}"

RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUN_LIST}")

if [[ -z "${RUN}" ]]; then
    echo "ERROR: No run accession found for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

FASTQ="${PROJECT_ROOT}/01_raw_fastq/${RUN}.fastq.gz"
OUTDIR="${PROJECT_ROOT}/02_qc/fastqc_raw"

if [[ ! -s "${FASTQ}" ]]; then
    echo "ERROR: FASTQ is missing or empty: ${FASTQ}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

echo "Job started: $(date)"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Run accession: ${RUN}"
echo "Input FASTQ: ${FASTQ}"
echo "Computer: $(hostname)"

fastqc \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --outdir "${OUTDIR}" \
    "${FASTQ}"

echo "FastQC results:"
ls -lh \
    "${OUTDIR}/${RUN}_fastqc.html" \
    "${OUTDIR}/${RUN}_fastqc.zip"

echo "Job finished: $(date)"
