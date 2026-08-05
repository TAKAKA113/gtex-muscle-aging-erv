#!/bin/bash

#SBATCH --job-name=sra_all
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=1-53%5
#SBATCH --output=logs/sra_all_%A_%a.out
#SBATCH --error=logs/sra_all_%A_%a.err

set -euo pipefail

module purge
module load bluebear
module load SRA-Toolkit/2.10.9-gompi-2020b

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
RUN_LIST="${PROJECT_ROOT}/00_metadata/run_accessions.txt"

cd "${PROJECT_ROOT}"

RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUN_LIST}")

if [[ -z "${RUN}" ]]; then
    echo "ERROR: No run accession found for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

FINAL_FASTQ="${PROJECT_ROOT}/01_raw_fastq/${RUN}.fastq.gz"

echo "Job started: $(date)"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Computer: $(hostname)"
echo "Run accession: ${RUN}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

WORKDIR=$(mktemp -d \
    "/scratch/${USER}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX")

cleanup() {
    rm -rf "${WORKDIR}"
}

trap cleanup EXIT

mkdir -p \
    "${WORKDIR}/sra" \
    "${WORKDIR}/fastq" \
    "${WORKDIR}/tmp"

echo "Temporary directory: ${WORKDIR}"

prefetch \
    --max-size u \
    --output-directory "${WORKDIR}/sra" \
    "${RUN}"

fasterq-dump \
    "${WORKDIR}/sra/${RUN}/${RUN}.sra" \
    --outdir "${WORKDIR}/fastq" \
    --temp "${WORKDIR}/tmp" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --progress

gzip -1 "${WORKDIR}/fastq/${RUN}.fastq"

cp \
    "${WORKDIR}/fastq/${RUN}.fastq.gz" \
    "${FINAL_FASTQ}"

echo "Final FASTQ:"
ls -lh "${FINAL_FASTQ}"

echo "Job finished: $(date)"
