#!/bin/bash

#SBATCH --job-name=sra_pilot
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/sra_pilot_%j.out
#SBATCH --error=logs/sra_pilot_%j.err

set -euo pipefail

module purge
module load bluebear
module load SRA-Toolkit/2.10.9-gompi-2020b

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
RUN="SRR13388735"

cd "${PROJECT_ROOT}"

echo "Job started: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Computer: $(hostname)"
echo "Run accession: ${RUN}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

WORKDIR=$(mktemp -d /scratch/${USER}_${SLURM_JOB_ID}.XXXXXX)

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
  "${PROJECT_ROOT}/01_raw_fastq/"

echo "Final FASTQ:"
ls -lh "${PROJECT_ROOT}/01_raw_fastq/${RUN}.fastq.gz"

rm -rf "${WORKDIR}"

echo "Job finished: $(date)"
