#!/bin/bash

#SBATCH --job-name=fastqc_pilot
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fastqc_pilot_%j.out
#SBATCH --error=logs/fastqc_pilot_%j.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2019b/live
module load FastQC/0.11.9-Java-11

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
RUN="SRR13388735"

cd "${PROJECT_ROOT}"

echo "Job started: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Computer: $(hostname)"
echo "Input FASTQ: ${RUN}.fastq.gz"

fastqc \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --outdir "${PROJECT_ROOT}/02_qc/fastqc_raw" \
  "${PROJECT_ROOT}/01_raw_fastq/${RUN}.fastq.gz"

echo "FastQC output:"
ls -lh "${PROJECT_ROOT}/02_qc/fastqc_raw"

echo "Job finished: $(date)"
