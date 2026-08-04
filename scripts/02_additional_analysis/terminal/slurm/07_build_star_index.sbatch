#!/bin/bash

#SBATCH --job-name=star_index
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/star_index_%j.out
#SBATCH --error=logs/star_index_%j.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a/live
module load STAR/2.7.11b-GCC-13.3.0

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

REF_DIR="${PROJECT_ROOT}/03_reference/gencode_v26"
FASTA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
GTF="${REF_DIR}/gencode.v26.primary_assembly.annotation.gtf"

INDEX_DIR="${PROJECT_ROOT}/03_reference/star_index_gencode_v26_sjdb140"

cd "${PROJECT_ROOT}"

if [[ ! -s "${FASTA}" ]]; then
    echo "ERROR: FASTA is missing or empty: ${FASTA}" >&2
    exit 1
fi

if [[ ! -s "${GTF}" ]]; then
    echo "ERROR: GTF is missing or empty: ${GTF}" >&2
    exit 1
fi

mkdir -p "${INDEX_DIR}"

echo "Job started: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Computer: $(hostname)"
echo "STAR version: $(STAR --version)"
echo "Genome FASTA: ${FASTA}"
echo "Gene annotation: ${GTF}"
echo "Index directory: ${INDEX_DIR}"
echo "sjdbOverhang: 140"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

STAR \
    --runMode genomeGenerate \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${INDEX_DIR}" \
    --genomeFastaFiles "${FASTA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang 140

echo "STAR index files:"
ls -lh "${INDEX_DIR}"

echo "Job finished: $(date)"
