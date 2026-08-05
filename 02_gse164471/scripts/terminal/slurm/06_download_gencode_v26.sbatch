#!/bin/bash

#SBATCH --job-name=download_ref
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/download_ref_%j.out
#SBATCH --error=logs/download_ref_%j.err

set -euo pipefail

module purge
module load bluebear

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
REF_DIR="${PROJECT_ROOT}/03_reference/gencode_v26"

FASTA_GZ="${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz"
FASTA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"

GTF_GZ="${REF_DIR}/gencode.v26.primary_assembly.annotation.gtf.gz"
GTF="${REF_DIR}/gencode.v26.primary_assembly.annotation.gtf"

FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.primary_assembly.annotation.gtf.gz"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "Job started: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Computer: $(hostname)"
echo "Reference directory: ${REF_DIR}"

echo "Downloading GRCh38 genome FASTA..."

curl \
    --fail \
    --location \
    --retry 5 \
    --retry-delay 10 \
    --output "${FASTA_GZ}.part" \
    "${FASTA_URL}"

gzip -t "${FASTA_GZ}.part"
mv "${FASTA_GZ}.part" "${FASTA_GZ}"

echo "Downloading GENCODE v26 gene annotation..."

curl \
    --fail \
    --location \
    --retry 5 \
    --retry-delay 10 \
    --output "${GTF_GZ}.part" \
    "${GTF_URL}"

gzip -t "${GTF_GZ}.part"
mv "${GTF_GZ}.part" "${GTF_GZ}"

echo "Decompressing FASTA..."

gzip -dc "${FASTA_GZ}" > "${FASTA}.part"
mv "${FASTA}.part" "${FASTA}"

echo "Decompressing GTF..."

gzip -dc "${GTF_GZ}" > "${GTF}.part"
mv "${GTF}.part" "${GTF}"

echo "Downloaded reference files:"
ls -lh "${FASTA_GZ}" "${FASTA}" "${GTF_GZ}" "${GTF}"

echo "First FASTA header:"
grep -m 1 '^>' "${FASTA}"

echo "First GTF annotation record:"
awk '!/^#/ {print; exit}' "${GTF}"

echo "Job finished: $(date)"
