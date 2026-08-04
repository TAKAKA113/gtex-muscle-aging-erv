#!/bin/bash

#SBATCH --job-name=feature_annotation
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/feature_annotation_%j.out
#SBATCH --error=logs/feature_annotation_%j.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2023a/live
module load R/4.4.1-gfbf-2023a

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

cd "${PROJECT_ROOT}"

echo "Job started: $(date)"
echo "Rscript: $(which Rscript)"
Rscript --version

Rscript scripts/R/14_build_feature_annotation.R

echo "Job finished: $(date)"
