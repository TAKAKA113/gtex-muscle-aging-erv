#!/bin/bash
#SBATCH --job-name=install_telescope
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/install_telescope_%j.out
#SBATCH --error=logs/install_telescope_%j.err

set -euo pipefail

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

ENV_DIR="${PROJECT_ROOT}/software/conda_envs/telescope-1.0.4.1"
SPEC_DIR="${PROJECT_ROOT}/software/environment_specs"

mkdir -p "${PROJECT_ROOT}/software/conda_envs"
mkdir -p "${SPEC_DIR}"
mkdir -p "${PROJECT_ROOT}/logs"

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

eval "$(${EBROOTMINIFORGE3}/bin/conda shell.bash hook)"
source "${EBROOTMINIFORGE3}/etc/profile.d/mamba.sh"

# Package download cacheをcompute nodeの一時領域へ置く
export CONDA_PKGS_DIRS="${TMPDIR}/conda_pkgs"
mkdir -p "${CONDA_PKGS_DIRS}"

# Telescope専用environmentを作成
mamba create --yes \
    --prefix "${ENV_DIR}" \
    --override-channels \
    --channel conda-forge \
    --channel bioconda \
    --strict-channel-priority \
    telescope=1.0.4.1

# Installation確認
"${ENV_DIR}/bin/telescope" --help \
    > "${SPEC_DIR}/telescope_help.txt"

mamba list --prefix "${ENV_DIR}" \
    > "${SPEC_DIR}/telescope_environment_packages.txt"

mamba list --prefix "${ENV_DIR}" telescope \
    > "${SPEC_DIR}/telescope_version.txt"

echo "Telescope environment created"
echo "Environment: ${ENV_DIR}"
echo
cat "${SPEC_DIR}/telescope_version.txt"
