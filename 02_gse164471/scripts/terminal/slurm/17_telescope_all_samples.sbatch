#!/bin/bash
#SBATCH --job-name=telescope
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --array=0-52%8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/telescope_%A_%a.out
#SBATCH --error=logs/telescope_%A_%a.err

set -euo pipefail

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
SAMPLE_LIST="${PROJECT_ROOT}/00_metadata/run_accessions.txt"
GTF="${PROJECT_ROOT}/03_reference/telescope_annotation/HERV_rmsk.hg38.v2/transcripts.gtf"
TELESCOPE="${PROJECT_ROOT}/software/conda_envs/telescope-1.0.4.1/bin/telescope"
OUTPUT_ROOT="${PROJECT_ROOT}/05_counts/telescope"

cd "${PROJECT_ROOT}"

module purge
module load bluebear
module load SAMtools

mapfile -t SAMPLES < <(
    grep -v '^[[:space:]]*$' "${SAMPLE_LIST}" |
    tr -d '\r'
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

INPUT_BAM="${PROJECT_ROOT}/04_alignment/bam_coord/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
SAMPLE_OUT="${OUTPUT_ROOT}/${SAMPLE}"

if [[ ! -s "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM not found: ${INPUT_BAM}" >&2
    exit 1
fi

if [[ ! -s "${GTF}" ]]; then
    echo "ERROR: Telescope annotation not found: ${GTF}" >&2
    exit 1
fi

if [[ ! -x "${TELESCOPE}" ]]; then
    echo "ERROR: Telescope executable not found: ${TELESCOPE}" >&2
    exit 1
fi

if [[ -f "${SAMPLE_OUT}/.telescope.complete" ]]; then
    echo "Already completed: ${SAMPLE}"
    exit 0
fi

# 前回の不完全な出力だけを削除
rm -rf "${SAMPLE_OUT}"
mkdir -p "${SAMPLE_OUT}"

# BlueBEAR compute nodeのlocal /scratchを使用
BB_WORKDIR=$(
    mktemp -d \
    "/scratch/${USER}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX"
)

export TMPDIR="${BB_WORKDIR}"

QNAME_BAM="${BB_WORKDIR}/${SAMPLE}.qname.bam"
TELESCOPE_TMP="${BB_WORKDIR}/telescope_tmp"

mkdir -p "${TELESCOPE_TMP}"

cleanup() {
    test -d "${BB_WORKDIR}" && rm -rf "${BB_WORKDIR}"
}
trap cleanup EXIT

echo "Sample: ${SAMPLE}"
echo "Input BAM: ${INPUT_BAM}"
echo "Scratch directory: ${BB_WORKDIR}"
echo "Output directory: ${SAMPLE_OUT}"
echo "Scratch availability:"
df -h "${BB_WORKDIR}"
echo "Start: $(date --iso-8601=seconds)"

echo "Creating query-name-sorted BAM"

samtools sort \
    -n \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -m 4G \
    -T "${BB_WORKDIR}/${SAMPLE}.sorttmp" \
    -o "${QNAME_BAM}" \
    "${INPUT_BAM}"

samtools quickcheck -v "${QNAME_BAM}"

echo "Running Telescope"

"${TELESCOPE}" assign \
    "${QNAME_BAM}" \
    "${GTF}" \
    --attribute locus \
    --ncpu 1 \
    --tempdir "${TELESCOPE_TMP}" \
    --outdir "${SAMPLE_OUT}" \
    --exp_tag "${SAMPLE}" \
    --stranded_mode None \
    --reassign_mode exclude \
    --overlap_mode threshold \
    --overlap_threshold 0.2 \
    --annotation_class intervaltree

touch "${SAMPLE_OUT}/.telescope.complete"

echo "Telescope completed: ${SAMPLE}"
echo "End: $(date --iso-8601=seconds)"
