#!/bin/bash

#SBATCH --job-name=fc_strand_all
#SBATCH --account=zhoujz-gnn-chem-mixture
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --array=1-53%10
#SBATCH --output=logs/fc_strand_all_%A_%a.out
#SBATCH --error=logs/fc_strand_all_%A_%a.err

set -euo pipefail

# ============================================================
# Module
# ============================================================

module purge
module load bluebear
module load bear-apps/2022a/live
module load Subread/2.0.4-GCC-11.3.0

# ============================================================
# Project paths
# ============================================================

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

RUN_LIST="${PROJECT_ROOT}/00_metadata/run_accessions.txt"

ALIGN_DIR="${PROJECT_ROOT}/04_alignment/bam_coord"

GTF="${PROJECT_ROOT}/03_reference/gencode_v26/gencode.v26.primary_assembly.annotation.gtf"

OUT_DIR="${PROJECT_ROOT}/02_qc/strandedness_featurecounts"

cd "${PROJECT_ROOT}"

mkdir -p "${OUT_DIR}"

# ============================================================
# Get sample accession
# ============================================================

RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUN_LIST}")

if [[ -z "${RUN}" ]]; then
    echo "ERROR: No run accession found for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

BAM="${ALIGN_DIR}/${RUN}/${RUN}_Aligned.sortedByCoord.out.bam"

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: BAM missing or empty: ${BAM}" >&2
    exit 1
fi

if [[ ! -s "${GTF}" ]]; then
    echo "ERROR: GTF missing or empty: ${GTF}" >&2
    exit 1
fi

# ============================================================
# Output summary for this sample
# ============================================================

ASSIGNED_TABLE="${OUT_DIR}/${RUN}_strand_assigned_reads.tsv"

printf "sample\tstrand_setting\tassigned_reads\n" > "${ASSIGNED_TABLE}"

echo "Job started: $(date)"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${RUN}"
echo "BAM: ${BAM}"
echo "GTF: ${GTF}"

# ============================================================
# Test all three strandedness settings
#
# s0 = unstranded
# s1 = forward-stranded
# s2 = reverse-stranded
# ============================================================

for STRAND in 0 1 2
do
    OUT_FILE="${OUT_DIR}/${RUN}_s${STRAND}.txt"

    echo
    echo "Running featureCounts with -s ${STRAND}"

    featureCounts \
        -T "${SLURM_CPUS_PER_TASK}" \
        -a "${GTF}" \
        -o "${OUT_FILE}" \
        -t exon \
        -g gene_id \
        -s "${STRAND}" \
        "${BAM}"

    SUMMARY_FILE="${OUT_FILE}.summary"

    if [[ ! -s "${OUT_FILE}" ]]; then
        echo "ERROR: featureCounts output missing: ${OUT_FILE}" >&2
        exit 1
    fi

    if [[ ! -s "${SUMMARY_FILE}" ]]; then
        echo "ERROR: featureCounts summary missing: ${SUMMARY_FILE}" >&2
        exit 1
    fi

    ASSIGNED_READS=$(awk -F '\t' '$1 == "Assigned" {print $2}' "${SUMMARY_FILE}")

    if [[ -z "${ASSIGNED_READS}" ]]; then
        echo "ERROR: Could not extract Assigned reads from ${SUMMARY_FILE}" >&2
        exit 1
    fi

    printf "%s\t%s\t%s\n" \
        "${RUN}" \
        "${STRAND}" \
        "${ASSIGNED_READS}" \
        >> "${ASSIGNED_TABLE}"
done

# ============================================================
# Show result
# ============================================================

echo
echo "Assigned-read comparison:"
column -t -s $'\t' "${ASSIGNED_TABLE}"

echo
echo "Job finished: $(date)"
