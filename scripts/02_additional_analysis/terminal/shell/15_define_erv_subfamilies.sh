#!/bin/bash

set -euo pipefail

PROJECT_ROOT="/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

ANNOTATION="${PROJECT_ROOT}/05_counts/tecount/gene_te_feature_annotation.tsv"

OUTPUT="${PROJECT_ROOT}/03_reference/te_annotation/erv_subfamily_whitelist.tsv"

SUMMARY="${PROJECT_ROOT}/03_reference/te_annotation/erv_subfamily_whitelist_summary.txt"

cd "${PROJECT_ROOT}"

if [[ ! -s "${ANNOTATION}" ]]; then
    echo "ERROR: Missing annotation file: ${ANNOTATION}" >&2
    exit 1
fi

# Column positions in gene_te_feature_annotation.tsv:
# 1 feature_id
# 2 feature_type
# 7 te_subfamily
# 8 te_family
# 9 te_class

awk -F '\t' '
BEGIN {
    OFS = "\t"
    print "feature_id", "te_subfamily", "te_family", "te_class", \
          "erv_definition", "include_primary", "include_broad"
}

NR > 1 && $2 == "TE" {
    if ($9 == "LTR" && ($8 == "ERV1" || $8 == "ERVK" || $8 == "ERVL")) {
        print $1, $7, $8, $9, "strict_ERV", "TRUE", "TRUE"
    }
    else if ($9 == "LTR" && $8 == "ERVL-MaLR") {
        print $1, $7, $8, $9, "MaLR_broad_only", "FALSE", "TRUE"
    }
}
' "${ANNOTATION}" > "${OUTPUT}"

STRICT_COUNT=$(
    awk -F '\t' '
        NR > 1 && $6 == "TRUE" {n++}
        END {print n + 0}
    ' "${OUTPUT}"
)

BROAD_COUNT=$(
    awk -F '\t' '
        NR > 1 && $7 == "TRUE" {n++}
        END {print n + 0}
    ' "${OUTPUT}"
)

if [[ "${STRICT_COUNT}" -ne 478 ]]; then
    echo "ERROR: Expected 478 strict ERV subfamilies, found ${STRICT_COUNT}" >&2
    exit 1
fi

if [[ "${BROAD_COUNT}" -ne 564 ]]; then
    echo "ERROR: Expected 564 broad ERV/LTR subfamilies, found ${BROAD_COUNT}" >&2
    exit 1
fi

{
    echo "File: ${OUTPUT}"
    echo "Created: $(date --iso-8601=seconds)"
    echo "Created by: scripts/shell/15_define_erv_subfamilies.sh"
    echo "Input: ${ANNOTATION}"
    echo
    echo "Primary strict ERV definition:"
    echo "  te_class = LTR"
    echo "  te_family = ERV1, ERVK, or ERVL"
    echo "  Subfamilies: ${STRICT_COUNT}"
    echo
    echo "Broad sensitivity definition:"
    echo "  strict ERV set plus ERVL-MaLR"
    echo "  Subfamilies: ${BROAD_COUNT}"
    echo
    echo "Family counts:"

    awk -F '\t' '
        NR > 1 {
            family[$3]++
        }
        END {
            for (x in family) {
                print "  " x ": " family[x]
            }
        }
    ' "${OUTPUT}" | sort

} > "${SUMMARY}"

md5sum "${OUTPUT}" > "${OUTPUT}.md5"

echo "ERV whitelist created"
echo "Output: ${OUTPUT}"
echo "Summary: ${SUMMARY}"
echo "Primary strict ERV subfamilies: ${STRICT_COUNT}"
echo "Broad ERV/LTR subfamilies: ${BROAD_COUNT}"
