#!/usr/bin/env bash
set -euo pipefail

#Check whether jq for processing JSON data is installed before start

command -v jq >/dev/null 2>&1 || {
  echo "ERROR: This script requires jq for processing JSON data files" >&2
  exit 1
}


# where your Kallisto results live
RESULTS_DIR="/mnt/vol1/RNA-Seq/kallisto_tutorial/results"

# 1) Summarise pseudo-alignment metrics
ALIGN_SUM="${RESULTS_DIR}/alignment_summary.tsv"

echo -e "Sample\tProcessed\tPseudoaligned\tRate(%)" > "${ALIGN_SUM}"
for D in "${RESULTS_DIR}"/SRR*; do
  SAMPLE=$(basename "$D")
  INFO="$D/run_info.json"
  NPRO=$(jq '.n_processed'  "$INFO")
  NPSE=$(jq '.n_pseudoaligned' "$INFO")
  RATE=$(jq -r '.n_pseudoaligned / .n_processed * 100' "$INFO")
  printf "%s\t%d\t%d\t%.2f\n" \
         "$SAMPLE" "$NPRO" "$NPSE" "$RATE" \
    >> "${ALIGN_SUM}"
done
echo "Aligment summary file complete"

# 2) Generate sample-to-file manifest
MANIFEST="${RESULTS_DIR}/samples.txt"
echo -e "sample\tpath" > "${MANIFEST}"
for D in "${RESULTS_DIR}"/SRR*; do
  SAMPLE=$(basename "$D")
  printf "%s\t%s/abundance.tsv\n" \
         "$SAMPLE" "$D" \
    >> "${MANIFEST}"
done
echo "MANIFEST file complete"
