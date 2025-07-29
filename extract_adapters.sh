#!/usr/bin/env bash
set -euo pipefail

QC_DIR="/mnt/vol1/RNA-Seq/qc_raw"

for zip in "$QC_DIR"/*_fastqc.zip; do
#strip path, e.g. “SRR4000502_fastqc.zip
sample="${zip##*/}" 
#strip suffix, e.g. “SRR4000502”
sample="${sample%_fastqc.zip}"
echo "== $sample =="

  unzip -p "$zip" '*/fastqc_data.txt' \
    | sed -n '/^>>Adapter Content/,/^>>END_MODULE/p'
done
