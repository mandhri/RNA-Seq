#!/usr/bin/env bash
set -euo pipefail

GTF_DIR="/mnt/vol1/RNA-Seq/kallisto_tutorial/references"
OUT="/mnt/vol1/RNA-Seq/kallisto_tutorial/results/tx2gene.tsv"

GTF_GZ="${GTF_DIR}/Mus_musculus.GRCm39.114.gtf.gz"
GTF_FILE="${GTF_DIR}/Mus_musculus.GRCm39.114.gtf"
GTF_URL="ftp://ftp.ensembl.org/pub/release-114/gtf/mus_musculus/Mus_musculus.GRCm39.114.gtf.gz"

if [[ ! -f "$GTF_FILE" ]]; then
    wget -O "$GTF_GZ" "$GTF_URL"
    gunzip -f "$GTF_GZ"
fi

echo "GTF file created"

grep -P "\texon\t" "$GTF_FILE" \
  | sed -nE 's/.*gene_id "([^"]+)";.*transcript_id "([^"]+)".*/\2\t\1/p' \
  | sort -u \
  > "$OUT"
echo "transcript to gene mapping completed"
