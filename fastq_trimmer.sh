#!/usr/bin/env bash
set -euo pipefail
export PATH="/mnt/vol1/FastQC:$HOME/bin:$PATH"
export PATH="/mnt/vol1/jre11/bin:$PATH"
export PATH="$HOME/.local/bin:$PATH"

FASTQ_DIR="/mnt/vol1/RNA-Seq/fastq_files"
TRIM_DIR="/mnt/vol1/RNA-Seq/trimmed"
QC_DIR="/mnt/vol1/RNA-Seq/qc_trimmed"
THREADS=16

#QCUT=20,trim bases below Phred 20
#LCUT=31, discard any reads < 31 nt so they’ll still contain 31-mers

QCUT=20
LCUT=31
FCUT=5
mkdir -p "$TRIM_DIR" "$QC_DIR"

for fq in "$FASTQ_DIR"/*.fastq.gz; do
sample=$(basename "$fq" .fastq.gz)
echo "Trimming $sample …"

#1.Trim adapters + low‑quality tails, drop reads < LCUT
skewer \
-q $QCUT \
-l $LCUT \
-f $FCUT \
-t $THREADS \
-o "$TRIM_DIR/$sample" \
"$fq"

# Skewer outputs:
  #$TRIM_DIR/${sample}-trimmed.fastq
  #$TRIM_DIR/${sample}-untrimmed.fastq

#2.Gzip the trimmed FASTQ
gzip -f "$TRIM_DIR/${sample}-trimmed.fastq"

#3.QC the trimmed reads
fastqc -t 4 \
-o "$QC_DIR" \
"$TRIM_DIR/${sample}-trimmed.fastq.gz"
done

#4.Aggregate all QC reports
multiqc "$QC_DIR" -o "$QC_DIR"

echo "Trimming complete: kept reads ≥${LCUT} nt. QC reports in $QC_DIR."
