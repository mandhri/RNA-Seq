#!/usr/bin/env bash
set -euo pipefail

# Paths & env
export PATH="/mnt/vol1/FastQC:$HOME/bin:$PATH"
export PATH="/mnt/vol1/jre11/bin:$PATH"
export PATH="$HOME/.local/bin:$PATH"

FASTQ_DIR="/mnt/vol1/RNA-Seq/fastq_files"     # contains *_R1.fastq.gz and *_R2.fastq.gz
TRIM_DIR="/mnt/vol1/RNA-Seq/trimmed"
QC_DIR="/mnt/vol1/RNA-Seq/qc_trimmed"
THREADS=16


# Skewer thresholds (your choices)
QCUT=20     # trim bases with Phred < 20 from 3' ends
LCUT=31     # drop reads shorter than 31 nt after trimming
FCUT=5      # max # of Ns allowed (kept from your script)

mkdir -p "$TRIM_DIR" "$QC_DIR"

# Loop over paired-end samples by discovering R1 files
shopt -s nullglob
for R1 in "$FASTQ_DIR"/*_R1.fastq.gz; do
  base="$(basename "$R1" _R1.fastq.gz)"                 # e.g., BB336
  R2="$FASTQ_DIR/${base}_R2.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "Skipping $base: missing mate file $R2"
    continue
  fi

  echo "Trimming paired-end sample $base …"

  # 1) Trim adapters & low-quality tails; discard very short reads
  # Skewer writes: ${TRIM_DIR}/${base}-trimmed-pair1.fastq and ...pair2.fastq
  skewer \
    -q "$QCUT" \
    -l "$LCUT" \
    -f "$FCUT" \
    -t "$THREADS" \
    -o "$TRIM_DIR/$base" \
    "$R1" "$R2"

  # 2) Gzip trimmed outputs
  gzip -f "$TRIM_DIR/${base}-trimmed-pair1.fastq"
  gzip -f "$TRIM_DIR/${base}-trimmed-pair2.fastq"

  # 3) Rename to your preferred convention: *_R{1,2}_trimmed.fastq.gz
  mv -f "$TRIM_DIR/${base}-trimmed-pair1.fastq.gz" "$TRIM_DIR/${base}_R1_trimmed.fastq.gz"
  mv -f "$TRIM_DIR/${base}-trimmed-pair2.fastq.gz" "$TRIM_DIR/${base}_R2_trimmed.fastq.gz"

  # 4) FastQC on the trimmed mates (optional: bump -t if you like)
  fastqc -t 4 -o "$QC_DIR" \
    "$TRIM_DIR/${base}_R1_trimmed.fastq.gz" \
    "$TRIM_DIR/${base}_R2_trimmed.fastq.gz"
done

# 5) Aggregate all QC reports
multiqc "$QC_DIR" -o "$QC_DIR"

echo "Paired trimming complete. Trimmed FASTQs in $TRIM_DIR  |  QC in $QC_DIR"
