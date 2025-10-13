#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# -------- config --------
export PATH="/mnt/vol1/kallisto:$PATH"   # optional if kallisto not in PATH
PROJ_DIR=/mnt/vol1/RNA-Seq

DATA_DIR="$PROJ_DIR/trimmed"                     # folder with paired trimmed FASTQs
RESULTS_DIR="$PROJ_DIR/kallisto_tutorial/results"
LOG_DIR="$RESULTS_DIR/logs"

INDEX_DIR="/mnt/vol1/Mouse_model_RNA_Seq/index"  # shared index location
IDX="$INDEX_DIR/mouse_transcriptome.k31.idx"     # use one consistent index path/name

THREADS=16

PAIRS_FILE="$RESULTS_DIR/filenames_pe.txt"       # two columns: R1<TAB>R2
MISSING_REPORT="$RESULTS_DIR/missing_pairs.txt"
SAMPLES_MANIFEST="$RESULTS_DIR/samples.txt"

mkdir -p "$INDEX_DIR" "$RESULTS_DIR" "$LOG_DIR"

# -------- sanity checks --------
command -v kallisto >/dev/null 2>&1 || { echo "ERROR: kallisto not found in PATH"; exit 1; }
[[ -d "$DATA_DIR" ]] || { echo "ERROR: DATA_DIR not found: $DATA_DIR"; exit 1; }

# -------- build index once (Ensembl r115 mouse cDNA) --------
if [[ ! -f "$IDX" ]]; then
  echo "[INFO] Kallisto index not found. Building at: $IDX"
  cd "$INDEX_DIR"
  if [[ ! -f Mus_musculus.GRCm39.cdna.all.fa.gz ]]; then
    wget -O Mus_musculus.GRCm39.cdna.all.fa.gz \
      https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
  fi
  gunzip -c Mus_musculus.GRCm39.cdna.all.fa.gz > Mus_musculus.GRCm39.cdna.all.fa
  kallisto index -i "$IDX" Mus_musculus.GRCm39.cdna.all.fa
  echo "[INFO] Index built."
else
  echo "[INFO] Using existing index: $IDX"
fi

# -------- preflight: build pairs list & validate --------
: > "$PAIRS_FILE"
: > "$MISSING_REPORT"

echo "[INFO] Scanning for paired trimmed FASTQs in: $DATA_DIR"
shopt -s nullglob

# Support common patterns:
#   *_R1_trimmed.fastq.gz   ↔ *_R2_trimmed.fastq.gz
#   *-R1_trimmed.fastq.gz   ↔ *-R2_trimmed.fastq.gz
R1_FILES=( "$DATA_DIR"/*_R1_trimmed.fastq.gz "$DATA_DIR"/*-R1_trimmed.fastq.gz )

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No R1 files found matching *_R1_trimmed.fastq.gz or *-R1_trimmed.fastq.gz in $DATA_DIR"
  exit 1
fi

MISSING=0
for R1 in "${R1_FILES[@]}"; do
  fname=$(basename "$R1")
  case "$fname" in
    *_R1_trimmed.fastq.gz) R2="${R1/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}" ;;
    *-R1_trimmed.fastq.gz) R2="${R1/-R1_trimmed.fastq.gz/-R2_trimmed.fastq.gz}" ;;
    *) R2="";;
  esac

  if [[ -n "$R2" && -f "$R2" ]]; then
    printf "%s\t%s\n" "$R1" "$R2" >> "$PAIRS_FILE"
  else
    echo "Missing mate for: $(basename "$R1")  (expected $(basename "${R2:-UNKNOWN}"))" >> "$MISSING_REPORT"
    MISSING=1
  fi
done

if [[ "$MISSING" -ne 0 ]]; then
  echo "[ERROR] Pairing check failed. Missing mates detected:"
  cat "$MISSING_REPORT"
  echo "[HINT] Fix file names/locations, then rerun. Aborting before quantification."
  exit 2
fi

PAIR_COUNT=$(wc -l < "$PAIRS_FILE" | tr -d ' ')
echo "[INFO] Pairing check OK. Found $PAIR_COUNT pairs."
cp "$PAIRS_FILE" "$RESULTS_DIR/filenames_pe.checked.txt"

# -------- quantify all pairs (basic options) --------
: > "$SAMPLES_MANIFEST"
echo "[INFO] Starting paired-end quantification with $THREADS threads…"

while IFS=$'\t' read -r FQZ1 FQZ2; do
  # derive sample name by stripping the R1 suffix
  base=$(basename "$FQZ1")
  sample=${base%_R1_trimmed.fastq.gz}
  if [[ "$sample" == "$base" ]]; then
    sample=${base%-R1_trimmed.fastq.gz}
  fi

  OUTDIR="$RESULTS_DIR/$sample"
  mkdir -p "$OUTDIR"
  echo "$sample" >> "$SAMPLES_MANIFEST"

  echo "[INFO] $sample"
  {
    echo "kallisto quant (paired-end)"
    echo "sample:  $sample"
    echo "R1:      $FQZ1"
    echo "R2:      $FQZ2"
    echo "threads: $THREADS"
    echo "index:   $IDX"
    date
    kallisto quant \
      -i "$IDX" \
      -o "$OUTDIR" \
      -t "$THREADS" \
      "$FQZ1" "$FQZ2"
    date
    echo "done $sample"
  } &> "$LOG_DIR/${sample}.kallisto.pe.log"

done < "$PAIRS_FILE"

echo "[INFO] Quantification complete."
echo "[INFO] Results in: $RESULTS_DIR"
echo "[INFO] Logs in:    $LOG_DIR"
echo "[INFO] Manifest:   $SAMPLES_MANIFEST"
