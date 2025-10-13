#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ---- config ----
export PATH="/mnt/vol1/kallisto:$PATH"   # optional if kallisto isn't already in PATH
PROJ_DIR=/mnt/vol1/RNA-Seq

DATA_DIR="$PROJ_DIR/trimmed"                   # single-end trimmed FASTQs live here
RESULTS_DIR="$PROJ_DIR/kallisto_tutorial/results"
LOG_DIR="$RESULTS_DIR/logs"

INDEX_DIR="/mnt/vol1/Mouse_model_RNA_Seq/index" # shared index location
IDX="$INDEX_DIR/mouse_transcriptome.k31.idx"     # use one consistent index path/name

# ---- SE quant params ------------------

#THREADS=8 means “use up to 8 cores at once
#Fewer threads= jobs run more slowly (because less parallel work), but will leave more spare CPU capacity if you are on a shared server.
#More threads= faster runtime (up to a point), but you need those cores free.
THREADS=8

#IN Single-end (SE):, you must provide both -l (FRAG_LEN, mean fragment length) and -s (FRAG_SD, sd).
# IN case you dont know the -l and -s, then;
#A safe starting point for standard Illumina RNA-seq: -l 200 -s 30 (or -s 40–50 if your library is known to be more variable).

FRAG_LEN=200 
FRAG_SD=20

mkdir -p "$INDEX_DIR" "$RESULTS_DIR" "$LOG_DIR"

# ---- sanity checks ----
command -v kallisto >/dev/null 2>&1 || { echo "ERROR: kallisto not found in PATH"; exit 1; }
[[ -d "$DATA_DIR" ]] || { echo "ERROR: DATA_DIR not found: $DATA_DIR"; exit 1; }


#Conditional reference build

# ---- build index once (Ensembl r115 mouse cDNA) ----
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
  echo "[INFO] Kallisto index already exists: $IDX"
fi

# ---- find SE files (support *-trimmed.fastq.gz and *_trimmed.fastq.gz) ----
shopt -s nullglob
SE_FILES=( "$DATA_DIR"/*-trimmed.fastq.gz "$DATA_DIR"/*_trimmed.fastq.gz )
if [[ ${#SE_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No trimmed single-end FASTQs found in $DATA_DIR"
  echo "Expect patterns: *-trimmed.fastq.gz or *_trimmed.fastq.gz"
  exit 1
fi

echo "[INFO] Found ${#SE_FILES[@]} single-end FASTQ files to quantify."

# ---- quantify each file ----
count=0
for FASTQ in "${SE_FILES[@]}"; do
  ((count++))
  # sample name = filename without the trimmed suffix
  fname=$(basename "$FASTQ")
  sample=${fname%-trimmed.fastq.gz}
  if [[ "$sample" == "$fname" ]]; then
    sample=${fname%_trimmed.fastq.gz}
  fi

  OUTDIR="$RESULTS_DIR/$sample"
  mkdir -p "$OUTDIR"

  echo "[INFO] ($count/${#SE_FILES[@]}) Quantifying $sample"
  {
    echo "kallisto quant (single-end)"
    echo "sample:  $sample"
    echo "fastq:   $FASTQ"
    echo "threads: $THREADS"
    echo "index:   $IDX"
    echo "fraglen: $FRAG_LEN, s.d.: $FRAG_SD"
    date
    kallisto quant --single \
      -i "$IDX" \
      -o "$OUTDIR" \
      -t "$THREADS" \
      -l "$FRAG_LEN" -s "$FRAG_SD" \
      "$FASTQ"
    date
    echo "done $sample"
  } &> "$LOG_DIR/${sample}.kallisto.se.log"
done

echo "[INFO] Quantification complete. Results in: $RESULTS_DIR"
echo "[INFO] Logs in: $LOG_DIR"
