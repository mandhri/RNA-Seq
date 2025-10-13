#!/usr/bin/env bash
set -euo pipefail

export PATH="/mnt/vol1/kallisto:$PATH"
PROJ_DIR=/mnt/vol1/RNA-Seq

#Where to build/index and store results
REF_DIR=$PROJ_DIR/kallisto_tutorial/references
RESULTS_DIR=$PROJ_DIR/kallisto_tutorial/results

#Point here at your already-downloaded trimmed FASTQs
DATA_DIR=$PROJ_DIR/trimmed

# ---- index lives in the Mouse_model_RNA_Seq area
INDEX_DIR="/mnt/vol1/Mouse_model_RNA_Seq/index"
IDX="$INDEX_DIR/mouse_transcriptome.k31.idx"

#THREADS=8 means “use up to 8 cores at once
#Fewer threads= jobs run more slowly (because less parallel work), but will leave more spare CPU capacity if you are on a shared server.
#More threads= faster runtime (up to a point), but you need those cores free.
THREADS=8

#IN Single-end (SE):, you must provide both -l (FRAG_LEN, mean fragment length) and -s (FRAG_SD, sd).
# IN case you dont know the -l and -s, then;
#A safe starting point for standard Illumina RNA-seq: -l 200 -s 30 (or -s 40–50 if your library is known to be more variable).

FRAG_LEN=200 
FRAG_SD=20

mkdir -p $REF_DIR $RESULTS_DIR



#Conditional reference build

if [ ! -f $REF_DIR/mouse_transcriptome.idx ]; then
echo "Kallisto index is not available, building the index"
cd $REF_DIR
[ ! -f Mus_musculus.GRCm39.cdna.all.fa.gz ] && \

wget https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

#keeping the original .fa.gz around for safekeeping or future re‑indexing
gunzip -c Mus_musculus.GRCm39.cdna.all.fa.gz > Mus_musculus.GRCm39.cdna.all.fa
kallisto index -i mouse_transcriptome.idx Mus_musculus.GRCm39.cdna.all.fa

else
echo "Kalisto index is already existing.Skipping the build"
fi

#Sample data download and processing

echo "Starting to download the relavant fastq file for processing and aligment"

# Check if any files exist
if ! ls "$DATA_DIR"/*-trimmed.fastq.gz 1> /dev/null 2>&1; then
    echo "ERROR: No trimmed.fastq.gz files found in $DATA_DIR"
    exit 1
fi

TOTAL_FILES=$(ls "$DATA_DIR"/*-trimmed.fastq.gz | wc -l)
CURRENT=0

echo "Found $TOTAL_FILES FASTQ files to process"

for FASTQ in "$DATA_DIR"/*-trimmed.fastq.gz; do
CURRENT=$((CURRENT + 1))
SAMPLE=$(basename "$FASTQ" -trimmed.fastq.gz)
OUTDIR="$RESULTS_DIR/$SAMPLE"
echo "Processing $SAMPLE ($CURRENT/$TOTAL_FILES)"

mkdir -p "$OUTDIR"

echo "Quantifying $SAMPLE"
kallisto quant \
-i "$REF_DIR/mouse_transcriptome.idx" \
-o "$OUTDIR" \
-t "$THREADS" \
--single -l "$FRAG_LEN" -s "$FRAG_SD" \
"$FASTQ"

echo "finished quanification of $SAMPLE"

done

echo "Results of alignment stored in $RESULTS_DIR"
