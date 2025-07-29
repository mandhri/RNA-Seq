#!/usr/bin/env bash
set -euo pipefail

#Defining variables
PROJ_DIR=/mnt/vol1/RNA-Seq/kallisto_tutorial
REF_DIR=$PROJ_DIR/references
DATA_DIR=$PROJ_DIR/fastq_data
RESULTS_DIR=$PROJ_DIR/results

#THREADS=8 means “use up to 8 cores at once
#Fewer threads= jobs run more slowly (because less parallel work), but will leave more spare CPU capacity if you are on a shared server.
#More threads= faster runtime (up to a point), but you need those cores free.
THREADS=8

mkdir -p $REF_DIR $DATA_DIR $RESULTS_DIR

#Conditional reference build
if [ ! -f $REF_DIR/mouse_transcriptome.idx ]; then 
echo " Kalisto index is not available,building the index"
cd $REF_DIR
[ ! -f mus_musculus.cdna.all.fa.gz ] && \
wget https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
#keeping the original .fa.gz around for safekeeping or future re‑indexing
gunzip -c mus_musculus.cdna.all.fa.gz \
> mus_musculus.cdna.all.fa
kallisto index -i mouse_transcriptome.idx mus_musculus.cdna.all.fa
else
echo "Kalisto index is already existing.Skipping the build"
fi

#Sample data download and processing

echo "Starting to download the relavant fastq file for processing and aligment"

SAMPLES="SRR4000502"
for SAMPLE in $SAMPLES; do
echo "processing fastq sample: $SAMPLE"
if [ ! -f $DATA_DIR/${SAMPLE}_1.fastq.gz ]; then 
echo "fastq file $SAMPLE  DOWNLOADING"

#Single‑end SRR
fasterq-dump --split-files --threads "$THREADS" -O "$DATA_DIR" "$SAMPLE" \
  && gzip -f "$DATA_DIR/${SAMPLE}"_{1,2}.fastq


fi 
done

#Quantify each sample with kallisto

echo "Starting quantification of samples with Kallisto"

for SAMPLE in $SAMPLES; do
echo "Quantifying the sample: $SAMPLE"

#Define the input fastq pathways
R1="$DATA_DIR/${SAMPLE}_1.fastq.gz"
R2="$DATA_DIR/${SAMPLE}_2.fastq.gz"

#Output folder for these samples

OUTDIR="$RESULTS_DIR/$SAMPLE"
mkdir -p "$OUTDIR"

#which index to use(i), # where to put the results(o), how many CPU threads(t) # the two mates ($R1,$R2)
kallisto quant \
-i "$REF_DIR/mouse_transcriptome.idx" \
-o "$OUTDIR" \
-t "$THREADS" \
"$R1" "$R2"


echo "finished quanification of $SAMPLE"

done

echo "Results of alignment stored in $RESULTS Directory folder"
