#!/usr/bin/env bash
set -euo pipefail

ENA_STUDY="SRP080947"
OUTDIR="/mnt/vol1/RNA-Seq/fastq_files"
LINKFILE="/mnt/vol1/RNA-Seq/ena_ftp_links.txt"
N_CORES=4

mkdir -p "$OUTDIR"
echo "Output directory: $OUTDIR"

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_STUDY}&result=read_run&fields=fastq_ftp&format=tsv&download=true&limit=0" \
  | tail -n +2 | cut -f2 | tr ';' '\n' | sed '/^$/d' | sed 's|^|ftp://|' > "$LINKFILE"

if [[ ! -s "$LINKFILE" ]]; then
  echo "ERROR: No URLs found in $LINKFILE. Check ENA_STUDY value or network."
  exit 1
fi

echo "Found $(wc -l < "$LINKFILE") files to download."
echo "Downloading..."
xargs -a "$LINKFILE" -n1 -P"$N_CORES" wget -c -P "$OUTDIR"

echo "Done: FASTQs are in $OUTDIR"

