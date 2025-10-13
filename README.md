# RNA-Seq Analysis Workflow

A reproducible, modular workflow for bulk RNA-Seq analysis using public FASTQ data, FastQC/MultiQC, Skewer and Kallisto.


## Table of Contents

- [Overview](#overview)
- [Workflow Architecture](#workflow-architecture)
- [Prerequisites](#prerequisites)
- [Installation Requirements](#installation-requirements)
- [Project Structure](#project-structure)
- [Data Download](#data-download)
- [Quality Control](#quality-control)
- [Adapter Detection and Removal](#adapter-detection-and-removal)
- [Read Trimming](#read-trimming)
- [Transcript Quantification](#transcript-quantification)
- [Results Reporting](#results-reporting)
- [Usage Examples](#usage-examples)
- [Project Structure](#project-structure)
- [Contributing](#contributing)

## Overview

This repository demonstrates a complete end-to-end bulk RNA-Seq analysis workflow, encompassing data acquisition through to quantified transcript abundance measurements.

### Key Features

- **Automated data retrieval** from ENA/SRA databases
- **Comprehensive quality assessment** using FastQC and MultiQC
- **Intelligent adapter detection** from FastQC reports
- **Robust read trimming** using Skewer algorithms
- **Efficient quantification** via Kallisto pseudoalignment
- **Automated reporting** and manifest generation

## Workflow Architecture

The analysis pipeline follows a systematic approach:

1. **Data Acquisition Phase**: Raw FASTQ files are retrieved from public repositories
2. **Quality Assessment Phase**: Initial quality metrics are calculated and visualised
3. **Preprocessing Phase**: Adapters are detected and reads are trimmed accordingly
4. **Quantification Phase**: Transcript abundances are estimated using Kallisto
5. **Compilation Phase**: Results are summarised and manifests are generated

## Prerequisites

### System Requirements

- Unix-like operating system (Linux0
- Minimum 8GB RAM (16GB recommended for large datasets)
- Sufficient storage space (datasets can be substantial)
- Internet connectivity for data downloads

### Software Dependencies

The following software packages must be installed and accessible via PATH:

- **Git** (version 2.20 or later)
- **Bash shell** (version 4.0 or later)
- **FastQC** (latest stable version)
- **MultiQC** (latest stable version)
- **Skewer** (adapter trimming utility)
- **Kallisto** (version 0.46 or later)

## Installation Requirements

Detailed installation procedures for dependencies:

### FastQC Installation
```bash
sudo apt-get update
sudo apt-get install fastqc
```

### MultiQC Installation
```bash
pip install --user multiqc
```

### Skewer Installation
```bash
# Compile from source
git clone https://github.com/relipmoc/skewer.git
cd skewer
make
sudo make install
```

### Kallisto Installation
```bash
sudo apt-get install kallisto
```

## Project Structure

```
RNA-Seq single_end/
├── download_sra.sh             # ENA/SRA bulk download utility
├── extract_adapters.sh         # FastQC adapter extraction module
├── fastq_trimmer_single_end.sh # Skewer trimming wrapper script
├── single_end_kallisto.sh      # Single-end Kallisto pipeline
├── generate_manifests.sh       # Output summarisation and manifest builder
├── SraRunTable.csv             # Sample metadata and phenotype information
└── README.md                   # Project documentation
```

## Step 1: Data Download 

### Preparing Download Links

The data download process has been automated to retrieve all FASTQ files associated with a specific ENA study accession. The download script performs the following operations:

1. **API Query**: ENA's file report API is queried to retrieve all FASTQ download URLs
2. **URL Processing**: FTP URLs are extracted and formatted appropriately
3. **Parallel Download**: Multiple files are downloaded concurrently using wget

### Configuration Parameters

The download script utilises several configurable parameters:

```bash
ENA_STUDY="SRP080947"                           # Target ENA study accession
OUTDIR="/mnt/vol1/RNA-Seq/fastq_files"         # Output directory path
LINKFILE="/mnt/vol1/RNA-Seq/ena_ftp_links.txt" # Generated URL list file
N_CORES=6                                       # Parallel download threads
```

### Executing Download

```bash
./download_sra.sh
```

The script automatically:
- Creates the output directory structure
- Queries ENA API for all FASTQ files in the specified study
- Generates `ena_ftp_links.txt` containing all download URLs
- Downloads files in parallel using the specified number of cores
- Provides progress feedback and error handling

### API Integration Details

The ENA API query format used:
```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_STUDY}&result=read_run&fields=fastq_ftp&format=tsv&download=true&limit=0"
```

This query retrieves tab-separated values containing FTP URLs for all FASTQ files associated with the study accession. The response is processed to extract individual URLs and formatted for wget compatibility.


## Step 1.1: Quality Assessment

### Initial Quality Control

Quality metrics are generated for all raw FASTQ files:

```bash
# Generate individual FastQC reports
fastqc -o qc_raw/ fastq_files/*.fastq.gz

# Compile aggregated MultiQC report
multiqc -o qc_raw/ qc_raw/
```

Quality reports will be available in the `qc_raw/` directory for examination.

## Step 2: Adapter Removal


### Read Trimming Execution
Most trimmers detect adapters for single-end reads. 

```bash
# Single-end trimming with fastp (auto-detect adapters)
./fastq_trimmer.sh \
  --in-dir fastq_files/ \
  --out-dir fastq_trimmed/
```

Trimmed reads will be written to the `fastq_trimmed/` directory.

### Step 2.1: Post-trimming Quality Assessment

```bash
# Optional: Re-assess quality after trimming
fastqc -o qc_trimmed/ fastq_trimmed/*.fastq.gz
multiqc -o qc_trimmed/ qc_trimmed/
```

## Step 3: Transcript Quantification

### Reference Preparation (transcriptome)

Kallisto builds an index from a **transcriptome FASTA** (e.g., Ensembl cDNA), not from a genome FASTA.
Use a transcriptome that matches your GTF **and release version**.  
In this project we use **mouse (Mus musculus), GRCm39, Ensembl r115**.

```bash
# Create an index once (example paths)
PROJ_DIR=/mnt/vol1/Mouse_model_RNA_Seq
INDEX_DIR="$PROJ_DIR/index"
mkdir -p "$INDEX_DIR"
cd "$INDEX_DIR"

# Download Ensembl r115 mouse cDNA transcriptome
wget -O Mus_musculus.GRCm39.cdna.all.fa.gz \
  https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

# Build Kallisto index from cDNA (default k=31)
gunzip -c Mus_musculus.GRCm39.cdna.all.fa.gz > transcripts.fa
kallisto index -i mouse_transcriptome.r115.k31.idx transcripts.fa

```

After preparing the index, quantify your samples with the relevant index and reads.

### Single-End Quantification

Kallisto requires the fragment length mean (-l) and sd (-s) for single-end.

```bash
kallisto quant --single \
  -i "$INDEX_DIR/mouse_transcriptome.r115.k31.idx" \
  -o /mnt/vol1/Mouse_model_RNA_Seq/kallisto_results/SAMPLE_ID_SE \
  -t 16 \
  -l 200 -s 20 \
  /path/to/trimmed/SAMPLE_ID_trimmed.fastq.gz
```


## Step 4: Results Compilation

### Step 4.1: Manifest Generation

Final summarisation and manifest creation:

```bash
./generate_manifests.sh \
  --qc-raw qc_raw/ \
  --qc-trimmed qc_trimmed/ \
  --kallisto-se kallisto_se_results/ \
  --kallisto-pe kallisto_pe_results/ \
  --manifest output_manifest.tsv
```

The output manifest (`output_manifest.tsv`) provides a comprehensive summary of all processed samples.

## Usage Examples

### Complete Workflow Example

```bash
# 1) Data acquisition
bash ./download_sra.sh ena_ftp_links.txt fastq_files/

# 2) Initial quality assessment
mkdir -p qc_raw
fastqc -o qc_raw/ fastq_files/*.fastq.gz
multiqc -o qc_raw/ qc_raw/

# 3) Adapter trimming

./fastq_trimmer.sh \
  --in-dir fastq_files/ \
  --out-dir fastq_trimmed/ \
  --adapters qc_raw/extracted_adapters.fa \
  --single-end

# 4) Post-trimming quality check
mkdir -p qc_trimmed
fastqc -o qc_trimmed/ fastq_trimmed/*.fastq.gz
multiqc -o qc_trimmed/ qc_trimmed/

# 5) Quantification (SINGLE-END requires -l and -s)
READS_DIR=fastq_trimmed
OUT_BASE=kallisto_se_results
INDEX=index/mouse_transcriptome.r115.k31.idx
FRAG_MEAN=200
FRAG_SD=20
mkdir -p "$OUT_BASE"

for FQ in "$READS_DIR"/*.fastq.gz; do
  sample=$(basename "$FQ" .fastq.gz)
  out="$OUT_BASE/$sample"; mkdir -p "$out"
  kallisto quant --single -l "$FRAG_MEAN" -s "$FRAG_SD" \
    -i "$INDEX" -o "$out" -t 16 "$FQ"
done

# 6) Results compilation
./generate_manifests.sh \
  --qc-raw qc_raw/ \
  --qc-trimmed qc_trimmed/ \
  --kallisto-se kallisto_se_results/ \
  --manifest final_manifest.tsv
```

## Output Interpretation

### Quality Control Outputs

- **FastQC reports**: Per-sample quality metrics including base quality, GC content, and adapter contamination
- **MultiQC reports**: Aggregated quality summaries across all samples

### Quantification Outputs

- **abundance.h5**: Kallisto's native HDF5 format containing quantification results
- **abundance.tsv**: Tab-separated quantification values (transcript ID, length, estimated counts, TPM)
- **run_info.json**: Quantification run metadata and statistics

### Manifest Contents

The final manifest includes:
- Sample identifiers and metadata
- Quality metrics pre- and post-trimming
- Quantification statistics
- File paths and checksums


### Performance Optimisation

- Utilise multiple threads where supported (`-t` parameter in Kallisto)
- Consider using compressed intermediate files to conserve storage
- Implement parallel processing for multiple samples

---

**Note**: For production-scale RNA-Seq analyses, additional considerations regarding statistical analysis, batch effects, and experimental design should be incorporated.
