# RNA-Seq

> Learning RNA-Seq from scratch using public FASTQ datasets and Kallisto quantification.

---

## 📋 Table of Contents

1. [Overview](#overview)  
2. [Project Structure](#project-structure)  
3. [Prerequisites](#prerequisites)  
4. [Data Download](#data-download)  
5. [Quality Control & Adapter Extraction](#quality-control--adapter-extraction)  
6. [Read Trimming](#read-trimming)  
7. [Quantification (Kallisto)](#quantification-kallisto)  
8. [Generate Manifests & Reports](#generate-manifests--reports)  
9. [Usage Examples](#usage-examples)  
10. [Contributing](#contributing)  
11. [License](#license)  

---

## 🔍 Overview

This repository demonstrates a simple end-to-end bulk RNA-Seq workflow:

1. **Download** raw FASTQ from ENA/SRA  
2. **Run** FastQC + MultiQC to assess read quality  
3. **Extract** adapter sequences  
4. **Trim** reads using `skewer`  
5. **Quantify** transcript abundance with Kallisto  
6. **Compile** manifests and summary reports  

---

## 🗂 Project Structure

.
├── download_sra.sh             # ENA/SRA bulk‐download helper
├── ena_ftp_links.txt           # List of SRR FTP URLs
├── extract_adapters.sh         # Pull adapter‐content module from FastQC
├── fastq_trimmer.sh            # Skewer trimming wrapper
├── paired_end_kallisto.sh      # Kallisto pipeline for paired‐end data
├── single_end_kallisto.sh      # Kallisto pipeline for single‐end data
├── generate_manifests.sh       # Summarise outputs and build sample table
└── SraRunTable.csv             # Sample metadata / phenotype table

---
## ⚙️ Prerequisites
Git (v2.20+), Bash shell
FastQC & MultiQC
Skewer (for adapter trimming)
Kallisto (v0.46+)


---
## ⬇️ Data Download
Edit ena_ftp_links.txt with your SRR URLs.

Run:
bash
./download_sra.sh ena_ftp_links.txt fastq_files/

This will fetch all .fastq.gz into fastq_files/.


---
## ✅ Quality Control & Adapter Extraction

bash
# Raw QC
fastqc -o qc_raw/ fastq_files/*.fastq.gz
multiqc -o qc_raw/ qc_raw/

# Extract adapter content
./extract_adapters.sh qc_raw/*_fastqc.zip qc_raw/

✂️ Read Trimming

./fastq_trimmer.sh \
  --in-dir fastq_files/ \
  --out-dir fastq_trimmed/ \
  --adapters qc_raw/extracted_adapters.fa

Then re-run FastQC/MultiQC on fastq_trimmed/ if you like.


---
## 🎯 Quantification (Kallisto)

Single-end:
bash
./single_end_kallisto.sh \
  --reads fastq_trimmed/SRRXXXXX.fastq.gz \
  --index ref/transcriptome.idx \
  --out-dir kallisto_se_results/

Paired-end:
bash
./paired_end_kallisto.sh \
  --reads1 fastq_trimmed/*_R1.fastq.gz \
  --reads2 fastq_trimmed/*_R2.fastq.gz \
  --index ref/transcriptome.idx \
  --out-dir kallisto_pe_results/


---
## 📊 Generate Manifests & Reports
Finally, summarise everything:

bash
./generate_manifests.sh \
  --qc-raw qc_raw/ \
  --qc-trimmed qc_trimmed/ \
  --kallisto-se kallisto_se_results/ \
  --kallisto-pe kallisto_pe_results/ \
  --manifest output_manifest.tsv








