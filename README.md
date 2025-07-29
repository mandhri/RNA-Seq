RNA‑Seq Pipeline

A simple, end‑to‑end bash workflow for bulk RNA‑Seq analysis: from raw FASTQ download through QC, trimming and transcript quantification with Kallisto.

Repository structure:


.
├── download_sra.sh            # download FASTQ files via ENA FTP
├── ena_ftp_links.txt          # list of SRR accession FTP links
├── extract_adapters.sh        # extract adapter statistics with FastQC
├── fastq_trimmer.sh           # trim adapters using Skewer
├── generate_manifests.sh      # build Kallisto manifest files
├── single_end_kallisto.sh     # Kallisto quantification for single‑end reads
├── paired_end_kallisto.sh     # Kallisto quantification for paired‑end reads

Prerequisites

bash (version 4 or higher)
FastQC
MultiQC
Skewer
Kallisto
git (for version control)

Ensure these are installed and available in your \$PATH.
