# GATK4 Snakemake Pipeline

A reproducible Snakemake workflow for GATK4-based variant calling, including reference preparation, SRA FASTQ retrieval, alignment, variant calling, and modular rule organization.

---

## 📌 Aim of the Project
This workflow automates a complete GATK4 preprocessing and variant-calling pipeline. It supports:
- Automatic download of chr15 reference (FASTA + index + dict)
- Automatic FASTQ download from SRA based on a list of accessions
- Alignment and BAM processing
- Variant calling and annotation (future work)
- Modular Snakemake rule structure
- Reproducible conda environments

---

## 📁 Folder Structure
gatk4-snakemake/  
├── LICENSE  
├── README.md  
├── Snakefile  
├── config.yaml  
├── environment.yml  
├── .gitignore  
├── accessions.txt  
├── envs/  
│   ├── bwa.yml  
│   ├── gatk.yml  
│   ├── reference.yml  
│   ├── snpeff.yml  
│   └── sra-tools.yml  
├── scripts/  
│   └── download_fastq.sh  
├── rules/  
│   ├── align.smk  
│   ├── download_fastq.smk  
│   └── reference.smk  
├── reference/  
│   └── Homo_sapiens.GRCh38.dna.chromosome.15.fa  
├── fastq/  
│   ├── SRRxxxxxx_1.fastq.gz  
│   └── SRRxxxxxx_2.fastq.gz  
├── results/  
│   ├── bam/  
│   ├── vcfs/  
│   └── logs/  


---

## 🚀 Installing on a Clean WSL2 + Ubuntu System

### 1️⃣ Install WSL2 + Ubuntu
wsl --install

### 2️⃣ Clone the Repository Inside WSL
cd ~
git clone git@github.com:rcrefcoeur/gatk4-snakemake.git
cd gatk4-snakemake

### 3️⃣ Install Miniconda
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

### 4️⃣ Install Mamba + Create Workflow Environment
conda install -n base -c conda-forge mamba -y

mamba create -n gatk-tools -c conda-forge -c bioconda -y \
  bwa samtools sra-tools gatk4 picard snpeff snakemake wget curl python

conda activate gatk-tools

---

## 🔐 GitHub Authentication in WSL

### SSH (Recommended)
ssh-keygen -t ed25519 -C "your_email@example.com"
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
cat ~/.ssh/id_ed25519.pub

Add this key at: https://github.com/settings/ssh/new

Test:
ssh -T git@github.com

Clone via SSH:
git clone git@github.com:rcrefcoeur/gatk4-snakemake.git

---

## ⚙ Running the Workflow

Run entire pipeline:
snakemake --use-conda --cores 4

Dry run:
snakemake -n

Run specific target:
snakemake --use-conda reference/chr15.fa

---

## 🧬 Workflow Stages

### ✔ Completed
- Download chr15 reference from Ensembl
- Indexing + dictionary creation
- Download FASTQ from SRA using accessions.txt
- Auto-generate samples.tsv

### ⏳ To Be Implemented
- BWA-MEM2 alignment
- Sorting + indexing BAMs
- Duplicate marking (Picard)
- BQSR
- HaplotypeCaller (per-sample GVCF)
- Joint genotyping
- Variant filtering
- Annotation with snpEff
- Final reports + QC summaries

---

## 📝 Configuration File
Example `config.yaml`:

reference_url: "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz"
reference_dir: "reference"
reference_fa: "reference/Homo_sapiens.GRCh38.dna.chromosome.15.fa"

accessions_file: "accessions.txt"
samples_tsv: "samples.tsv"
fastq_dir: "fastq"

threads: 6
platform: "ILLUMINA"
outdir: "results"

---

## 🧾 FASTQ Download
`accessions.txt` example:
SRR2584863
SRR2584866
SRR2584868

`scripts/download_fastq.sh` runs automatically during:
snakemake --use-conda

---

## 👤 Contact
For issues or contributions, open a GitHub issue or pull request.
remco.crefcoeur@students.fhnw.ch
---
