# GATK3 Snakemake Pipeline (Legacy)

A reproducible Snakemake workflow for the **legacy GATK3 preprocessing path**, implemented in a modular Snakemake project (WSL2/Ubuntu-friendly).

This repository **does NOT implement the GATK4 Best Practices pipeline**. It intentionally **stops after**:

✅ **Create Realignment Targets (GATK3 RealignerTargetCreator)**

---

## 📌 Aim of the Project

This workflow automates the *legacy* preprocessing steps up to realignment targets:

- Download and prepare a reference (Ensembl chr15 in this repo)
- Download FASTQs from SRA based on `accessions.txt`
- Align reads (BWA)
- Sort BAM
- Collect basic metrics
- Mark duplicates (Picard)
- Index dedup BAM (Picard)
- **Create realignment targets** (GATK3 RealignerTargetCreator)

---

## ⚠️ Important: GATK3 Tarball Required (Manual Step)

Because GATK3 is **license-restricted**, Bioconda may not ship the `GenomeAnalysisTK.jar`.  
This workflow expects you to provide the Broad tarball:

- **GenomeAnalysisTK-3.8-0.tar.bz2**
[Link](https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2)

### Place the file here (do not commit it)
resources/gatk/GenomeAnalysisTK-3.8-0.tar.bz2

Notes:
- The workflow extracts the tarball and caches the jar under `.snakemake/` for reuse.
- Ensure `resources/gatk/*.tar.bz2` is ignored in `.gitignore`.

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
├── samples.tsv  
├── envs/  
│ ├── bwa.yml  
│ ├── gatk.yml  
│ ├── gatk3.yml  
│ ├── picard.yml  
│ ├── reference.yml  
│ ├── snpeff.yml  
│ └── sra-tools.yml  
├── scripts/  
│ └── download_fastq.sh  
├── rules/  
│ ├── align.smk  
│ ├── download_fastq.smk  
│ ├── index_bam.smk  
│ ├── mark_duplicates.smk  
│ ├── metrics.smk  
│ ├── realigner_target_creator.smk  
│ ├── reference.smk  
│ └── sort_bam.smk  
├── reference/  
│ ├── chr_15.fa  
│ ├── chr_15.dict  
│ ├── chr_15.fa.fai  
│ └── (bwa index files)  
├── fastq/  
│ ├── SRRxxxxxx_1.fastq.gz  
│ └── SRRxxxxxx_2.fastq.gz  
├── results/  
│ ├── bam/  
│ ├── metrics/  
│ ├── dedup/  
│ └── realign/  
└── resources/  
└── gatk/  
└── GenomeAnalysisTK-3.8-0.tar.bz2 # manual, ignored by git  

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

### 4️⃣ Create a Snakemake “driver” environment

This env is only to run Snakemake. Rule-specific envs come from envs/*.yml via --use-conda.

conda install -n base -c conda-forge mamba -y  
mamba create -n snakemake -c conda-forge -c bioconda -y snakemake  
conda activate snakemake  

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

## ⚙ Configuration
Edit config.yaml as needed. Key entries:  
- accessions_file: list of SRA accessions  
- samples_tsv: mapping of sample IDs to FASTQ paths  
- reference_dir, reference_canonical  
- metrics_dir, bam_dir, dedup_dir, realign_dir  
- gatk3_tarball: path to the manual tarball (see above)  

---

## 🧾 FASTQ Download
`accessions.txt` example:  
SRR2584863  
SRR2584866  
SRR2584868  

FASTQs are downloaded by scripts/download_fastq.sh via the Snakemake rule.  

---

## ▶️ Running the Workflow
### Dry run (validate DAG)  
snakemake -n -p --use-conda --cores 1

### Run a single realignment targets output (recommended first)  
snakemake -p --use-conda --cores 4 --rerun-incomplete --keep-going \  
&nbsp;&nbsp;&nbsp;&nbsp;results/realign/SRR2584863.realignment_targets.list

### Run all targets defined in rule all
snakemake -p --use-conda --cores 4 --rerun-incomplete --keep-going

---

## 🧬 Workflow Stages (Legacy GATK3)
### ✔ Implemented in this branch
- Reference download + canonicalization + indices (chr15)
- FASTQ download from SRA (from accessions.txt)
- Alignment (BWA)
- Sort BAM
- Collect metrics
- Mark duplicates
- Build BAM index
- **Create Realignment Targets (GATK3 RealignerTargetCreator)**

### 🛑 Stopping point

This branch intentionally stops after RealignerTargetCreator.  
GATK4 Best Practices will be implemented on a separate branch based off the last shared preprocessing step.

---

## 👤 Contact
For issues or contributions, open a GitHub issue or pull request.
remco.crefcoeur@students.fhnw.ch
---
