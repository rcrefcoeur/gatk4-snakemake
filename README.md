# GATK4 Snakemake Pipeline (gatk4-mainline)

A reproducible Snakemake workflow that implements the **shared preprocessing baseline** for a GATK4-based variant calling pipeline:
**download → reference prep → FASTQ retrieval → alignment → sort → mark duplicates → BAM index**.

This branch intentionally **stops at**:

✅ `results/dedup/<SAMPLE>.dedup.bai`

GATK4 variant calling steps (BQSR, HaplotypeCaller, joint genotyping, filtering, annotation) will be added after this checkpoint.
---

## 📌 Project Goals

- Minimal required inputs: **`config.yaml`** + **`accessions.txt`**
- Fully reproducible toolchain via **Snakemake + per-rule conda environments**
- Modular rules under `rules/`

---

## ✅ What’s Implemented (Current Baseline)

For each accession / sample:

1. Download Ensembl chr15 reference (`.fa.gz` → `.fa`)
2. Create:
   - `.fai` (samtools faidx)
   - `.dict` (GATK CreateSequenceDictionary)
   - **BWA index** (bwa index) for canonical `reference/chr_15.fa`
3. Download FASTQs (SRA Toolkit)
4. Align reads (bwa mem)
5. Sort alignments (samtools/picard rule)
6. Mark duplicates (Picard)
7. Build BAM index (Picard)

---

## 🧾 accessions.txt (Comments Supported)

Lines starting with `#` are ignored, so you can keep optional accessions ready:

```text
SRR2584863
# SRR2584866
# SRR2584868

Uncomment later if time/disk allow.

---

## 📁 Folder Structure
gatk4-snakemake/  
├── Snakefile  
├── config.yaml  
├── accessions.txt  
├── scripts/  
│   └── download_fastq.sh  
├── envs/  
│   ├── bwa.yml  
│   ├── picard.yml  
│   ├── reference.yml  
│   ├── sra-tools.yml  
│   └── gatk.yml              # reserved for upcoming GATK4 steps  
├── rules/  
│   ├── download_fastq.smk  
│   ├── reference.smk  
│   ├── align.smk  
│   ├── sort_bam.smk  
│   ├── metrics.smk  
│   ├── mark_duplicates.smk  
│   └── index_bam.smk  
├── reference/                # generated  
├── fastq/                    # generated  
├── results/                  # generated  
│   ├── bam/  
│   ├── metrics/  
│   └── dedup/  
└── .snakemake/               # generated (conda env cache, logs, metadata)  

---

## 🚀 Installing on a Clean WSL2 + Ubuntu System

### 1️⃣ Install WSL2 + Ubuntu
wsl --install

### 2️⃣ Clone the Repository Inside WSL
cd ~
git clone git@github.com:rcrefcoeur/gatk4-snakemake.git
cd gatk4-snakemake
git checkout gatk4-mainline  

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

## ▶️ Running the Workflow
### Dry run (validate DAG)  
snakemake -n -p --use-conda --cores 1 results/dedup/SRR2584863.dedup.bai

### Run the baseline for one sample (recommended)  
snakemake -p --use-conda --cores 4 --rerun-incomplete --keep-going \  
&nbsp;&nbsp;&nbsp;&nbsp;results/dedup/SRR2584863.dedup.bai

### Run everything in rule all
snakemake -p --use-conda --cores 4 --rerun-incomplete --keep-going

---

## 👤 Contact
For issues or contributions, open a GitHub issue or pull request.
remco.crefcoeur@students.fhnw.ch
---