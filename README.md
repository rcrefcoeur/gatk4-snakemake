# GATK4 Snakemake Pipeline (gatk4-mainline)

A reproducible Snakemake workflow implementing a **GATK4-based pipeline up through first-pass variant calling**:

**download → reference prep → FASTQ retrieval → alignment → sort → mark duplicates → BAM index → HaplotypeCaller**

Minimal required inputs on a clean clone:

✅ `config.yaml` + `accessions.txt`

---

## 📌 Project Goals

- Minimal inputs: **`config.yaml`** + **`accessions.txt`**
- Reproducible execution via **Snakemake + per-rule conda envs** (`--use-conda`)
- Modular rule organization under `rules/`
- WSL2 + Ubuntu friendly

---

## ✅ Implemented Steps (Current)

For each accession / sample:

1. **Reference prep (chr15)**
   - Download Ensembl chr15 FASTA (`.fa.gz` → `.fa`)
   - Create `.fai` (samtools faidx)
   - Create `.dict` (GATK CreateSequenceDictionary)
   - Build **BWA index** for canonical `reference/chr_15.fa`

2. **FASTQ download**
   - SRA Toolkit download + conversion to gzipped paired FASTQs
   - Auto-generate `samples.tsv` (for inspection; not required as an input)

3. **Alignment + BAM processing**
   - Align with `bwa mem` → SAM
   - Sort → `results/bam/<sample>.sorted.bam`
   - Mark duplicates (Picard) → `results/dedup/<sample>.dedup.bam`
   - Build BAM index (Picard) → `results/dedup/<sample>.dedup.bai`

4. **Step 4 — Call Variants**
   - GATK4 `HaplotypeCaller`
   - Output: `results/vcfs/<sample>.raw_variants.vcf` (first-pass variant calling)

---

## 🧾 accessions.txt (Comments Supported)

Lines starting with `#` are ignored, so you can keep optional accessions ready:

```text
SRR2584863
# SRR2584866
# SRR2584868
```
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
│   └── gatk.yml    
├── rules/  
│   ├── download_fastq.smk  
│   ├── reference.smk  
│   ├── align.smk  
│   ├── sort_bam.smk  
│   ├── metrics.smk  
│   ├── mark_duplicates.smk  
│   └── index_bam.smk  
│   └── haplotypecaller.smk  
├── reference/                # generated  
├── fastq/                    # generated  
├── results/                  # generated  
│   ├── bam/  
│   ├── metrics/  
│   └── dedup/  
│   └── vcfs/  
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
Edit config.yaml as needed. Common keys:
- accessions_file
- fastq_dir
- reference_base_url, references, reference_canonical, reference_dir
- bam_dir, metrics_dir, dedup_dir, vcf_dir
- samples_tsv (generated) 

---

## ▶️ Running the Workflow
### Dry run (validate DAG)  
snakemake -n -p --use-conda --cores 1 results/dedup/SRR2584863.dedup.bai

### Run Step 4: HaplotypeCaller for one sample 
snakemake -p --use-conda --cores 4 --rerun-incomplete --keep-going \  
&nbsp;&nbsp;&nbsp;&nbsp;results/vcfs/SRR2584863.raw_variants.vcf

### Run everything in rule all
snakemake -p --use-conda --cores 4 --rerun-incomplete --keep-going

---

## 👤 Contact
For issues or contributions, open a GitHub issue or pull request.
remco.crefcoeur@students.fhnw.ch
---