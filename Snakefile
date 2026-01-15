import os
import pandas as pd
from pathlib import Path

# Load config
configfile: "config.yaml"

# ------------------------------------------------------------------
# FASTQ / accessions (download stage)
# ------------------------------------------------------------------

def read_accessions(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip() and not l.startswith("#")]

ACCESSIONS = read_accessions(config["accessions_file"])

FASTQ_R1 = expand(
    "{fastq_dir}/{acc}_1.fastq.gz",
    fastq_dir=config["fastq_dir"],
    acc=ACCESSIONS,
)
FASTQ_R2 = expand(
    "{fastq_dir}/{acc}_2.fastq.gz",
    fastq_dir=config["fastq_dir"],
    acc=ACCESSIONS,
)

# ------------------------------------------------------------------
# Samples (alignment stage)
# ------------------------------------------------------------------

samples_df = pd.read_csv(config["samples_tsv"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

# ------------------------------------------------------------------
# Reference (raw + canonical)
# ------------------------------------------------------------------

REF_DIR = config["reference_dir"]
REF_FA = os.path.join(REF_DIR, config["references"][0].replace(".gz", ""))
CANON_FA = os.path.join(REF_DIR, config["reference_canonical"][0])
REF_FAIDX = CANON_FA + ".fai"
REF_DICT = CANON_FA.replace(".fa", ".dict")

# ------------------------------------------------------------------
# BAM outputs (from align.smk)
# ------------------------------------------------------------------

BAM_DIR = config.get("bam_dir", "results/bam")
BAMS = expand(
    os.path.join(BAM_DIR, "{sample}.sam"),
    sample=SAMPLES,
)

# ------------------------------------------------------------------
# Sorted BAMS (from sort_bam.smk)
# -----------------------------------------------------------------

SORTED_BAMS = expand(
    "{bam_dir}/{sample}.sorted.bam",
    bam_dir=config["bam_dir"],
    sample=SAMPLES
)

# ------------------------------------------------------------------
# Metrics (from metrics.smk)
# ------------------------------------------------------------------

METRICS = expand(
    "{metrics_dir}/{sample}.alignment_metrics.txt",
    metrics_dir=config["metrics_dir"],
    sample=SAMPLES
)

# ------------------------------------------------------------------
# Deduplicated BAMs (from mark_duplicates.smk)
# ------------------------------------------------------------------

DEDUP_DIR = config.get("dedup_dir", "results/dedup")

DEDUP_BAMS = expand(
    "{dedup_dir}/{sample}.dedup.bam",
    dedup_dir=DEDUP_DIR,
    sample=SAMPLES
)

DEDUP_METRICS = expand(
    "{dedup_dir}/{sample}.metrics.txt",
    dedup_dir=DEDUP_DIR,
    sample=SAMPLES
)

# NOTE: small fix approved — do NOT mark these as temp
DEDUP_BAI = expand(
    "results/dedup/{sample}.dedup.bai", sample=SAMPLES
)

# ------------------------------------------------------------------
# Realignment targets (Step 6)
# ----------------------------------------------------------------

REALIGN_DIR = config.get("realign_dir", "results/realign")

REALIGN_TARGETS = expand(
    "{realign_dir}/{sample}.realignment_targets.list",
    realign_dir=REALIGN_DIR,
    sample=SAMPLES
)

# ------------------------------------------------------------------
# Include modular rules
# ------------------------------------------------------------------

include: "rules/download_fastq.smk"
include: "rules/reference.smk"
include: "rules/align.smk"
include: "rules/sort_bam.smk"
include: "rules/metrics.smk"
include: "rules/mark_duplicates.smk"
include: "rules/index_bam.smk"
include: "rules/realigner_target_creator.smk"

# ------------------------------------------------------------------
# Rule all
# ------------------------------------------------------------------

rule all:
    input:
        config["samples_tsv"],
        FASTQ_R1,
        FASTQ_R2,
        REF_FA,
        CANON_FA,
        REF_FAIDX,
        REF_DICT,
        SORTED_BAMS,
        METRICS,
        DEDUP_BAMS,
        DEDUP_METRICS,
        DEDUP_BAI,
        REALIGN_TARGETS