# Snakefile

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

REF_FILE = config["references"][0].replace(".gz", "")
CANON_FILE = config.get("reference_canonical", [REF_FILE])[0]

REF_FA = os.path.join(REF_DIR, REF_FILE)
CANON_FA = os.path.join(REF_DIR, CANON_FILE)

REF_FAIDX = [
    REF_FA + ".fai",
    CANON_FA + ".fai",
]

REF_DICT = [
    REF_FA.replace(".fa", ".dict"),
    CANON_FA.replace(".fa", ".dict"),
]

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
# Include rules
# ------------------------------------------------------------------

include: "rules/download_fastq.smk"
include: "rules/reference.smk"
include: "rules/align.smk"
include: "rules/sort_bam.smk"


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
        SORTED_BAMS
