import pandas as pd

configfile: "config.yaml"

# Include only rules needed for downloads + reference preparation
include: "rules/download_fastq.smk"
include: "rules/reference.smk"

# Load samples table
samples_df = pd.read_csv(config["samples_tsv"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

# Default target: FASTQ + chr15 reference
rule all:
    input:
        # FASTQ files
        expand("fastq/{sample}_1.fastq.gz", sample=SAMPLES),
        expand("fastq/{sample}_2.fastq.gz", sample=SAMPLES),

        # chr15-only reference (adjust to match your reference.smk outputs)
        "reference/chr15.fa",
        "reference/chr15.fa.fai",
        "reference/chr15.dict"
