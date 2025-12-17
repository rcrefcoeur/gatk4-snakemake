import os

# Load config
configfile: "config.yaml"

# Accessions
def read_accessions(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip() and not l.startswith("#")]

ACCESSIONS = read_accessions(config["accessions_file"])

include: "rules/download_fastq.smk"
include: "rules/reference.smk"

# FASTQ outputs
FASTQ_R1 = expand("{fastq_dir}/{acc}_1.fastq.gz", fastq_dir=config["fastq_dir"], acc=ACCESSIONS)
FASTQ_R2 = expand("{fastq_dir}/{acc}_2.fastq.gz", fastq_dir=config["fastq_dir"], acc=ACCESSIONS)

# Reference outputs (use canonical name from config)
REF_DIR = config["reference_dir"]
REF_FILE = config["references"][0].replace(".gz", "")
CANON_FILE = config.get("reference_canonical", [REF_FILE])[0]  # ensure list access

REF_FA = os.path.join(REF_DIR, REF_FILE)
CANON_FA = os.path.join(REF_DIR, CANON_FILE)
REF_FAIDX = [REF_FA + ".fai", CANON_FA + ".fai"]
REF_DICT = [REF_FA.replace(".fa", ".dict"), CANON_FA.replace(".fa", ".dict")]

# Rule all
rule all:
    input:
        config["samples_tsv"],
        FASTQ_R1,
        FASTQ_R2,
        REF_FA,
        CANON_FA,
        REF_FAIDX,
        REF_DICT,