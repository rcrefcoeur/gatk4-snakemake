from pathlib import Path

# ------------------------------------------------------------------
# Step 1 — Download FASTQ + Generate samples.tsv
# Input: accessions.txt
# Output: fastq/<ACC>_1.fastq.gz, fastq/<ACC>_2.fastq.gz, samples.tsv
# Notes: comments (#) and empty lines in accessions.txt are ignored by the script
# ------------------------------------------------------------------
LOG_DIR = "results/logs"
Path(LOG_DIR).mkdir(parents=True, exist_ok=True)


rule download_fastq:
    input:
        accessions=config["accessions_file"]
    output:
        r1=expand("{fastq_dir}/{sample}_1.fastq.gz", fastq_dir=config["fastq_dir"], sample=SAMPLES),
        r2=expand("{fastq_dir}/{sample}_2.fastq.gz", fastq_dir=config["fastq_dir"], sample=SAMPLES),
        samples=config["samples_tsv"]
    threads: 4
    conda: "../envs/sra-tools.yml"
    log:
        f"{LOG_DIR}/download_fastq.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {LOG_DIR}
        bash scripts/download_fastq.sh \
            {input.accessions} \
            {config[fastq_dir]} \
            {output.samples} \
            > {log} 2>&1
        """