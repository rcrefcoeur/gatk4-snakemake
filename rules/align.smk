import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 3 — Align Reads
# Tool: bwa mem
# Input: fastq/<sample>_1.fastq.gz, fastq/<sample>_2.fastq.gz, reference/chr_15.fa (+ BWA index)
# Output: results/bam/<sample>.sam
# ------------------------------------------------------------------
BAM_DIR = config.get("bam_dir", "results/bam")
Path(BAM_DIR).mkdir(parents=True, exist_ok=True)


rule align:
    """
    Align paired-end reads to canonical reference with bwa mem.
    """
    input:
        r1=os.path.join(config["fastq_dir"], "{sample}_1.fastq.gz"),
        r2=os.path.join(config["fastq_dir"], "{sample}_2.fastq.gz"),
        ref=CANON_FA,
        idx=BWA_INDEX,  # ensures bwa index exists before alignment
    output:
        sam=os.path.join(BAM_DIR, "{sample}.sam"),
    threads: config.get("threads", 6)
    conda:
        "../envs/bwa.yml"
    shell:
        r"""
        set -euo pipefail

        bwa mem -M \
          -t {threads} \
          -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tLB:{wildcards.sample}\tPL:{config[platform]}' \
          {input.ref} \
          {input.r1} {input.r2} \
          > {output.sam}
        """