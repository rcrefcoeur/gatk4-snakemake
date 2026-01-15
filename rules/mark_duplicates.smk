import os
from pathlib import Path

DEDUP_DIR = config.get("dedup_dir", "results/dedup")
Path(DEDUP_DIR).mkdir(parents=True, exist_ok=True)

PICARD_ENV = "../envs/picard.yml"


rule mark_duplicates:
    """
    Mark duplicates in sorted BAM using Picard MarkDuplicates.
    """
    input:
        bam=os.path.join(BAM_DIR, "{sample}.sorted.bam"),
    output:
        bam=os.path.join(DEDUP_DIR, "{sample}.dedup.bam"),
        metrics=os.path.join(DEDUP_DIR, "{sample}.metrics.txt"),
    threads: 2
    conda: PICARD_ENV
    shell:
        r"""
        set -euo pipefail

        # Avoid double indexing; indexing is done by rules/index_bam.smk
        picard MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            CREATE_INDEX=false \
            VALIDATION_STRINGENCY=STRICT
        """
