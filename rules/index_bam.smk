import os

PICARD_ENV = "../envs/picard.yml"  # reuse Picard conda env
DEDUP_DIR = config.get("dedup_dir", "results/dedup")

rule build_bam_index:
    """
    Build BAM index for deduplicated BAMs using Picard.
    """
    input:
        bam = os.path.join(DEDUP_DIR, "{sample}.dedup.bam")
    output:
        bai = os.path.join(DEDUP_DIR, "{sample}.dedup.bai")
    threads: 1
    conda: PICARD_ENV
    shell:
        r"""
        set -euo pipefail
        picard BuildBamIndex INPUT={input.bam} OUTPUT={output.bai}
        """