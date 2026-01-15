import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 3b — Sort Alignments
# Tool: samtools sort
# Input: results/bam/<sample>.sam
# Output: results/bam/<sample>.sorted.bam
# ------------------------------------------------------------------
BAM_DIR = config["bam_dir"]
Path(BAM_DIR).mkdir(parents=True, exist_ok=True)

PICARD_ENV = "../envs/picard.yml"

rule sort_sam:
    """
    Sort SAM by coordinate and convert to BAM using Picard SortSam.
    """
    input:
        sam = os.path.join(BAM_DIR, "{sample}.sam")
    output:
        bam = os.path.join(BAM_DIR, "{sample}.sorted.bam")
    threads: 2
    conda: f"{PICARD_ENV}"
    shell:
        r"""
        set -euo pipefail

        picard SortSam \
            INPUT={input.sam} \
            OUTPUT={output.bam} \
            SORT_ORDER=coordinate \
            CREATE_INDEX=false \
            VALIDATION_STRINGENCY=STRICT
        """
