import os
from pathlib import Path
import pandas as pd

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------

BAM_DIR = config.get("bam_dir", "results/bam")
Path(BAM_DIR).mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Samples
# ------------------------------------------------------------------

samples_df = pd.read_csv(config["samples_tsv"], sep="\t")

# Expect columns: sample, r1, r2
SAMPLE_PATHS = {
    row["sample"]: (row["r1"], row["r2"])
    for _, row in samples_df.iterrows()
}

SAMPLES = list(SAMPLE_PATHS.keys())

# ------------------------------------------------------------------
# Canonical reference
# ------------------------------------------------------------------

REF_DIR = config["reference_dir"]
CANON_FILE = config.get("reference_canonical")[0]
CANON_FA = os.path.join(REF_DIR, CANON_FILE)

# BWA index files required by bwa mem
BWA_INDEX_EXTS = ["amb", "ann", "bwt", "pac", "sa"]
BWA_INDEX_FILES = expand("{fa}.{ext}", fa=CANON_FA, ext=BWA_INDEX_EXTS)

# ------------------------------------------------------------------
# Alignment
# ------------------------------------------------------------------

rule bwa_align:
    """
    Align paired-end reads using BWA MEM.
    Produces one SAM file per sample.
    """

    input:
        r1 = lambda wc: SAMPLE_PATHS[wc.sample][0],
        r2 = lambda wc: SAMPLE_PATHS[wc.sample][1],
        ref = CANON_FA,
        index = BWA_INDEX_FILES

    output:
        sam = os.path.join(BAM_DIR, "{sample}.sam")

    threads: config["threads"]

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
