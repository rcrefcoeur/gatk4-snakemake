# =========================
# FILE: rules/reference.smk
# =========================
import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 2 — Reference Preparation
# Input: Ensembl reference FASTA (.fa.gz URL), specified in config.yaml
# Output:
#   - reference/<downloaded>.fa
#   - reference/<downloaded>.fa.fai
#   - reference/<canonical>.fa (symlink/copy)
#   - reference/<canonical>.fa.fai
#   - reference/<*.dict> dictionaries
#   - BWA index for canonical FASTA (.amb/.ann/.bwt/.pac/.sa)
# -----------------------------------------------------------------
REF_DIR = config["reference_dir"]
REF_GZ = config["references"][0]
REF_FILE = REF_GZ.replace(".gz", "")
REF_URL = config["reference_base_url"] + REF_GZ

CANON_FILE = config.get("reference_canonical", [REF_FILE])[0]

REF_FA = os.path.join(REF_DIR, REF_FILE)
CANON_FA = os.path.join(REF_DIR, CANON_FILE)

REF_FAI = REF_FA + ".fai"
CANON_FAI = CANON_FA + ".fai"

REF_DICT = os.path.splitext(REF_FA)[0] + ".dict"
CANON_DICT = os.path.splitext(CANON_FA)[0] + ".dict"

BWA_INDEX = [CANON_FA + ext for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]

CONDA_ENV = "../envs/reference.yml"


rule download_reference:
    """
    Download and unzip the reference FASTA from Ensembl.
    """
    output:
        fa=REF_FA
    conda:
        CONDA_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {REF_DIR}
        wget -c -O {output.fa}.gz {REF_URL}
        gunzip -f {output.fa}.gz
        """


rule faidx:
    """
    Create FASTA index (.fai) and canonicalize reference filename.
    """
    input:
        fa=REF_FA
    output:
        ref_fai=REF_FAI,
        canon_fa=CANON_FA,
        canon_fai=CANON_FAI,
    threads: 4
    conda:
        CONDA_ENV
    run:
        # Index the downloaded reference
        shell("samtools faidx {input.fa}")

        # Create canonical FASTA as a symlink (fallback to copy on Windows filesystems)
        src = os.path.abspath(input.fa)
        dst = os.path.abspath(output.canon_fa)
        if os.path.lexists(dst):
            os.remove(dst)
        try:
            os.symlink(src, dst)
        except OSError:
            import shutil
            shutil.copy2(src, dst)

        # Index the canonical FASTA too (separate .fai)
        shell("samtools faidx {output.canon_fa}")


rule create_dict:
    """
    Create sequence dictionary (.dict) for both reference and canonical FASTA.
    """
    input:
        fa=REF_FA,
        canon_fa=CANON_FA
    output:
        ref_dict=REF_DICT,
        canon_dict=CANON_DICT
    conda:
        CONDA_ENV
    shell:
        r"""
        set -euo pipefail
        gatk CreateSequenceDictionary -R {input.fa} -O {output.ref_dict}
        gatk CreateSequenceDictionary -R {input.canon_fa} -O {output.canon_dict}
        """


rule bwa_index:
    """
    Build BWA index for the canonical FASTA used by alignment.
    """
    input:
        fa=CANON_FA
    output:
        idx=BWA_INDEX
    threads: 2
    conda:
        CONDA_ENV
    shell:
        r"""
        set -euo pipefail
        bwa index {input.fa}
        """