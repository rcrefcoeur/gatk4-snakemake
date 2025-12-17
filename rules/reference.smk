import os
from pathlib import Path

REF_DIR = config["reference_dir"]

# pick the first active reference from the config
REF_FILE = config["references"][0].replace(".gz", "")
REF_URL = config["reference_base_url"] + config["references"][0]

# pick canonical name from config if defined, else default to REF_FILE
CANON_FILE = config.get("reference_canonical", [REF_FILE])[0]

REF_FA = os.path.join(REF_DIR, REF_FILE)
CANON_FA = os.path.join(REF_DIR, CANON_FILE)
CANON_FAI = CANON_FA + ".fai"
CANON_DICT = CANON_FA.replace(".fa", ".dict")

def ensure_ref_dir():
    Path(REF_DIR).mkdir(parents=True, exist_ok=True)

CONDA_ENV = "../envs/reference.yml"

rule download_reference:
    """Download reference fasta and gunzip."""
    output:
        fa = REF_FA
    threads: 1
    conda: CONDA_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {REF_DIR}
        wget -c -O {output.fa}.gz {REF_URL}
        gunzip -f {output.fa}.gz
        """

rule faidx:
    """Index reference with samtools and create canonical names."""
    input:
        fa = REF_FA
    output:
        fai = REF_FA + ".fai",
        canon_fa = CANON_FA,
        canon_fai = CANON_FAI
    threads: 6
    conda: CONDA_ENV
    run:
        ensure_ref_dir()
        shell("samtools faidx {input.fa}")

        # canonical fasta
        src_fa = os.path.abspath(input.fa)
        dst_fa = os.path.abspath(output.canon_fa)
        try:
            if os.path.exists(dst_fa):
                os.remove(dst_fa)
            os.symlink(src_fa, dst_fa)
        except Exception:
            import shutil
            shutil.copy2(src_fa, dst_fa)

        # canonical fai
        src_fai = src_fa + ".fai"
        dst_fai = os.path.abspath(output.canon_fai)
        try:
            if os.path.exists(dst_fai):
                os.remove(dst_fai)
            os.symlink(src_fai, dst_fai)
        except Exception:
            import shutil
            shutil.copy2(src_fai, dst_fai)

rule bwa_index:
    """Index reference with BWA."""
    input:
        fa = REF_FA
    output:
        bwt = REF_FA + ".bwt"  # representative BWA output
    threads: 6
    conda: CONDA_ENV
    shell:
        "bwa index {input.fa}"

rule dict:
    """Create GATK sequence dictionary and canonical dict."""
    input:
        fa = REF_FA
    output:
        dict = REF_FA.replace(".fa", ".dict"),
        canon_dict = CANON_DICT
    threads: 1
    conda: CONDA_ENV
    run:
        shell("gatk CreateSequenceDictionary -R {input.fa} -O {output.dict}")

        # canonical dict
        src_dict = os.path.abspath(output.dict)
        dst_dict = os.path.abspath(output.canon_dict)
        try:
            if os.path.exists(dst_dict):
                os.remove(dst_dict)
            os.symlink(src_dict, dst_dict)
        except Exception:
            import shutil
            shutil.copy2(src_dict, dst_dict)
