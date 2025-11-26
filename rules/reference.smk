# rules/reference.smk
import os
from pathlib import Path

REF_DIR = config["reference_dir"]
REF_FA = config["reference_fa"]         # e.g. reference/Homo_sapiens.GRCh38.dna.chromosome.15.fa
REF_URL = config["reference_url"]

# canonical names expected by rule all
CANON_FA = os.path.join(REF_DIR, "chr15.fa")
CANON_FAI = CANON_FA + ".fai"
CANON_DICT = CANON_FA.replace(".fa", ".dict")

# Helper: ensure reference dir exists
def ensure_ref_dir():
    Path(REF_DIR).mkdir(parents=True, exist_ok=True)

CONDA_ENV = "../envs/reference.yml"

rule download_chr15:
    output:
        fa = REF_FA
    threads: 1
    conda: CONDA_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {REF_DIR}
        # download to REF_FA.gz (if remote) then gunzip -> REF_FA
        wget -c -O {output.fa}.gz {REF_URL}
        gunzip -f {output.fa}.gz
        """

# create samtools faidx for REF_FA and also create canonical names (fa + fai)
rule faidx:
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

        # Create canonical fasta (symlink if possible, else copy)
        src_fa = os.path.abspath(str(input.fa))
        dest_fa = os.path.abspath(str(output.canon_fa))
        try:
            if os.path.islink(dest_fa) or os.path.exists(dest_fa):
                os.remove(dest_fa)
            os.symlink(src_fa, dest_fa)
        except Exception:
            import shutil
            shutil.copy2(src_fa, dest_fa)

        # Create canonical fai (symlink to the produced .fai)
        src_fai = os.path.abspath(str(input.fa)) + ".fai"
        dest_fai = os.path.abspath(str(output.canon_fai))
        try:
            if os.path.islink(dest_fai) or os.path.exists(dest_fai):
                os.remove(dest_fai)
            os.symlink(src_fai, dest_fai)
        except Exception:
            import shutil
            shutil.copy2(src_fai, dest_fai)

# create bwa index files (creates many files; use one representative output)
rule bwa_index:
    input:
        fa = REF_FA
    output:
        bwt = REF_FA + ".bwt"
    threads: 6
    conda: CONDA_ENV
    shell:
        "bwa index {input.fa}"

# create sequence dictionary for REF_FA and also for canonical name if needed
rule dict:
    input:
        fa = REF_FA
    output:
        dict = REF_FA.replace(".fa", ".dict"),
        canon_dict = CANON_DICT
    threads: 1
    conda: "../envs/reference.yml"
    run:
        shell("gatk CreateSequenceDictionary -R {input.fa} -O {output.dict}")
        # ensure a canonical dict for chr15.fa name too
        if os.path.abspath(str(input.fa)) != os.path.abspath(CANON_FA):
            try:
                if os.path.exists(output.canon_dict):
                    os.remove(output.canon_dict)
                os.symlink(os.path.abspath(output.dict), output.canon_dict)
            except Exception:
                import shutil
                shutil.copy2(os.path.abspath(output.dict), output.canon_dict)
