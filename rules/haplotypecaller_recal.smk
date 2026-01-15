import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 13 — Call Variants (GATK4 HaplotypeCaller, post-BQSR)
# Input: recal_reads.bam (+ index), reference
# Output: raw_variants_recal.vcf (second-pass callset)
# ------------------------------------------------------------------

BQSR_DIR = config["bqsr_dir"]
VCF_DIR = config.get("vcf_dir", "results/vcfs")

Path(BQSR_DIR).mkdir(parents=True, exist_ok=True)
Path(VCF_DIR).mkdir(parents=True, exist_ok=True)


rule build_recal_bam_index:
    input:
        bam=os.path.join(BQSR_DIR, "{sample}.recal_reads.bam"),
    output:
        bai=os.path.join(BQSR_DIR, "{sample}.recal_reads.bai"),
    threads: 2
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail
        gatk BuildBamIndex \
          -I {input.bam} \
          -O {output.bai}
        """


rule haplotypecaller_recal:
    input:
        bam=os.path.join(BQSR_DIR, "{sample}.recal_reads.bam"),
        bai=os.path.join(BQSR_DIR, "{sample}.recal_reads.bai"),
        ref=CANON_FA,
        fai=CANON_FA + ".fai",
        dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_variants_recal.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_variants_recal.vcf.idx"),
    threads: 4
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk HaplotypeCaller \
          -R {input.ref} \
          -I {input.bam} \
          -O {output.vcf}
        """