import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 14 — Extract SNPs & Indels (post-BQSR)
# Tool: GATK4 SelectVariants
# Input: raw_variants_recal.vcf (from Step 13), reference
# Output: raw_snps_recal.vcf, raw_indels_recal.vcf (+ .idx)
#
# Why:
# - Split SNPs and INDELs so they can be filtered/handled independently.
# ------------------------------------------------------------------

VCF_DIR = config.get("vcf_dir", "results/vcfs")
Path(VCF_DIR).mkdir(parents=True, exist_ok=True)

rule select_snps_recal:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_variants_recal.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_variants_recal.vcf.idx"),
        ref="reference/chr_15.fa",
        fai="reference/chr_15.fa.fai",
        dict="reference/chr_15.dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_snps_recal.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_snps_recal.vcf.idx"),
    threads: 2
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk SelectVariants \
          -R {input.ref} \
          -V {input.vcf} \
          --select-type-to-include SNP \
          -O {output.vcf}
        """

rule select_indels_recal:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_variants_recal.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_variants_recal.vcf.idx"),
        ref="reference/chr_15.fa",
        fai="reference/chr_15.fa.fai",
        dict="reference/chr_15.dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_indels_recal.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_indels_recal.vcf.idx"),
    threads: 2
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk SelectVariants \
          -R {input.ref} \
          -V {input.vcf} \
          --select-type-to-include INDEL \
          -O {output.vcf}
        """