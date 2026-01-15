import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 8 — Exclude Filtered Variants (GATK4 SelectVariants)
# Input: filtered_snps.vcf, filtered_indels.vcf
# Output: bqsr_snps.vcf, bqsr_indels.vcf (PASS-only sets for BQSR known-sites)
# Notes:
#   - `--exclude-filtered` drops variants that failed hard filters (keeps PASS).
# ------------------------------------------------------------------
VCF_DIR = config.get("vcf_dir", "results/vcfs")
Path(VCF_DIR).mkdir(parents=True, exist_ok=True)


rule bqsr_snps:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.filtered_snps.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.filtered_snps.vcf.idx"),
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.bqsr_snps.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.bqsr_snps.vcf.idx"),
    threads: 2
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk SelectVariants \
          --exclude-filtered \
          -V {input.vcf} \
          -O {output.vcf}
        """


rule bqsr_indels:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.filtered_indels.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.filtered_indels.vcf.idx"),
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.bqsr_indels.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.bqsr_indels.vcf.idx"),
    threads: 2
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk SelectVariants \
          --exclude-filtered \
          -V {input.vcf} \
          -O {output.vcf}
        """