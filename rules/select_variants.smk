import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 5 — Split Variants (GATK4 SelectVariants)
# Input: raw_variants.vcf (first-pass call), reference
# Output: raw_snps.vcf and raw_indels.vcf
# Notes: This is used for downstream bootstrapping steps (e.g., BQSR).
# ------------------------------------------------------------------
VCF_DIR = config.get("vcf_dir", "results/vcfs")
Path(VCF_DIR).mkdir(parents=True, exist_ok=True)


def vcf_base(sample: str) -> str:
    return os.path.join(VCF_DIR, f"{sample}.raw_variants.vcf")


rule select_snps:
    input:
        vcf=lambda wc: vcf_base(wc.sample),
        idx=lambda wc: vcf_base(wc.sample) + ".idx",
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_snps.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_snps.vcf.idx"),
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


rule select_indels:
    input:
        vcf=lambda wc: vcf_base(wc.sample),
        idx=lambda wc: vcf_base(wc.sample) + ".idx",
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_indels.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_indels.vcf.idx"),
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