import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 4 — Call Variants (GATK4 HaplotypeCaller)
# Input: dedup + indexed BAM, reference
# Output: raw_variants.vcf (first-pass call)
# ------------------------------------------------------------------
VCF_DIR = config.get("vcf_dir", "results/vcfs")
Path(VCF_DIR).mkdir(parents=True, exist_ok=True)


rule haplotypecaller_raw:
    input:
        bam=os.path.join(config.get("dedup_dir", "results/dedup"), "{sample}.dedup.bam"),
        bai=os.path.join(config.get("dedup_dir", "results/dedup"), "{sample}.dedup.bai"),
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_variants.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_variants.vcf.idx"),
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