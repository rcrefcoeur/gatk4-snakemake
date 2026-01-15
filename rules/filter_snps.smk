import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 6 — Filter SNPs (GATK4 VariantFiltration)
# Input: raw_snps.vcf, reference
# Output: filtered_snps.vcf (+ filtered_snps.vcf.idx)
# Notes:
#   - Variants remain in the output VCF; failing variants are tagged in FILTER,
#     passing variants are tagged PASS. PASS-only extraction happens later.
#   - Filter thresholds follow Broad starter recommendations (as in NYU pipeline).
# ------------------------------------------------------------------
VCF_DIR = config.get("vcf_dir", "results/vcfs")
Path(VCF_DIR).mkdir(parents=True, exist_ok=True)


rule filter_snps:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_snps.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_snps.vcf.idx"),
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.filtered_snps.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.filtered_snps.vcf.idx"),
    threads: 2
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk VariantFiltration \
          -R {input.ref} \
          -V {input.vcf} \
          -O {output.vcf} \
          --filter-name "QD_filter" --filter "QD < 2.0" \
          --filter-name "FS_filter" --filter "FS > 60.0" \
          --filter-name "MQ_filter" --filter "MQ < 40.0" \
          --filter-name "SOR_filter" --filter "SOR > 4.0" \
          --filter-name "MQRankSum_filter" --filter "MQRankSum < -12.5" \
          --filter-name "ReadPosRankSum_filter" --filter "ReadPosRankSum < -8.0"
        """