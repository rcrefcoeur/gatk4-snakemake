import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 7 — Filter INDELs (GATK4 VariantFiltration)
# Input: raw_indels.vcf, reference
# Output: filtered_indels.vcf (+ filtered_indels.vcf.idx)
# Notes:
#   - Variants remain in the output VCF; failing variants are tagged in FILTER,
#     passing variants are tagged PASS. PASS-only extraction happens next.
#   - Filter thresholds follow Broad starter recommendations (as in NYU pipeline).
# ------------------------------------------------------------------
VCF_DIR = config.get("vcf_dir", "results/vcfs")
Path(VCF_DIR).mkdir(parents=True, exist_ok=True)


rule filter_indels:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_indels.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.raw_indels.vcf.idx"),
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.filtered_indels.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.filtered_indels.vcf.idx"),
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
          --filter-name "FS_filter" --filter "FS > 200.0" \
          --filter-name "SOR_filter" --filter "SOR > 10.0"
        """