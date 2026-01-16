import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 16 — Filter Indels (FINAL; post-BQSR)
# Tool: GATK4 VariantFiltration
# Input: raw_indels_recal.vcf, reference
# Output: filtered_indels_final.vcf (+ .idx)
#
# Notes:
# - Variants remain in the output VCF; failing variants are tagged in FILTER,
#   passing variants are tagged PASS.
# - Filter thresholds follow Broad starter recommendations.
# ------------------------------------------------------------------
VCF_DIR = config.get("vcf_dir", "results/vcfs")
REF_DIR = config["reference_dir"]
CANON_FA = os.path.join(REF_DIR, config["reference_canonical"][0])

Path(VCF_DIR).mkdir(parents=True, exist_ok=True)


rule filter_indels_final:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.raw_indels_recal.vcf"),
        vcf_idx=os.path.join(VCF_DIR, "{sample}.raw_indels_recal.vcf.idx"),
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        vcf=os.path.join(VCF_DIR, "{sample}.filtered_indels_final.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.filtered_indels_final.vcf.idx"),
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