import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 9 — Base Quality Score Recalibration (BQSR) #1
# Tool: GATK4 BaseRecalibrator
# Input: dedup + indexed BAM, bqsr_snps.vcf, bqsr_indels.vcf, reference
# Output: recal_data.table
#
# Step 10 — Apply BQSR
# Tool: GATK4 ApplyBQSR
# Input: recal_data.table, dedup + indexed BAM, reference
# Output: recal_reads.bam (analysis-ready BAM)
# ------------------------------------------------------------------
BQSR_DIR = config.get("bqsr_dir", "results/bqsr")
Path(BQSR_DIR).mkdir(parents=True, exist_ok=True)

DEDUP_DIR = config.get("dedup_dir", "results/dedup")
VCF_DIR = config.get("vcf_dir", "results/vcfs")


rule base_recalibrator:
    input:
        bam=os.path.join(DEDUP_DIR, "{sample}.dedup.bam"),
        bai=os.path.join(DEDUP_DIR, "{sample}.dedup.bai"),
        snps=os.path.join(VCF_DIR, "{sample}.bqsr_snps.vcf"),
        snps_idx=os.path.join(VCF_DIR, "{sample}.bqsr_snps.vcf.idx"),
        indels=os.path.join(VCF_DIR, "{sample}.bqsr_indels.vcf"),
        indels_idx=os.path.join(VCF_DIR, "{sample}.bqsr_indels.vcf.idx"),
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        table=os.path.join(BQSR_DIR, "{sample}.recal_data.table"),
    threads: 4
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk BaseRecalibrator \
          -R {input.ref} \
          -I {input.bam} \
          --known-sites {input.snps} \
          --known-sites {input.indels} \
          -O {output.table}
        """


rule apply_bqsr:
    input:
        bam=os.path.join(DEDUP_DIR, "{sample}.dedup.bam"),
        bai=os.path.join(DEDUP_DIR, "{sample}.dedup.bai"),
        table=os.path.join(BQSR_DIR, "{sample}.recal_data.table"),
        ref=CANON_FA,
        ref_fai=CANON_FA + ".fai",
        ref_dict=os.path.splitext(CANON_FA)[0] + ".dict",
    output:
        bam=os.path.join(BQSR_DIR, "{sample}.recal_reads.bam"),
    threads: 4
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail

        gatk ApplyBQSR \
          -R {input.ref} \
          -I {input.bam} \
          --bqsr-recal-file {input.table} \
          -O {output.bam}
        """