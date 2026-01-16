import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 18 — Compile Statistics (parse_metrics.sh)
# Tool: bin/parse_metrics.sh
# Input: sample_id_* files from metrics/dedup/vcfs outputs
# Output:
#   - per-sample: results/reports/{sample}.report.csv
#   - combined:   results/reports/report.csv
#
# Notes:
# - On WSL/Windows clones, executable bits can be lost. We therefore run the
#   script explicitly via `bash` and only require that the file exists.
# - `insert` is reserved in Snakemake, so we use `insert_metrics`.
# ------------------------------------------------------------------

REPORT_DIR = config.get("report_dir", "results/reports")
PARSE_METRICS = config.get("parse_metrics_script", "bin/parse_metrics.sh")

METRICS_DIR = config["metrics_dir"]
DEDUP_DIR = config["dedup_dir"]
VCF_DIR = config["vcf_dir"]

Path(REPORT_DIR).mkdir(parents=True, exist_ok=True)

rule compile_sample_report:
    input:
        align_metrics=os.path.join(METRICS_DIR, "{sample}.alignment_metrics.txt"),
        insert_metrics=os.path.join(METRICS_DIR, "{sample}.insert_metrics.txt"),
        dedup_metrics=os.path.join(DEDUP_DIR, "{sample}.metrics.txt"),
        raw_snps=os.path.join(VCF_DIR, "{sample}.raw_snps.vcf"),
        filtered_snps=os.path.join(VCF_DIR, "{sample}.filtered_snps.vcf"),
        raw_snps_recal=os.path.join(VCF_DIR, "{sample}.raw_snps_recal.vcf"),
        filtered_snps_final=os.path.join(VCF_DIR, "{sample}.filtered_snps_final.vcf"),
        depth_out=os.path.join(METRICS_DIR, "{sample}.depth_out.txt"),
    output:
        csv=os.path.join(REPORT_DIR, "{sample}.report.csv"),
    params:
        script=PARSE_METRICS,
        metrics_dir=METRICS_DIR,
        dedup_dir=DEDUP_DIR,
        vcf_dir=VCF_DIR,
        report_dir=REPORT_DIR,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.report_dir}"

        if [ ! -f "{params.script}" ]; then
          echo "ERROR: parse_metrics script not found: {params.script}" >&2
          exit 2
        fi

        bash "{params.script}" "{wildcards.sample}" "{params.metrics_dir}" "{params.dedup_dir}" "{params.vcf_dir}" > "{output.csv}"
        """

rule compile_reports:
    input:
        lambda wildcards: expand(os.path.join(REPORT_DIR, "{sample}.report.csv"), sample=SAMPLES),
    output:
        os.path.join(REPORT_DIR, "report.csv"),
    params:
        reports=lambda wildcards, input: " ".join(input),
    shell:
        r"""
        set -euo pipefail

        head -n 1 "{input[0]}" > "{output}"

        for f in {params.reports}; do
          tail -n +2 "$f" >> "{output}"
        done
        """