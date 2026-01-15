import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 12 — AnalyzeCovariates (BQSR recalibration report)
# Tool: GATK4 AnalyzeCovariates
# Input: recal_data.table (before), post_recal_data.table (after)
# Output: recalibration_plots.pdf + recalibration_plots.csv
# ------------------------------------------------------------------

BQSR_DIR = config["bqsr_dir"]
LOG_DIR = config.get("log_dir", "results/logs")

Path(BQSR_DIR).mkdir(parents=True, exist_ok=True)
Path(LOG_DIR).mkdir(parents=True, exist_ok=True)


rule analyze_covariates:
    input:
        before=os.path.join(BQSR_DIR, "{sample}.recal_data.table"),
        after=os.path.join(BQSR_DIR, "{sample}.post_recal_data.table"),
    output:
        pdf=os.path.join(BQSR_DIR, "{sample}.recalibration_plots.pdf"),
        csv=os.path.join(BQSR_DIR, "{sample}.recalibration_plots.csv"),
    log:
        os.path.join(LOG_DIR, "{sample}.analyze_covariates.log"),
    threads: 2
    conda:
        "../envs/gatk.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"

        gatk AnalyzeCovariates \
          -before {input.before} \
          -after {input.after} \
          -plots {output.pdf} \
          -csv {output.csv} \
          > {log} 2>&1
        """