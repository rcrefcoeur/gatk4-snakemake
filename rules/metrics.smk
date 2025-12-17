import os
from pathlib import Path

METRICS_DIR = config.get("metrics_dir", "results/metrics")
Path(METRICS_DIR).mkdir(parents=True, exist_ok=True)

PICARD_ENV = "../envs/picard.yml"

rule collect_metrics:
    """
    Collect alignment and insert size metrics using Picard and compute depth with samtools
    """
    input:
        bam = os.path.join(config["bam_dir"], "{sample}.sorted.bam"),
        ref = os.path.join(config["reference_dir"], config.get("reference_canonical", ["chr_15"])[0])
    output:
        align_metrics = os.path.join(METRICS_DIR, "{sample}.alignment_metrics.txt"),
        insert_metrics = os.path.join(METRICS_DIR, "{sample}.insert_metrics.txt"),
        insert_hist = os.path.join(METRICS_DIR, "{sample}.insert_size_histogram.pdf"),
        depth = os.path.join(METRICS_DIR, "{sample}.depth_out.txt")
    threads: 2
    conda: PICARD_ENV
    shell:
        r"""
        picard CollectAlignmentSummaryMetrics \
            R={input.ref} \
            I={input.bam} \
            O={output.align_metrics}

        picard CollectInsertSizeMetrics \
            INPUT={input.bam} \
            OUTPUT={output.insert_metrics} \
            HISTOGRAM_FILE={output.insert_hist}

        samtools depth -a {input.bam} > {output.depth}
        """
