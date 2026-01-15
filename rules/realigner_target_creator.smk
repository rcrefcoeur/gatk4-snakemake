import os
from pathlib import Path

GATK3_ENV = "../envs/gatk3.yml"
DEDUP_DIR = config.get("dedup_dir", "results/dedup")
REALIGN_DIR = config.get("realign_dir", "results/realign")

Path(REALIGN_DIR).mkdir(parents=True, exist_ok=True)


rule realigner_target_creator:
    """
    Create realignment targets using GATK 3.x (RealignerTargetCreator).

    We do NOT rely on bioconda's gatk-register (it can fail depending on tarball layout).
    Instead, we extract the user-provided tarball and cache GenomeAnalysisTK.jar locally.
    """
    input:
        bam=os.path.join(DEDUP_DIR, "{sample}.dedup.bam"),
        bai=os.path.join(DEDUP_DIR, "{sample}.dedup.bai"),
        ref=os.path.join(config["reference_dir"], config["reference_canonical"][0]),
    output:
        targets=os.path.join(REALIGN_DIR, "{sample}.realignment_targets.list"),
    threads: 2
    conda: GATK3_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {REALIGN_DIR}

        # 1) If the jar is already present in the conda env, use it.
        JAR="$(find "$CONDA_PREFIX" -name GenomeAnalysisTK.jar -print -quit 2>/dev/null || true)"

        # 2) Otherwise, extract from Broad tarball and cache it under .snakemake/
        if [ -z "$JAR" ]; then
          TARBALL="{config[gatk3_tarball]}"
          if [ ! -f "$TARBALL" ]; then
            echo "ERROR: GenomeAnalysisTK.jar not found in env and tarball not found:" >&2
            echo "  $TARBALL" >&2
            exit 2
          fi

          CACHE_DIR=".snakemake/gatk3-cache"
          JAR_CACHE="$CACHE_DIR/GenomeAnalysisTK.jar"

          if [ ! -f "$JAR_CACHE" ]; then
            mkdir -p "$CACHE_DIR"
            TMPDIR="$(mktemp -d)"
            # tarball may unpack into a subdirectory; we locate the jar after extraction
            tar -xjf "$TARBALL" -C "$TMPDIR"

            FOUND="$(find "$TMPDIR" -name GenomeAnalysisTK.jar -print -quit 2>/dev/null || true)"
            if [ -z "$FOUND" ]; then
              echo "ERROR: Could not find GenomeAnalysisTK.jar inside tarball:" >&2
              echo "  $TARBALL" >&2
              rm -rf "$TMPDIR"
              exit 2
            fi

            cp "$FOUND" "$JAR_CACHE"
            rm -rf "$TMPDIR"
          fi

          JAR="$JAR_CACHE"
        fi

        # Run RealignerTargetCreator
        java -Xmx4g -jar "$JAR" \
          -T RealignerTargetCreator \
          -R {input.ref} \
          -I {input.bam} \
          -o {output.targets}
        """