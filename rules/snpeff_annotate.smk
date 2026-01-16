import os
from pathlib import Path

# ------------------------------------------------------------------
# Step 17 — Annotate SNPs and Predict Effects (SnpEff)
# ------------------------------------------------------------------

VCF_DIR = config["vcf_dir"]
SNPEFF_DB = config["snpeff_db"]
SNPEFF_DATA_DIR = config.get("snpeff_data_dir", "resources/snpeff")

Path(VCF_DIR).mkdir(parents=True, exist_ok=True)
Path(SNPEFF_DATA_DIR).mkdir(parents=True, exist_ok=True)

rule snpeff_annotate_snps_final:
    input:
        vcf=os.path.join(VCF_DIR, "{sample}.filtered_snps_final.vcf"),
        idx=os.path.join(VCF_DIR, "{sample}.filtered_snps_final.vcf.idx"),
    output:
        ann_vcf=os.path.join(VCF_DIR, "{sample}.filtered_snps_final.ann.vcf"),
        html=os.path.join(VCF_DIR, "{sample}.snpeff_summary.html"),
        genes=os.path.join(VCF_DIR, "{sample}.snpEff_genes.txt"),
    conda:
        "../envs/snpeff.yml"
    shell:
        r"""
        set -euo pipefail

        # Give Java more heap; GRCh38 databases can require several GB
        snpEff -Xmx8g -dataDir {SNPEFF_DATA_DIR} -v \
          -stats {output.html} \
          {SNPEFF_DB} \
          {input.vcf} \
          > {output.ann_vcf}

        # Extract unique gene names from ANN field (Allele|Annotation|Impact|Gene_Name|...)
        awk -F'\t' '
          BEGIN {{ }}
          /^#/ {{ next }}
          {{
            info=$8
            if (match(info, /ANN=([^;]+)/, a)) {{
              n=split(a[1], ann, ",")
              for (i=1; i<=n; i++) {{
                m=split(ann[i], f, "\\|")
                if (m>=4 && f[4] != "") print f[4]
              }}
            }}
          }}
        ' {output.ann_vcf} | sort -u > {output.genes}
        """