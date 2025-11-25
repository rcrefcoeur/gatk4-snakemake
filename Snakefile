# Snakefile (top-level)
import os
configfile: "config.yaml"

# default target
rule all:
    input:
        expand("{outdir}/{sample}.vcf.gz", outdir=config["outdir"], sample=list(config["samples"].keys()))

# include modular rules
include: "rules/align.smk"
include: "rules/markdup.smk"
include: "rules/bqsr.smk"
include: "rules/call_variants.smk"
include: "rules/snpeff.smk"