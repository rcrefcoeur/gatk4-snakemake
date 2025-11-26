# rules/download_fastq.smk

ACCESSIONS_FILE = config.get("accessions_file", "accessions.txt")
OUTDIR = config.get("fastq_dir", "fastq")
OUT_SAMPLES = config.get("samples_tsv", "samples.tsv")

rule download_fastq:
    input:
        accessions = ACCESSIONS_FILE
    output:
        samples = OUT_SAMPLES
    threads: 4
    conda: "../envs/reference.yml"  # or a minimal env with curl, sra-tools, gzip
    shell:
        """
        bash scripts/download_fastq.sh {input.accessions} {OUTDIR} {output.samples}
        """
