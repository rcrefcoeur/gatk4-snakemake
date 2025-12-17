# rules/download_fastq.smk

rule download_fastq:
    input:
        accessions = config["accessions_file"]
    output:
        r1 = expand("{fastq_dir}/{acc}_1.fastq.gz", fastq_dir=config["fastq_dir"], acc=[l.strip() for l in open(config["accessions_file"]) if l.strip()]),
        r2 = expand("{fastq_dir}/{acc}_2.fastq.gz", fastq_dir=config["fastq_dir"], acc=[l.strip() for l in open(config["accessions_file"]) if l.strip()]),
        samples = config["samples_tsv"]
    threads: 4
    conda: "../envs/sra-tools.yml"
    log:
        "../logs/download_fastq.log"
    shell:
        r"""
        bash scripts/download_fastq.sh \
            {input.accessions} \
            {config[fastq_dir]} \
            {output.samples} \
            > {log} 2>&1
        """
