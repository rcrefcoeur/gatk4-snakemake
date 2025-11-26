rule align:
    input:
        fa=config["reference_fa"],
        fai=config["reference_fa"] + ".fai",
        dict=config["reference_fa"].replace(".fa", ".dict"),
        bwt=config["reference_fa"] + ".bwt",
        r1=lambda wildcards: SAMPLE_PATHS[wildcards.sample][0],
        r2=lambda wildcards: SAMPLE_PATHS[wildcards.sample][1]  # can be None
    output:
        bam=lambda wildcards: f"{config['outdir']}/{wildcards.sample}.sorted.bam",
        bai=lambda wildcards: f"{config['outdir']}/{wildcards.sample}.sorted.bam.bai"
    threads: 6
    conda: "../envs/reference.yml"
    shell:
        r"""
        mkdir -p {config['outdir']}
        if [ -z "{input.r2}" ]; then
            bwa mem -t {threads} {input.fa} {input.r1} | \
            samtools sort -@ {threads} -o {output.bam}
        else
            bwa mem -t {threads} {input.fa} {input.r1} {input.r2} | \
            samtools sort -@ {threads} -o {output.bam}
        fi
        samtools index {output.bam}
        """
