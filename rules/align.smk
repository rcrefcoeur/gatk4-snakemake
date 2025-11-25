# rules/align.smk
rule bwa_mem:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample][0],
        r2 = lambda wildcards: config["samples"][wildcards.sample][1],
        ref = config["reference"]
    output:
        temp("work/{sample}.aligned.bam")
    threads: 6
    conda:
        "envs/bwa.yml"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -b - > {output}
        """
