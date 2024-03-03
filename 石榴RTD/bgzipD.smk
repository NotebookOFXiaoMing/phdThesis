SRR, = glob_wildcards("clean.fastq/{srr}_1.fq.gz")

print("Total sample size: ",len(SRR))

rule all:
    input:
        expand("clean.fastq/{srr}_1.fq",srr=SRR)

rule bgzip:
    input:
        r1 = "clean.fastq/{srr}_1.fq.gz",
        r2 = "clean.fastq/{srr}_2.fq.gz"
    output:
        r1 = "clean.fastq/{srr}_1.fq",
        r2 = "clean.fastq/{srr}_2.fq"
    threads:
        8
    resources:
        mem_mb = 24000
    shell:
        """
        bgzip -d {input.r1}
        bgzip -d {input.r2}
        """