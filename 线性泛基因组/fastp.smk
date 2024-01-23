SAMPLES, = glob_wildcards("01.raw/{sample}_1.fq.gz")
print("Total samples: ",len(SAMPLES))

rule all:
    input:
        expand("01.clean.fq/{sample}_1.fq.gz",sample=SAMPLES)

rule run_fastp:
    input:
        r1 = "01.raw/{sample}_1.fq.gz",
        r2 = "01.raw/{sample}_2.fq.gz"
    output:
        r1 = "01.clean.fq/{sample}_1.fq.gz",
        r2 = "01.clean.fq/{sample}_2.fq.gz",
        html = "01.fastp.report/{sample}.html",
        json = "01.fastp.report/{sample}.json"
    threads:
        8
    resources:
        mem = 16000
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -w {threads} -h {output.html} -j {output.json}
        """