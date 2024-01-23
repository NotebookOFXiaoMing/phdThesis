SAMPLES, = glob_wildcards("../01.clean.fq/{sample}_1.fq.gz")

print("Total sample: ",len(SAMPLES))

rule all:
    input:
        expand("{sample}.log",sample=SAMPLES)

rule run_megahit:
    input:
        r1 = "../01.clean.fq/{sample}_1.fq.gz",
        r2 = "../01.clean.fq/{sample}_2.fq.gz"
    output:
        "{sample}.log"
    threads:
        12
    resources:
        mem = 24000
    params:
        output_folder = "{sample}",
        prefix = "{sample}",
        mem = "24000000000"
    log:
        "{sample}.log"
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -o {params.output_folder} --out-prefix {params.prefix} -t {threads} -m {params.mem} --min-contig-len 500 1>{log} 2>&1
        """