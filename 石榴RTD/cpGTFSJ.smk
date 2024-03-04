SRR, = glob_wildcards("clean.fastq/{srr}_2.fq")

#SRR = ["SRR5279386","SRR5279387"]
print("Total sample size: ",len(SRR))


rule all:
    input:
        expand("04.assemblies/{srr}.gtf",srr=SRR)


rule cp:
    input:
        gtf = "03.stringtie.output/{srr}/{srr}.gtf",
        SJ = "02.alignment/{srr}/{srr}_SJ.out.tab"
    output:
        gtf = "04.assemblies/{srr}.gtf",
        SJ = "04.SJdata/{srr}_SJ.out.tab"
    threads:
        4
    resources:
        mem_mb = 16000
    shell:
        """
        cp {input.gtf} {output.gtf}
        cp {input.SJ} {output.SJ}
        """