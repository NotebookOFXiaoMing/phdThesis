SRR, = glob_wildcards("clean.fastq/{srr}_2.fq")

print("Total sample size: ",len(SRR))


rule all:
    input:
        expand("03.stringtie.output/{srr}/{srr}.gtf",srr=SRR)


rule stringtie:
    input:
        "02.alignment/{srr}/{srr}_Aligned.sortedByCoord.out.bam"
    output:
        gtf = "03.stringtie.output/{srr}/{srr}.gtf",
        gene_abund = "03.stringtie.output/{srr}/{srr}_gene_abund.tab"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        stringtie {input} -o {output.gtf} -A {output.gene_abund} \
        -a 5 -j 0.1 -f 0.3 -g 50 -M 1 -c 2.5 -p {threads}
        """
