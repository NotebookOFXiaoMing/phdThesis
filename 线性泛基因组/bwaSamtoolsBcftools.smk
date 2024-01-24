#SAMPLES, = glob_wildcards("/data/yewx/grape_new/NGS_public/{sample}_1.fastq.gz")

SAMPLES, = glob_wildcards(config['fq']['r1'])



print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        expand("05.vcf/{sample}.vcf",sample=SAMPLES)


rule bwa:
    input:
        r1 = config['fq']['r1'],
        r2 = config['fq']['r2'],
        ref = config['ref']['fa']
    output:
        "01.sam/{sample}.sam"
    threads:
        32
    resources:
        mem_mb = 96000
    shell:
        """
        bwa mem -t {threads} -M {input.ref} {input.r1} {input.r2} -o {output}
        """

rule samtools_sort:
    input:
        rules.bwa.output
    output:
        "02.sorted.bam/{sample}.sorted.bam"
    threads:
        16
    resources:
        mem_mb = 64000
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

rule picard:
    input:
        rules.samtools_sort.output
    output:
        bam = "03.sorted.md.bam/{sample}.sorted.md.bam",
        txt = "03.sorted.md.bam/{sample}.sorted_dup_metrics.txt"
    threads:
        32
    resources:
        mem_mb = 64000
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} Remove_Duplicates=true ASSUME_SORTED=true M={output.txt}
        """

rule bcftools_mpileup:
    input:
        bam = rules.picard.output.bam,
        ref = config['ref']['fa']
    output:
        "04.bcf/{sample}.bcf"
    threads:
        32
    resources:
        mem_mb = 96000
    shell:
        """
        bcftools mpileup -q 20 -Q 20 -C 50 --fasta-ref {input.ref} {input.bam} -O u -o {output} --threads {threads}
        """

rule bcftools_call:
    input:
        rules.bcftools_mpileup.output
    output:
        "05.vcf/{sample}.vcf"
    threads:
        32
    resources:
        mem_mb = 64000
    shell:
        """
        bcftools call -c -V indels -v {input} -O v -o {output} --threads {threads}
        """