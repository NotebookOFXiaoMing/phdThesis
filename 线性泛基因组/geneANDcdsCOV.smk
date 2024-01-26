SAMPLES, = glob_wildcards("02.sorted.bam/{sample}.sorted.bam")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("03.cdsANDgeneCOV/{sample}/summary_cds.cov",sample=SAMPLES)


rule CalculateCOV:
    input:
        bam = "02.sorted.bam/{sample}.sorted.bam",
        gtf = "ys.final.plus.nonRefSeq.masked.gtf"
    output:
        cov = "03.cdsANDgeneCOV/{sample}/summary_cds.cov",
        sta = "03.cdsANDgeneCOV/{sample}/{sample}.sta"
    threads:
        4
    resources:
        mem_mb = 32000
    shell:
        """
        /home/myan/biotools/software.package/EUPAN-v0.44/bin/ccov \
        {input.gtf} {input.bam} > {output.sta}
        """