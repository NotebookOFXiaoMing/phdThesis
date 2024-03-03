SRR, = glob_wildcards("clean.fastq/{srr}_2.fq")

#SRR = ["SRR5279386","SRR5279387"]
print("Total sample size: ",len(SRR))

import os 

rule all:
    input:
        expand(os.path.abspath("01.alignment/{srr}/{srr}_Aligned.sortedByCoord.out.bam"),srr=SRR)

rule star01:
    input:
        r1 = os.path.abspath("clean.fastq/{srr}_1.fq"),
        r2 = os.path.abspath("clean.fastq/{srr}_2.fq")
    output:
        os.path.abspath("01.alignment/{srr}/{srr}_Aligned.sortedByCoord.out.bam")
    threads:
        16
    resources:
        mem_mb = 24000
    params:
        ref_index = os.path.abspath("ref.index"),
        prefix = "{srr}_",
        output_folder = "01.alignment/{srr}"
    shell:
        """
        cd {params.output_folder}
        STAR --runThreadN {threads} \
        --readFilesIn {input.r1} {input.r2} \
        --genomeDir {params.ref_index} \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.prefix} \
        --outSAMunmapped Within
        """