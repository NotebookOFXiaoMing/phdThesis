CHRS = ["chr" + str(i) for i in range(1,9)]
SCAFS = ["scaffold" + str(i).zfill(2) for i in range(1,17)]


rule all:
    input:
        expand("{sample}/{sample}.fna",sample= CHRS + SCAFS),
        expand("{sample}/{sample}.trf",sample= CHRS + SCAFS),
        expand("{sample}/{sample}.trf.filter",sample=CHRS + SCAFS),
        "trf.output.bed"

rule seqkit_grep:
    input:
        "/data/myan/raw_data/pome/sour.pome/20231015.reanalysis/08.proteinCodingGenes/ys.final.masked.fna"
    output:
        "{sample}/{sample}.fna"
    threads:
        8
    resources:
        mem_mb = 16000
    params:
        "{sample}"
    shell:
        """
        seqkit grep -r -p {params} {input} -o {output}
        """

rule trf:
    input:
        rules.seqkit_grep.output
    output:
        "{sample}/{sample}.trf"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        trf {input} 2 7 7 80 10 5 500 -h -d -l 6 -ngs > {output}
        """

rule filter:
    input:
        rules.trf.output
    output:
        "{sample}/{sample}.trf.filter"
    threads:
        8
    resources:
        mem_mb = 8000
    shell:
        """
        python fix_trf_output.py {input} {output}
        """

rule bed:
    input:
        expand(rules.filter.output,sample=CHRS + SCAFS)
    output:
        "trf.output.bed"
    threads:
        8
    resources:
        mem_mb = 16000
    params:
        ','.join(expand(rules.filter.output,sample=CHRS + SCAFS))
    shell:
        """
        python trf_parser.py {params} > {output}
        """