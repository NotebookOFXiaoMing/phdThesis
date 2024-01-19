SAMPLES, = glob_wildcards("/data/myan/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/{sample}.Genome.fa.masked")

print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        expand("01.gaf/{sample}.gaf",sample=SAMPLES),
        expand("02.cov/{sample}Cov.tsv",sample=SAMPLES)

rule minigraph_align:
    input:
        qry = "/data/myan/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/{sample}.Genome.fa.masked",
        ref = "/data/myan/raw_data/pome/sour.pome/20231015.reanalysis/18.PanGenomeSize/pome.gfa"
    output:
        "01.gaf/{sample}.gaf"
    threads:
        16
    resources:
        mem_mb = 24000
    shell:
        """
        minigraph -t {threads} --cov -x asm {input.ref} {input.qry} > {output}
        """

def myfun(wildcards):
    if wildcards.sample == "ys":
        return "Y"
    else:
        return "N"


rule step02:
    input:
        rules.minigraph_align.output
    output:
        "02.cov/{sample}Cov.tsv"
    threads:
        4
    resources:
        mem_mb = 8000
    params:
        myfun,
        "{sample}"
    shell:
        """
        python comb_coverage01.py -g {input} -a {params[1]} -o {output} -r {params[0]}
        """
