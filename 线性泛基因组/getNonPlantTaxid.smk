SAMPLES, = glob_wildcards("blastn.output/{sample}.ntblast")

print("total sample size: ",len(SAMPLES))

plant_taxid = []
with open("/data/myan/database/nt/plant.taxid.txt",'r') as fr:
    for line in fr:
        plant_taxid.append(line.strip())

print(plant_taxid[0:5])


rule all:
    input:
        expand("nonPlantTaxid/{sample}.taxid",sample=SAMPLES)

rule rmDuplicate:
    input:
        "blastn.output/{sample}.ntblast"
    output:
        "blastn.output.rmDup/{sample}.ntblast"
    threads:
        4
    resources:
        mem_mb = 16000
    shell:
        """
        cat {input} | awk -v FS="\t" '{{print $1"\t"$15}}' | sort -u > {output}
        """

rule getNotPlantTaxid:
    input:
        rules.rmDuplicate.output
    output:
        "nonPlantTaxid/{sample}.taxid"
    threads:
        4
    resources:
        mem_mb = 96000
    run:
        fw = open(output[0],'w')
        with open(input[0],'r') as fr:
            for line in fr:
                if line.strip().split("\t")[-1] not in plant_taxid:
                    fw.write('%s\n'%(line.strip().split("\t")[0]))

        fw.close()