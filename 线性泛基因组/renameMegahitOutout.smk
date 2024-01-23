SAMPLES, = glob_wildcards("02.megahit/{sample}.contigs.fa")

SAMPLES = [SAMPLE.split("/")[0] for SAMPLE in SAMPLES if len(SAMPLE.split("/")) == 2]
#print(SAMPLES)
print("Total sample: ",len(SAMPLES))

rule all:
    input:
        expand("02.megahit.output.rename/{sample}.contigs.fa",sample=SAMPLES)

rule rename:
    input:
        "02.megahit/{sample}/{sample}.contigs.fa"
    output:
        "02.megahit.output.rename/{sample}.contigs.fa"
    threads:
        1
    resources:
        mem = 4000
    params:
        "{sample}"
    run:
        from Bio import SeqIO
        fw = open(output[0],'w')
        for rec in SeqIO.parse(input[0],'fasta'):
            new_id = params[0] + "_" + rec.id
            fw.write(">%s\n%s\n"%(new_id,str(rec.seq)))
        fw.close()