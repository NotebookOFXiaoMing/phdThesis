SAMPLES, = glob_wildcards("02.megahit.output.rename/{sample}.contigs.fa")
SAMPLES.remove("Zpt")
SAMPLES.remove("Ycdbft")
print("Total sample: ",len(SAMPLES))

rule all:
    input:
        expand("04.unalignContigs/{sample}.unalignContig.fa",sample=SAMPLES)

rule extractUnalignContig:
    input:
        "03.quast/{sample}/contigs_reports/contigs_report_{sample}-contigs.unaligned.info",
        "02.megahit.output.rename/{sample}.contigs.fa"
    output:
        "04.unalignContigs/{sample}.unalignContig.fa"
    threads:
        4
    resources:
        mem_mb = 16000
    run:
        from Bio import SeqIO
        i = 0
        contig_list = []

        with open(input[0],'r') as fr:
            for line in fr:
                i += 1
                if i > 1:
                    if line.strip().split()[3] == "full":
                        contig_list.append(line.strip().split()[0])
                    elif line.strip().split()[3] == "partial":
                        Unaligned_length = int(line.strip().split()[2])
                        Total_length = int(line.strip().split()[1])

                        if Unaligned_length/Total_length >= 0.5:
                            contig_list.append(line.strip().split()[0])

        #print(contig_list)

        fw = open(output[0],'w')

        for rec in SeqIO.parse(input[1],'fasta'):
            
            if rec.id in contig_list:
                #print(rec.id)
                fw.write(">%s\n%s\n"%(rec.id,str(rec.seq)))

        fw.close()
