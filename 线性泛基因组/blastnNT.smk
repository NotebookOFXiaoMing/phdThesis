SAMPLES, = glob_wildcards("split.seqs/{sample}.fa")

print("total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("blastn.output/{sample}.ntblast",sample=SAMPLES)

rule blastn:
    input:
        "split.seqs/{sample}.fa"
    output:
        "blastn.output/{sample}.ntblast"
    threads:
        13
    resources:
        mem_mb = 60000
    params:
        "/home/myan/my_data/database/nt/ntdb"
    shell:
        """
        blastn -query {input} -db {params} -out {output} \
        -evalue 1e-5 -perc_identity 0.8 \
        -task megablast -outfmt '6 std qcovs stitle staxid' \
        -max_target_seqs 5 -num_threads {threads}
        """