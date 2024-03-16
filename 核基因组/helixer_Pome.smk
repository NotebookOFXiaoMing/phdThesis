SAMPLES, = glob_wildcards("{sample}.fna")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("{sample}.gff",sample=SAMPLES)

def speciesName(wildcards):
    tmp_list = wildcards.sample.split("_")
    return tmp_list[0][0] + tmp_list[1][0:2]

rule helixer:
    input:
        "{sample}.fna"
    output:
        "{sample}.gff"
    threads:
        128
    resources:
        mem_mb = 220000
    params:
        speciesName
    shell:
        """
        singularity run ../helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py \
        --fasta-path {input} --lineage land_plant --species {params} \
        --gff-output-path {output} \
        --model-filepath /home/myan/biotools/Helixer-main/land_plant_model/land_plant_v0.3_a_0080.h5 \
        --subsequence-length 64152
        """