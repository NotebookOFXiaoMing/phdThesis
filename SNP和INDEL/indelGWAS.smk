SAMPLES, = glob_wildcards("/data/myan/raw_data/pome/20240602.reanalysis/01.snp.indel/06.gwas/00.pheno/{sample}.txt")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("{sample}/{sample}.ps",sample=SAMPLES)


rule gwas:
    input:
        pheno = "/data/myan/raw_data/pome/20240602.reanalysis/01.snp.indel/06.gwas/00.pheno/{sample}.txt",
        knif = "/data/myan/raw_data/pome/20240602.reanalysis/01.snp.indel/06.gwas/02.indel/02.emmex/pome.indel.26.aBN.kinf",
        vcf = "/data/myan/raw_data/pome/20240602.reanalysis/01.snp.indel/06.gwas/02.indel/pome.indel.26.filter.onlyChr.new.name.vcf"
    output:
        ps = "{sample}/{sample}.ps",
        signifi_id = "{sample}/{sample}.Signifi.ID"
    threads:
        8
    params:
        p1 = "/data/myan/raw_data/pome/20240602.reanalysis/01.snp.indel/06.gwas/02.indel/02.emmex/pome.indel.26",
        p2 = "{sample}",
        p3 = "{sample}.ps",
        p4 = "{sample}.Signifi.ID",
        p5 = "pome.indel.26.{sample}.Signifi.vcf"
    resources:
        mem_mb = 24000
    shell:
        """
        #mkdir {params.p2}
        cd {params.p2}
        ~/biotools/emmax/emmax-intel64 -v -d 10 -t {params.p1} \
        -k {input.knif} -p {input.pheno} -o {params.p2}

        Rscript ../../../01.snp/03.gwas.output/01.fruit.weight/filterSignificSites.R \
        {params.p3} {params.p4}
        bcftools view -i 'ID=@{params.p4}' {input.vcf} > {params.p5}
        ~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample \
        -withfreq {params.p5} > annovar.input
        ~/biotools/annovar/annotate_variation.pl -geneanno --neargene 2000 \
        -buildver genome -dbtype refGene -outfile annovar.output \
        -exonsort annovar.input ~/biotools/deepVariant/pome/annovar.output/
        """