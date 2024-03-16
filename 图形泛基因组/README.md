### 统计二代测序数据的数据量
```
snakemake -s fastp_report_stat.smk --cores 2 -p

library(tidyverse)
list.files("01.fastp.report.summary",pattern = "*.csv",full.names = TRUE)%>%map(read_csv)%>%bind_rows()%>%write_csv("total_base.csv")
```



### 二代测序数据检测snp

```
snakemake -s bwaSamtoolsBcftools.smk --configfiles=config_SNP.yaml --cores 128 -p
snakemake -s vcf_merged.smk --cores 128 -p
bgzip merged92.snp.vcf
tabix merged92.snp.vcf.gz
vcftools --gzvcf merged92.snp.vcf.gz --max-missing 0.8 --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out merged92.snp.filter
```

### 统计比对率

```
snakemake -s samtools_flagstat.smk --cores 128 -p
myfun<-function(x){
  dat<-jsonlite::fromJSON(x)
  return(data.frame(sample_id=x,
                    mapping_ratio=dat$`QC-passed reads`$`properly paired %`))
}
list.files("03.samtools.flagstat",pattern = "*.json",full.names = TRUE)%>%map(myfun)%>%bind_rows()%>%write_csv("flagstat.csv")
```

### 二代测序数据检测indels

```
snakemake -s bcftools_call_indels.smk --cores 128 -p
bcftools view -i 'abs(ILEN)<=50' merged92.indel.vcf -O v -o merged92.indel.50.vcf
vcftools --vcf merged92.indel.50vcf --max-missing 0.8 --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out merged92.indel.filter
```

### annovar 注释变异位点位置

```
#library(rtracklayer)
#export(import("ys.rename.gff3"),"ys.rename.gtf","gtf")
gffread -E ys.rename.gff3 -T -o ys.rename.gtf
~/biotools/annovar/gtfToGenePred -genePredExt ys.rename.gtf genome_refGene.txt
~/biotools/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile ys.final.masked.fna genome_refGene.txt --out genome_refGeneMrna.fa # out的文件名必须是这个
~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq ../../merged92.snp.filter.recode.vcf > all.snp.vcf.annovar.input
~/biotools/annovar/annotate_variation.pl -geneanno --neargene 3000 -buildver genome -dbtype refGene -outfile all.anno -exonsort all.snp.vcf.annovar.input ~/my_data/raw_data/practice/annovar/

cat all.anno.variant_function | awk '{print $1}' | sort | uniq -c
cat all.anno.exonic_variant_function | awk '{print $2}' | sort | uniq -c

~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq ../../merged92.indel.filter.recode.vcf > all.indel.vcf.annovar.input
~/biotools/annovar/annotate_variation.pl -geneanno --neargene 3000 -buildver genome -dbtype refGene -outfile all.anno -exonsort all.indel.vcf.annovar.input ~/my_data/raw_data/practice/annovar/
```
