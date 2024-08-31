## 表型数据整理

```
library(readxl)

read_excel("D:/000博士毕业论文/表型数据/01.单果重.xlsx") %>% 
  group_by(X1) %>% 
  summarise(singleFruitWeight=mean(重量,na.rm = TRUE),
            sd=sd(重量)) %>% 
  rename("Taxa"="X1") %>% 
  filter(Taxa!="Ch_TNS") %>% 
  mutate(单果重=paste(as.character(round(singleFruitWeight,2)),
                   as.character(round(sd,2)),
                   sep = "+")) %>% 
  select(1,4) -> tmp.dat.01

read_excel("D:/000博士毕业论文/表型数据/03.果皮厚度mm.xlsx") %>% 
  group_by(X1) %>% 
  summarise(FruitPeelThickness=mean(果皮厚mm,na.rm = TRUE),
            sd=sd(果皮厚mm)) %>% 
  rename("Taxa"="X1") %>% 
  filter(Taxa!="Ch_TNS") %>%
  mutate(果皮厚度=paste(as.character(round(FruitPeelThickness,2)),
                   as.character(round(sd,2)),
                   sep = "+")) %>% 
  select(1,4)-> tmp.dat.02

read_excel("D:/000博士毕业论文/表型数据/04.百粒重.xlsx") %>% 
  group_by(X1) %>% 
  summarise(HundredSeedWeight=mean(百粒重,na.rm = TRUE),
            sd=sd(百粒重)) %>% 
  rename("Taxa"="X1") %>% 
  filter(Taxa!="Ch_TNS") %>%
  mutate(百粒重=paste(as.character(round(HundredSeedWeight,2)),
                   as.character(round(sd,2)),
                   sep = "+")) %>% 
  select(1,4) -> tmp.dat.03

read_excel("D:/000博士毕业论文/表型数据/05.籽粒硬度.xlsx") %>% 
  group_by(X1) %>% 
  summarise(SeedHardness=mean(`kg/cm2`,na.rm = TRUE),
            sd=sd(`kg/cm2`)) %>% 
  rename("Taxa"="X1") %>% 
  filter(Taxa!="Ch_TNS") %>%
  mutate(籽粒硬度=paste(as.character(round(SeedHardness,2)),
                   as.character(round(sd,2)),
                   sep = "+")) %>% 
  select(1,4) -> tmp.dat.04

read_excel("D:/000博士毕业论文/表型数据/09.可滴定酸.xlsx")  %>% 
  group_by(X1) %>% 
  summarise(TitratableAcidity=mean(`可滴定酸（%）`,na.rm = TRUE),
            sd=sd(`可滴定酸（%）`)) %>% 
  rename("Taxa"="X1") %>% 
  filter(Taxa!="Ch_TNS") %>%
  mutate(可滴定酸含量=paste(as.character(round(TitratableAcidity,2)),
                    as.character(round(sd,2)),
                    sep = "+")) %>% 
  select(1,4) -> tmp.dat.05


tmp.dat.01 %>% 
  left_join(tmp.dat.03,by=c("Taxa"="Taxa"))%>% 
  left_join(tmp.dat.02,by=c("Taxa"="Taxa"))%>% 
  left_join(tmp.dat.04,by=c("Taxa"="Taxa"))%>% 
  left_join(tmp.dat.05,by=c("Taxa"="Taxa")) %>% 
  write_tsv("phdthesis/chapter2/tables/table1.txt")
```

## 表型数据的相关性

```
library(GGally)

ggpairs(pheno.dat,columns = 1:5,
        upper = list(continuous=wrap('cor',method="pearson")))+
  theme_bw()+
  theme(panel.grid = element_blank())

```

## 表型数据主成分分析

```
pheno.pca<-prcomp(pheno.dat,scale. = TRUE)

pheno.pca %>% 
  summary()

pheno.pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column("sampleid") %>% 
  mutate(group=str_sub(sampleid,1,2)) %>% 
  ggplot(aes(x=PC1,y=PC2))+
  geom_point(aes(color=group),size=5)+
  stat_ellipse(aes(x=PC1,y=PC2,fill=group),
               geom="polygon",
               alpha=0.2,
               lty="dashed",
               color="black",
               show.legend = FALSE)+
  labs(x="PC1 57.35%",y="PC2 22.36%")+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.1,0.8))+
  scale_fill_manual(values = c("#009f73","#0072b1"))+
  scale_color_manual(values = c("#009f73","#0072b1")) 
```

## snp和indels检测

```
bwa index ref/ys.final.masked.fna
sbatch run_bwa.slurm
sbatch run_picard.slurm
sbatch run_deepVariant.slurm

## 数据存储路径 biotools/deepVariant/pome
singularity run ../glnexus_v1.4.1.sif glnexus_cli --config DeepVariantWGS 02.g.vcf/*.vcf > pome.bcf
bcftools view pome.bcf -O v -o pome.vcf

bcftools view pome.vcf -v snps -H | wc -l ## snp位点的总数量 1423513
bcftools view pome.vcf -V snps -H | wc -l ## indel位点的总数量 279731

bcftools view pome.vcf -v snps -m2 -M2 -O v -o pome.snp.vcf

vcftools --vcf pome.snp.vcf --remove remove.4.samples --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out pome.snp.92.filter
## After filtering, kept 954669 out of a possible 1399357 Sites

bcftools view pome.vcf -V snps -m2 -M2 -e 'strlen(REF)>=50 || strlen(ALT)>=50' -O v -o pome.indel.vcf
vcftools --vcf pome.indel.vcf --remove remove.4.samples --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out pome.indel.92.filter

## After filtering, kept 116292 out of a possible 212752 Sites

bgzip -c pome.snp.92.filter.recode.vcf > pome.snp.92.filter.recode.vcf.gz
tabix pome.snp.92.filter.recode.vcf.gz

bgzip -c pome.indel.92.filter.recode.vcf > pome.indel.92.filter.recode.vcf.gz
tabix pome.indel.92.filter.recode.vcf.gz

bcftools view -H pome.snp.92.filter.recode.vcf.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8 | wc -l

bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8 pome.snp.92.filter.recode.vcf.gz -O v -o pome.snp.92.filter.onlyChr.vcf
bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8 pome.indel.92.filter.recode.vcf.gz -O v -o pome.indel.92.filter.onlyChr.vcf

bcftools stats pome.snp.92.filter.onlyChr.vcf # 转换和颠换的比例
```

## 计算基因区间和基因间区 的snp数量

```
grep "gene" ~/my_data/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/ys.rename.gff3 | awk '{print $1"\t"$4"\t"$5}' > gene.bed
cat ~/my_data/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/ys.Genome.fa.masked.fai | awk '{print $1"\t"$2}' > chr.len.bed

bedtools sort -i gene.bed -g chr.len.bed > gene.sorted.bed
bedtools complement -i gene.sorted.bed -g chr.len.bed > intergenic.bed
```

## annovar 变异位点相对于基因的位置

```
conda activate syri
conda install ucsc-gtftogenepred
conda install ucsc-gff3togenepred

gffread ~/my_data/raw_data/pome/sour.pome/20231015.reanalysis/08.proteinCodingGenes/05.evm/ys/ys.rename.gff3 -T -o ys.gtf
gtfToGenePred -genePredExt ys.gtf genome_refGene.txt
~/biotools/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile ~/my_data/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/ys.Genome.fa.masked --out genome_refGeneMrna.fa genome_refGene.txt

~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq ../pome.snp.92.filter.recode.vcf > pome.snp.92.filter.annovar.input

#NOTICE: Finished reading 954714 lines from VCF file
#NOTICE: A total of 954669 locus in VCF file passed QC threshold, representing 954667 SNPs (634892 transitions and 319775 transversions) and 2 indels/substitutions
#NOTICE: Finished writing allele frequencies based on 87829364 SNP genotypes (58410064 transitions and 29419300 transversions) and 184 indels/substitutions for 92 samples

 ~/biotools/annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver genome -dbtype refGene -outfile snp.anno -exonsort pome.snp.92.filter.annovar.input ./
## -buildver 后面跟的内容就是genome_refGene.txt 这个文件的前缀名


~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq ../pome.indel.92.filter.recode.vcf > pome.indel.92.filter.annovar.input

#NOTICE: Finished reading 116337 lines from VCF file
#NOTICE: A total of 116292 locus in VCF file passed QC threshold, representing 0 SNPs (0 transitions and 0 transversions) and 116292 indels/substitutions
#NOTICE: Finished writing allele frequencies based on 0 SNP genotypes (0 transitions and 0 transversions) and 10698864 indels/substitutions for 92 samples

~/biotools/annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver genome -dbtype refGene -outfile indel.anno -exonsort pome.indel.92.filter.annovar.input ./
```

## snpEff 变异位点相对于基因的位置

```
java -Xmx4G -jar snpEff.jar build -gtf22 ys

```

## 分染色体计算snp和indel的数量

```
library(data.table)
library(tidyverse)

fread("pome.snp.92.filter.recode.vcf",skip = "#CHROM") -> snp.dat



read_tsv("~/my_data/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/ys.Genome.fa.masked.fai",col_names = FALSE)%>%filter(str_starts(X1,"chr"))%>%select(1,2)%>%magrittr::set_colnames(c("X1","chrlen")) -> chr.len

snp.dat%>%select(`#CHROM`)%>%group_by(`#CHROM`)%>%summarise(snpcount=n())%>%filter(str_starts(`#CHROM`,"chr")) -> snp.count

fread("pome.indel.92.filter.recode.vcf",skip = "#CHROM") -> indel.dat
indel.dat%>%select(`#CHROM`)%>%group_by(`#CHROM`)%>%summarise(indelcount=n())%>%filter(str_starts(`#CHROM`,"chr")) -> indel.count

chr.len%>%left_join(snp.count,by=c("X1"="#CHROM"))%>%left_join(indel.count,by=c("X1"="#CHROM"))%>%write_tsv("chrlen_snp_indel_num.txt")
```
