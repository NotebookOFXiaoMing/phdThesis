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
bcftools view -H pome.snp.92.filter.onlyChr.vcf | awk '{print $1"\t"$2"\t"$2+1}' > snp.bed
bcftools view -H pome.indel.92.filter.onlyChr.vcf | awk '{print $1"\t"$2"\t"$2+1}' > indel.bed

bedtools coverage -a gene.sorted.bed -b snp.bed -counts > gene.snp.count
bedtools coverage -a intergenic.bed -b snp.bed -counts > intergenic.snp.count


read_tsv("phdthesis/chapter2/data/gene.snp.count",col_names = FALSE) %>% 
  mutate(X5=X4/(X3-X2)) %>% 
  mutate(X6="gene") %>% 
  bind_rows(
    read_tsv("phdthesis/chapter2/data/intergenic.snp.count",col_names = FALSE) %>% 
      mutate(X5=X4/(X3-X2)) %>% 
      mutate(X6="intergenic")
  ) %>% 
  ggplot(aes(x=X6,y=-log10(X5)))+
  geom_boxplot()

read_tsv("phdthesis/chapter2/data/gene.snp.count",col_names = FALSE) %>% 
  mutate(X5=X4/(X3-X2)) %>% 
  mutate(X6="gene") %>% 
  bind_rows(
    read_tsv("phdthesis/chapter2/data/intergenic.snp.count",col_names = FALSE) %>% 
      mutate(X5=X4/(X3-X2)) %>% 
      mutate(X6="intergenic")
  ) %>%
  t.test(X5~X6,data=.)

bedtools coverage -a gene.sorted.bed -b indel.bed -counts > gene.indel.count
bedtools coverage -a intergenic.bed -b indel.bed -counts > intergenic.indel.count
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

## snp和indel的连锁不平衡

```
## pwd 20240602.reanalysis/01.snp.indel/04.ldR2

plink --vcf ~/biotools/deepVariant/pome/pome.snp.92.filter.onlyChr.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 --vcf-half-call r -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out pome.snp.ld

plink --vcf ~/biotools/deepVariant/pome/pome.indel.92.filter.onlyChr.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 --vcf-half-call r -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out pome.indel.ld

~/anaconda3/envs/py2/bin/python ld_decay_calc.py -i pome.snp.ld.ld.gz -o pome.snp
~/anaconda3/envs/py2/bin/python ld_decay_calc.py -i pome.indel.ld.ld.gz -o pome.indel

## 以上是总的，一下是分群体
vcfsamplenames ~/biotools/deepVariant/pome/pome.snp.92.filter.onlyChr.vcf | grep "Ti" > Ti.samples
vcftools --vcf ~/biotools/deepVariant/pome/pome.snp.92.filter.onlyChr.vcf --max-missing 0.8 --maf 0.05 --keep Ti.samples --recode --recode-INFO-all --out pome.snp.Ti.sample
vcftools --vcf ~/biotools/deepVariant/pome/pome.snp.92.filter.onlyChr.vcf --max-missing 0.8 --maf 0.05 --remove Ti.samples --recode --recode-INFO-all --out pome.snp.nonTi.sample

plink --vcf pome.snp.Ti.sample.recode.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 --vcf-half-call r -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out pome.Ti.sample
~/anaconda3/envs/py2/bin/python ld_decay_calc.py -i pome.Ti.sample.ld.gz -o pome.Ti.sample

plink --vcf pome.snp.nonTi.sample.recode.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 --vcf-half-call r -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out pome.nonTi.sample

~/anaconda3/envs/py2/bin/python ld_decay_calc.py -i pome.nonTi.sample.ld.gz -o pome.nonTi.sample

~/biotools/PopLDdecay-3.42/bin/PopLDdecay -InVCF pome.snp.Ti.sample.recode.vcf -OutStat pome.snp.Ti.LDdecay
perl ~/biotools/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile pome.snp.Ti.LDdecay.stat.gz -output pome.snp.Ti.LDdecay.Fig

~/biotools/PopLDdecay-3.42/bin/PopLDdecay -InVCF pome.snp.nonTi.sample.recode.vcf -OutStat pome.snp.nonTi.LDdecay
perl ~/biotools/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile pome.snp.nonTi.LDdecay.stat.gz -output pome.snp.nonTi.LDdecay.Fig

## Ti地区的ld有点大，每个地点只用一个样本试试
cat Ti.samples | grep "01" > Ti.samples.01 #9个
vcftools --vcf pome.snp.Ti.sample.recode.vcf --keep Ti.samples.01 --recode --recode-INFO-all --out pome.snp.Ti.sample.01
plink --vcf pome.snp.Ti.sample.01.recode.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 --vcf-half-call r -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out pome.Ti.sample.01

## 这个也差不多

## Indel

vcftools --vcf ~/biotools/deepVariant/pome/pome.indel.92.filter.onlyChr.vcf --keep ../01.snp/Ti.samples --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out pome.indel.Ti.sample
vcftools --vcf ~/biotools/deepVariant/pome/pome.indel.92.filter.onlyChr.vcf --remove ../01.snp/Ti.samples --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out pome.indel.nonTi.sample

plink --vcf pome.indel.nonTi.sample.recode.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 --vcf-half-call r -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out pome.nonTi.sample.indel

plink --vcf pome.indel.Ti.sample.recode.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 --vcf-half-call r -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out pome.Ti.sample.indel

python 20240524_01.py pome.indel.Ti.sample.recode.vcf pome.indel.Ti.sample.recode.edited.vcf ## 把alt ref 改成单碱基的形式
python 20240524_01.py pome.indel.nonTi.sample.recode.vcf pome.indel.nonTi.sample.recode.edited.vcf


~/biotools/PopLDdecay-3.42/bin/PopLDdecay -InVCF pome.indel.Ti.sample.recode.vcf -OutStat pome.indel.Ti.LDdecay
perl ~/biotools/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile pome.indel.Ti.LDdecay.stat.gz -output pome.indel.Ti.LDdecay.Fig

~/biotools/PopLDdecay-3.42/bin/PopLDdecay -InVCF pome.indel.nonTi.sample.recode.vcf -OutStat pome.indel.nonTi.LDdecay
perl ~/biotools/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile pome.indel.nonTi.LDdecay.stat.gz -output pome.indel.nonTi.LDdecay.Fig

## snp 和 indel的连锁不平恒

python 20240524_02.py ~/biotools/deepVariant/pome/pome.snp.92.filter.onlyChr.vcf snp.edited.vcf
python ../02.indel/20240524_01.py ~/biotools/deepVariant/pome/pome.indel.92.filter.onlyChr.vcf indel.edited.vcf

bgzip snp.edited.vcf
tabix snp.edited.vcf.gz

bgzip indel.edited.vcf
tabix indel.edited.vcf.gz

vcfcat snp.edited.vcf.gz indel.edited.vcf.gz > merged.snp.indel.vcf
vcfsort merged.snp.indel.vcf > merged.snp.indel.sorted.vcf

plink --vcf merged.snp.indel.sorted.vcf --double-id --allow-extra-chr --maf 0.05 --geno 0.1 --mind 0.5 --thin 0.4 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --vcf-half-call r --out merged.snp.indel

```

## PCA 和 进化树
```
plink --vcf ~/biotools/deepVariant/pome/pome.snp.92.filter.onlyChr.vcf --recode12 --allow-extra-chr --allow-no-sex --vcf-half-call r --out pome.snp
plink --allow-extra-chr --file pome.snp --indep-pairwise 50 5 0.1 --recode vcf-iid --out pome.snp.LDpruned
plink --allow-extra-chr --file pome.snp --recode vcf-iid --extract pome.snp.LDpruned.prune.in --out pome.snp.LDpruned

~/biotools/VCF2PCACluster-1.40/bin/VCF2PCACluster -InVCF pome.snp.LDpruned.vcf -OutPut pome.snp.PCA

plink --vcf ~/biotools/deepVariant/pome/pome.indel.92.filter.onlyChr.vcf --recode12 --allow-extra-chr --allow-no-sex --vcf-half-call r --out pome.indel
plink --allow-extra-chr --file pome.indel --indep-pairwise 50 5 0.1 --recode vcf-iid --out pome.indel.LDpruned
plink --allow-extra-chr --file pome.indel --recode vcf-iid --extract pome.indel.LDpruned.prune.in --out pome.indel.LDpruned

~/biotools/VCF2PCACluster-1.40/bin/VCF2PCACluster -InVCF pome.indel.LDpruned.vcf -OutPut pome.indel.PCA
```

## 核苷酸多样性
```
vcfsamplenames ../01.snp/pome.snp.LDpruned.vcf | tail -n 36 > Ti.samples
vcfsamplenames ../01.snp/pome.snp.LDpruned.vcf | head -n 56 > nonTi.samples

vcftools --vcf ../01.snp/pome.snp.LDpruned.vcf --keep Ti.samples --window-pi 100000 --out Ti.pi
vcftools --vcf ../01.snp/pome.snp.LDpruned.vcf --keep nonTi.samples --window-pi 100000 --out nonTi.pi
```

## 全基因组关联分析

```
## 全基因组关联分析的表型整理

read_tsv("phdthesis/chapter2/tables/table1.txt") %>% 
  dplyr::select(Taxa,`单果重`) %>% 
  mutate(Taxa=str_replace(Taxa,"_",""),
         `单果重`=str_split_fixed(`单果重`,"\\+",2) %>% 
           as.data.frame() %>% 
           pull(V1)) %>% 
  mutate(newcol=Taxa) %>% 
  dplyr::select(1,3,2) %>% 
  write_tsv("phdthesis/chapter2/data/gwas.pheno/01.fruit.weight.txt",col_names = FALSE)



bcftools view -S 26.samples.with.pheno.txt ~/biotools/deepVariant/pome/pome.snp.92.filter.onlyChr.vcf | grep -v "scaffold" > pome.snp.26.filter.onlyChr.vcf
cat 26.samples.with.pheno.txt | awk 'gsub("_","")' > 26.samples.with.pheno.new.name

~/biotools/software.package/bcftools-1.17/bcftools reheader -s 26.samples.with.pheno.new.name pome.snp.26.filter.onlyChr.vcf > pome.snp.26.filter.onlyChr.new.name.vcf

## snp显著性于阈值
plink --vcf ../pome.snp.26.filter.onlyChr.new.name.vcf --make-bed --vcf-half-call r --out pome.snp.26
java -jar -Xmx8g ~/biotools/GEC/gec/gec.jar --effect-number --plink-binary pome.snp.26 --genome --out pome.snp.26


plink --vcf ../pome.snp.26.filter.onlyChr.new.name.vcf --allow-extra-chr --make-bed --vcf-half-call r --out pome.snp.26
plink --bfile pome.snp.26 --recode 12 transpose --out pome.snp.26

~/biotools/emmax/emmax-kin-intel64 -v -d 10 pome.snp.26

### 单果重关联
~/biotools/emmax/emmax-intel64 -v -d 10 -t ../../02.emmex/pome.snp.26 -k ../../02.emmex/pome.snp.26.aBN.kinf -p ../../../00.pheno/01.fruit.weight.txt -o 01.fruit.weigth
Rscript filterSignificSites.R 01.fruit.weigth.ps 01.fruit.weigth.Signifi.ID
bcftools view -i 'ID=@01.fruit.weigth.Signifi.ID' ../../pome.snp.26.filter.onlyChr.new.name.vcf > pome.snp.26.01.ruit.weight.Signifi.vcf

~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq pome.snp.26.01.ruit.weight.Signifi.vcf > annovar.input
~/biotools/annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver genome -dbtype refGene -outfile annovar.output -exonsort annovar.input ~/biotools/deepVariant/pome/annovar.output/

### 百粒重关联
~/biotools/emmax/emmax-intel64 -v -d 10 -t ../../02.emmex/pome.snp.26 -k ../../02.emmex/pome.snp.26.aBN.kinf -p ../../../00.pheno/02.hundred.seed.weight.txt -o 02.hundred.seed.weight
Rscript ../01.fruit.weight/filterSignificSites.R 02.hundred.seed.weight.ps 02.hundred.seed.weight.Signifi.ID
bcftools view -i 'ID=@02.hundred.seed.weight.Signifi.ID' ../../pome.snp.26.filter.onlyChr.new.name.vcf > pome.snp.26.02.hundred.seed.weight.Signifi.vcf

~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq pome.snp.26.02.hundred.seed.weight.Signifi.vcf > annovar.input
~/biotools/annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver genome -dbtype refGene -outfile annovar.output -exonsort annovar.input ~/biotools/deepVariant/pome/annovar.output/

### 果皮厚度关联
~/biotools/emmax/emmax-intel64 -v -d 10 -t ../../02.emmex/pome.snp.26 -k ../../02.emmex/pome.snp.26.aBN.kinf -p ../../../00.pheno/03.peelthickness.txt -o 03.peelthickness
Rscript ../01.fruit.weight/filterSignificSites.R 03.peelthickness.ps 03.peelthickness.Signifi.ID
bcftools view -i 'ID=@03.peelthickness.Signifi.ID' ../../pome.snp.26.filter.onlyChr.new.name.vcf > pome.snp.26.03.peelthickness.Signifi.vcf

~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq pome.snp.26.03.peelthickness.Signifi.vcf > annovar.input
~/biotools/annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver genome -dbtype refGene -outfile annovar.output -exonsort annovar.input ~/biotools/deepVariant/pome/annovar.output/

### 籽粒硬度相关

~/biotools/emmax/emmax-intel64 -v -d 10 -t ../../02.emmex/pome.snp.26 -k ../../02.emmex/pome.snp.26.aBN.kinf -p ../../../00.pheno/04.seedhardness.txt -o 04.seedhardness
Rscript ../01.fruit.weight/filterSignificSites.R 04.seedhardness.ps 04.seedhardness.Signifi.ID

### 可滴定酸
~/biotools/emmax/emmax-intel64 -v -d 10 -t ../../02.emmex/pome.snp.26 -k ../../02.emmex/pome.snp.26.aBN.kinf -p ../../../00.pheno/05.titratableacidity.txt -o 05.titratableacidity
Rscript ../01.fruit.weight/filterSignificSites.R 05.titratableacidity.ps 05.titratableacidity.Signifi.ID

bcftools view -i 'ID=@05.titratableacidity.Signifi.ID' ../../pome.snp.26.filter.onlyChr.new.name.vcf > pome.snp.26.05.titratableacidity.Signifi.vcf
~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq pome.snp.26.05.titratableacidity.Signifi.vcf > annovar.input
~/biotools/annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver genome -dbtype refGene -outfile annovar.output -exonsort annovar.input ~/biotools/deepVariant/pome/annovar.output/



## indel关联分析

### 先把26个样本挑出来，然后更改样本名，因为样本名里不能有下划线

bcftools view -S ../01.snp/26.samples.with.pheno.txt ~/biotools/deepVariant/pome/pome.indel.92.filter.onlyChr.vcf | grep -v 'scaffold' > pome.indel.26.filter.onlyChr.vcf
~/biotools/software.package/bcftools-1.17/bcftools reheader -s ../01.snp/26.samples.with.pheno.new.name pome.indel.26.filter.onlyChr.vcf > pome.indel.26.filter.onlyChr.new.name.vcf

## 显著性阈值
plink --vcf ../pome.indel.26.filter.onlyChr.new.name.vcf --make-bed --vcf-half-call r --out pome.indel.26
java -jar -Xmx8g ~/biotools/GEC/gec/gec.jar --effect-number --plink-binary pome.indel.26 --genome --out pome.indel.26

## emmax 前期准备

plink --vcf ../pome.indel.26.filter.onlyChr.new.name.vcf --allow-extra-chr --make-bed --vcf-half-call r --out pome.indel.26
plink --bfile pome.indel.26 --recode 12 transpose --out pome.indel.26

~/biotools/emmax/emmax-kin-intel64 -v -d 10 pome.indel.26

## GWAS关联 indel 写个批量吧

snakemake -s indelGWAS.smk --cores 8 -p


```
