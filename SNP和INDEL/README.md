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


```

## 变异位点相对于基因的位置

```

```
