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

### 基于基因组比对的Indels

```
python filterSyriVCFInDel.py bhagwa ../02.syri.output/bhagwa/bhagwa_syri.vcf bhagwa.InDels.vcf bhagwa.indel.len
python filterSyriVCFInDel.py tunisia ../02.syri.output/tunisia/tunisia_syri.vcf tunisia.InDels.vcf tunisia.indel.len

ls *.vcf > vcf.list
SURVIVOR merge vcf.list 1000 0 1 1 1 0 merged.SURVIVOR.vcf

## 数量
read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%filter(str_sub(X3,1,3)=="INS")
read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%filter(str_sub(X3,1,3)=="DEL")

## 长度
read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%filter(str_sub(X3,1,3)=="INS")%>%mutate(new_col=str_count(X5))%>%pull(new_col)%>%range()
read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%filter(str_sub(X3,1,3)=="DEL")%>%mutate(new_col=str_count(X4))%>%pull(new_col)%>%range()
## 按染色体统计
read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%mutate(vartype=str_sub(X3,1,3))%>%select(X1,vartype) -> indel.count
read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%mutate(vartype=str_sub(X3,1,3))%>%select(X1,vartype)%>%write_csv("indel.count.csv")

dat<-read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/indel.count.csv")
## 分染色体数量的柱形图
dat %>% 
  group_by(X1,vartype) %>% 
  summarise(value_count=n()) %>% 
  ungroup() %>% 
  ggplot(aes(x=X1,y=value_count))+
  geom_bar(stat="identity",aes(fill=vartype))

dat %>% 
  group_by(X1) %>% 
  summarise(value_count=n()) %>%
  ungroup()

read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%filter(str_sub(X3,1,3)=="DEL")%>%mutate(new_col=str_count(X4))%>%pull(new_col)%>%write_lines("del.len")
read_tsv("merged.SURVIVOR.vcf",col_names = FALSE,comment = "#")%>%filter(str_sub(X3,1,3)=="INS")%>%mutate(new_col=str_count(X5))%>%pull(new_col)%>%write_lines("ins.len")


## 不同长度的插入缺失变异数量
read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/ins.len",
         col_names = FALSE) %>% 
  bind_rows(read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/del.len",
                     col_names = FALSE) ) %>% 
  mutate(X2=case_when(
    X1<=100 ~ "A",
    X1>100 & X1<= 1000 ~ "B",
    X1>1000 & X1 <= 10000 ~ "C",
    X1> 10000 ~ "D"
  )) %>% 
  pull(X2) %>% table()

70485/(70485+810+657+159)
159/(70485+810+657+159)

## 相对于基因位置的注释
~/biotools/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq merged.SURVIVOR.vcf > merged.SURVIVOR.annovar.input
~/biotools/annovar/annotate_variation.pl -geneanno --neargene 3000 -buildver genome -dbtype refGene -outfile all.anno -exonsort merged.SURVIVOR.annovar.input ~/my_data/raw_data/practice/annovar/
cat all.anno.variant_function | awk '{print $1}' | sort | uniq -c

```

### 图形泛基因组合并92个样本的vcf

```
snakemake -s bgzip_vcf.smk --cores 24 -p
vcftools --vcf merged.vcf --max-missing 0.8 --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out merged92.vg.filter
#2103
```

### 纯合有利变异

```
read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/merged92.vg.filter.recode.vcf",
         comment = "##") %>% 
  mutate(across(contains("_"),function(x){str_sub(x,1,3)})) -> dat

mut2_freq_list<-list()
for(i in 1:nrow(dat)){
  print(i)
  dat[i,] %>% select(-c(1:9)) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("V2") %>% 
    mutate(V3=str_sub(V2,1,2)) %>% 
    filter(V1 != "./.") %>% 
    mutate(V4=case_when(
      V1 == "0/0" ~ "0",
      V1 == "1/0" | V1 == "0/1" ~ "1",
      V1 == "1/1" ~ "2"
    ))%>% 
    mutate(V5=case_when(
      V3 == "Ti" ~ "Ti",
      TRUE ~ "nonTi"
    )) %>% 
    group_by(V5,V4) %>% 
    summarise(count=n()) %>% 
    ungroup() %>% 
    pivot_wider(names_from = V4,values_from = count) %>% 
    pivot_longer(!V5,values_to = "count",names_to = "V4") %>% 
    mutate(count=replace_na(count,0)) %>% 
    filter(V4==0|V4==2) -> temp.dat
    if(nrow(temp.dat)==4){
      temp.dat %>% 
        pivot_wider(names_from = "V4",values_from = "count") %>% 
        column_to_rownames("V5") %>% 
        fisher.test() %>% .$p.value -> pvalue
      temp.dat %>% 
        group_by(V5) %>% 
        mutate(total=sum(count)) %>% 
        ungroup() %>% 
        mutate(freq=count/total) %>% 
        filter(V4=="2") %>% 
        select(V5,freq) %>% 
        pivot_wider(values_from = freq,names_from = V5) %>% 
        mutate(chrpos=paste(dat[i,1],dat[i,2],sep="_"),
               pval=pvalue) -> mut2_freq_list[[i]]
    }
    

}


mut2_freq_list
mut2_freq_list %>% 
  bind_rows() %>% 
  mutate(pval.adj=p.adjust(pval)) %>% 
  mutate(group=case_when(
    pval.adj > 0.001 ~ "A",
    pval.adj <= 0.001 & pval.adj >= 1e-10 ~ "B",
    pval.adj < 1e-10 ~ "D"
  )) %>% 
  ggplot(aes(x=nonTi,y=Ti))+
  geom_point(aes(color=group))

mut2_freq_list %>% 
  bind_rows() %>% 
  mutate(pval.adj=p.adjust(pval,method = "bonferroni")) %>% 
  mutate(group=case_when(
    pval > 0.001 ~ "A",
    pval <= 0.001 & pval >= 1e-10 ~ "B",
    pval < 1e-10 ~ "D"
  )) %>% 
  filter(pval.adj<=0.001) %>% 
  filter(nonTi>Ti)

mut2_freq_list %>% 
  bind_rows() %>% 
  mutate(pval.adj=p.adjust(pval,method = "bonferroni")) %>% 
  mutate(group=case_when(
    pval > 0.001 ~ "A",
    pval <= 0.001 & pval >= 1e-10 ~ "B",
    pval < 1e-10 ~ "D"
  )) %>% 
  filter(pval.adj<=0.001) %>% 
  filter(nonTi<Ti)
save(mut2_freq_list,file = "D:/Jupyter/panPome/Figures/结构变异图形泛基因组/mut2_freq_list.Rdata")
```

### 杂合

```
mut1_freq_list<-list()
for(i in 1:nrow(dat)){
  print(i)
  dat[i,] %>% select(-c(1:9)) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("V2") %>% 
    mutate(V3=str_sub(V2,1,2)) %>% 
    filter(V1 != "./.") %>% 
    mutate(V4=case_when(
      V1 == "0/0" ~ "0",
      V1 == "1/0" | V1 == "0/1" ~ "1",
      V1 == "1/1" ~ "2"
    ))%>% 
    mutate(V5=case_when(
      V3 == "Ti" ~ "Ti",
      TRUE ~ "nonTi"
    )) %>% 
    group_by(V5,V4) %>% 
    summarise(count=n()) %>% 
    ungroup() %>% 
    pivot_wider(names_from = V4,values_from = count) %>% 
    pivot_longer(!V5,values_to = "count",names_to = "V4") %>% 
    mutate(count=replace_na(count,0)) %>% 
    filter(V4==0|V4==1) -> temp.dat
  if(nrow(temp.dat)==4){
    temp.dat %>% 
      pivot_wider(names_from = "V4",values_from = "count") %>% 
      column_to_rownames("V5") %>% 
      fisher.test() %>% .$p.value -> pvalue
    temp.dat %>% 
      group_by(V5) %>% 
      mutate(total=sum(count)) %>% 
      ungroup() %>% 
      mutate(freq=count/total) %>% 
      filter(V4=="1") %>% 
      select(V5,freq) %>% 
      pivot_wider(values_from = freq,names_from = V5) %>% 
      mutate(chrpos=paste(dat[i,1],dat[i,2],sep="_"),
             pval=pvalue) -> mut1_freq_list[[i]]
  }
  
  
}

save(mut1_freq_list,file = "D:/Jupyter/panPome/Figures/结构变异图形泛基因组/mut1_freq_list.Rdata")
mut1_freq_list %>% length()
mut1_freq_list %>% 
  bind_rows() %>% 
  mutate(pval.adj=p.adjust(pval)) %>% 
  mutate(group=case_when(
    pval.adj > 0.001 ~ "A",
    pval.adj <= 0.001 & pval.adj >= 1e-10 ~ "B",
    pval.adj < 1e-10 ~ "D"
  )) %>% 
  ggplot(aes(x=nonTi,y=Ti))+
  geom_point(aes(color=group))

mut1_freq_list %>% 
  bind_rows() %>% 
  mutate(pval.adj=p.adjust(pval,method = "bonferroni")) %>% 
  mutate(group=case_when(
    pval > 0.001 ~ "A",
    pval <= 0.001 & pval >= 1e-10 ~ "B",
    pval < 1e-10 ~ "D"
  )) %>% 
  filter(pval.adj<=0.001) %>% 
  filter(nonTi>Ti)

mut1_freq_list %>% 
  bind_rows() %>% 
  mutate(pval.adj=p.adjust(pval,method = "bonferroni")) %>% 
  mutate(group=case_when(
    pval > 0.001 ~ "A",
    pval <= 0.001 & pval >= 1e-10 ~ "B",
    pval < 1e-10 ~ "D"
  )) %>% 
  filter(pval.adj<=0.001) %>% 
  filter(nonTi<Ti)
```
