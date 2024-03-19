setwd("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/")
library(tidyverse)
library(clusterProfiler)
library(readxl)
library(GGally)

read_tsv("bhagwa_inversion.txt",col_names = FALSE) %>% 
  mutate(X4=X3-X2+1,
         X8=X7-X6+1) %>% 
  write_csv("bhagwa_倒位变异.csv",col_names = FALSE)

read_tsv("bhagwa_translocation.txt",col_names = FALSE) %>% 
  mutate(X4=X3-X2+1,
         X8=X7-X6+1) %>% 
  write_csv("bhagwa_易位变异.csv",col_names = FALSE)
  

read_tsv("tunisia_inversion.txt",col_names = FALSE) %>% 
  mutate(X4=X3-X2+1,
         X8=X7-X6+1) %>% 
  write_csv("tunisia_倒位变异.csv",col_names = FALSE)


read_tsv("tunisia_translocation.txt",col_names = FALSE) %>% 
  mutate(X4=X3-X2+1,
         X8=X7-X6+1) %>% 
  write_csv("tunisia_易位变异.csv",col_names = FALSE)


## 二代测序数据数据量
read_csv("D:/Jupyter/panPome/total_base.csv") %>% 
  pull(r1) %>% sum() -> r1.total

read_csv("D:/Jupyter/panPome/total_base.csv") %>% 
  pull(r2) %>% sum() -> r2.total

(r1.total + r2.total)/1000000000

403043829841/1000000000

read_csv("D:/Jupyter/panPome/total_base.csv") %>% 
  mutate(depth=(r1+r2)/(337.22*1000000)) %>% 
  pull(depth) %>% mean()

read_csv("D:/Jupyter/panPome/total_base.csv") %>% 
  mutate(depth=(r1+r2)/(337.22*1000000)) %>% 
  #filter(depth<20) %>% 
  pull(depth) %>% range()
library(viridis)
library(extrafont)
fonts()

## 测序深度的图
cairo_pdf("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/重测序数据的测序深度.pdf",
    width=12,height = 8)
read_csv("D:/Jupyter/panPome/total_base.csv") %>% 
  mutate(depth=(r1+r2)/(337.22*1000000)) %>%
  ggplot(aes(x=sample_id,y=depth))+
  geom_segment(aes(x=sample_id,xend=sample_id,y=0,yend=depth),
               color="gray")+
  geom_point(aes(color=depth),size=5)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,hjust = 1,vjust=1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.title.y = element_text(size=20,family="KaiTi"),
        axis.text.y = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,16))+
  geom_hline(yintercept = 12.99,lty="dashed")+
  labs(x=NULL,y="测序深度")+
  scale_color_viridis_c(direction = -1)
dev.off()

## 比对率的图

read_csv("D:/Jupyter/panPome/flagstat.csv") %>% 
  pull(mapping_ratio) %>% mean()
read_csv("D:/Jupyter/panPome/flagstat.csv") %>% 
  pull(mapping_ratio) %>% range()
cairo_pdf("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/重测序数据的比对率.pdf",
          width=12,height = 8)
read_csv("D:/Jupyter/panPome/flagstat.csv") %>% 
  mutate(sample_id=str_replace(sample_id,pattern = "03.samtools.flagstat/","") %>% 
           str_replace(pattern = ".json","")) %>% 
  ggplot(aes(x=sample_id,y=mapping_ratio))+
  geom_segment(aes(x=sample_id,xend=sample_id,y=0,yend=mapping_ratio),
               color="gray")+
  geom_point(aes(color=mapping_ratio),size=5)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,hjust = 1,vjust=1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.title.y = element_text(size=20,family="KaiTi"),
        axis.text.y = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,100),
                     labels = function(x){paste0(x,"%")})+
  geom_hline(yintercept = 93.88,lty="dashed")+
  labs(x=NULL,y="比对率")+
  scale_color_viridis_c(direction = -1)
dev.off()


## 598175 transitions and 322584 transversions
598175/322584

read_delim("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/snp相对于基因的位置.txt",
           col_names = FALSE,
           delim = " ") %>% 
  mutate(X3=X1/sum(X1)) 

31569/21450

read_delim("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/indel相对于基因的位置.txt",
           col_names = FALSE,
           delim = " ") %>% 
  mutate(X3=X1/sum(X1))


read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/bhagwa_倒位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% range()

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/bhagwa_倒位变异.csv",
         col_names = FALSE) %>% 
  arrange(desc(X4))

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/bhagwa_倒位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% mean()

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/bhagwa_易位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% range()

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/bhagwa_易位变异.csv",
         col_names = FALSE) %>%
  arrange(desc(X4))

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/bhagwa_易位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% mean()


read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/tunisia_倒位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% range()

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/tunisia_倒位变异.csv",
         col_names = FALSE) %>% 
  arrange(desc(X4))

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/tunisia_倒位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% mean()

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/tunisia_易位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% range()

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/tunisia_易位变异.csv",
         col_names = FALSE) %>%
  arrange(desc(X4))

read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/tunisia_易位变异.csv",
         col_names = FALSE) %>% 
  pull(X4) %>% mean()

dat<-read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/indel.count.csv")

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

## 不同长度
70485/(70485+810+657+159)
159/(70485+810+657+159)

## 不同区间
28549/(70485+810+657+159)
1457/(70485+810+657+159) ## cds
(9315+13473+5539)/(70485+810+657+159) ## 上下游
9315+13473+5539

read_tsv("D:/Jupyter/panPome/ys.emapper.annotations",
         col_names = FALSE,
         comment = "#") %>% 
  select(1,10) %>% 
  filter(X10!="-") %>% 
  write_tsv("D:/Jupyter/panPome/pomeGene2Term.tsv",col_names = FALSE)
  
### GO富集分析

library(clusterProfiler)
library(tidyverse)
term2gene<-read_delim("D:/Jupyter/panPome/pomeTerm2Gene.txt",delim = "\t",col_names = FALSE)
term2name<-read_delim("D:/Jupyter/panPome/go.tb",delim = "\t")

core.gene.list<-read_lines("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/indel.impact.list")
enricher(gene = core.gene.list,
         TERM2NAME = term2name,
         TERM2GENE = term2gene,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) -> core.enrich

dotplot(core.enrich)
core.enrich@result %>% colnames()

core.enrich@result %>% select(-geneID) %>% 
  left_join(term2name,by=c("ID"="GO")) %>% 
  filter(pvalue<0.01) %>% pull(level) %>% table()
core.enrich@result %>% select(-geneID) %>% 
  left_join(term2name,by=c("ID"="GO")) %>% 
  filter(pvalue<0.05) %>% 
  filter(level=="MF")

core.enrich@result %>% select(-geneID) %>% 
  left_join(term2name,by=c("ID"="GO")) %>% 
  filter(pvalue<0.01) %>% 
  filter(level=="BP")


read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/merged92.vg.filter.recode.vcf",
         comment = "##") %>% 
  mutate(across(contains("_"),function(x){str_sub(x,1,3)})) -> dat

## 纯合突变
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

save(mut2_freq_list,file = "D:/Jupyter/panPome/Figures/结构变异图形泛基因组/mut2_freq_list.Rdata")
mut2_freq_list %>% length()
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
  filter(nonTi>Ti) -> mut2.dat

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



## 杂合突变
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
  filter(nonTi>Ti) -> mut1.dat

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


bind_rows(mut1.dat,mut2.dat) %>% 
  pull(chrpos) %>% 
  write_lines("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/有利的大片段插入缺失变异位点.txt")


read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/有利变异位点all.anno.variant_function",
         col_names = FALSE) %>% 
  filter(X1=="upstream"|X1=="downstream"|X1=="exonic"|X1=="upstream;downstream") %>% 
  pull(X2) %>% 
  str_split(pattern = ",|;") %>% 
  unlist() %>% 
  str_replace(pattern = "\\(dist=[0-9]+\\)","") %>% 
  unique() %>% 
  paste("mRNA1",sep = ".")-> gene.list

term2gene<-read_delim("D:/Jupyter/panPome/pomeTerm2Gene.txt",delim = "\t",col_names = FALSE)
term2name<-read_delim("D:/Jupyter/panPome/go.tb",delim = "\t")


enricher(gene = gene.list,
         TERM2NAME = term2name,
         TERM2GENE = term2gene,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) -> enrich.dat

enrich.dat@result %>% colnames()
enrich.dat@result %>% 
  select(-geneID)


?wilcox.test()

read_excel("D:/000博士毕业论文/表型数据/01.单果重.xlsx") %>% 
  pull("重量") %>% summary()
read_excel("D:/000博士毕业论文/表型数据/01.单果重.xlsx") %>% 
  pull("重量") %>% sd()

read_excel("D:/000博士毕业论文/表型数据/04.百粒重.xlsx") %>%   
  pull("百粒重") %>% summary()
read_excel("D:/000博士毕业论文/表型数据/04.百粒重.xlsx") %>%   
  pull("百粒重") %>% sd()

read_excel("D:/000博士毕业论文/表型数据/03.果皮厚度mm.xlsx") %>% 
  pull("果皮厚mm") %>% summary()
read_excel("D:/000博士毕业论文/表型数据/03.果皮厚度mm.xlsx") %>%   
  pull("果皮厚mm") %>% sd()

read_excel("D:/000博士毕业论文/表型数据/05.籽粒硬度.xlsx")  %>% 
  pull("kg/cm2") %>% summary()
read_excel("D:/000博士毕业论文/表型数据/05.籽粒硬度.xlsx") %>%   
  pull("kg/cm2") %>% sd()

read_excel("D:/000博士毕业论文/表型数据/09.可滴定酸.xlsx") %>% 
  pull("可滴定酸（%）") %>% summary()
read_excel("D:/000博士毕业论文/表型数据/09.可滴定酸.xlsx") %>%   
  pull("可滴定酸（%）") %>% sd()


read_excel("D:/000博士毕业论文/表型数据/01.单果重.xlsx") %>% 
  group_by(X1) %>% 
  summarise(sf=mean(重量)) -> dat.A 

read_excel("D:/000博士毕业论文/表型数据/04.百粒重.xlsx") %>% 
  group_by(X1) %>% 
  summarise(hgw=mean(百粒重)) -> dat.B

read_excel("D:/000博士毕业论文/表型数据/03.果皮厚度mm.xlsx") %>% 
  group_by(X1) %>% 
  summarise(fst=mean(果皮厚mm)) -> dat.C

read_excel("D:/000博士毕业论文/表型数据/05.籽粒硬度.xlsx") %>% 
  group_by(X1) %>% 
  summarise(gh=mean(`kg/cm2`))-> dat.D

read_excel("D:/000博士毕业论文/表型数据/09.可滴定酸.xlsx") %>% 
  group_by(X1) %>% 
  summarise(ta=mean(`可滴定酸（%）`)) -> dat.E

dat.A %>% 
  left_join(dat.B,by=c("X1"="X1"))%>% 
  left_join(dat.C,by=c("X1"="X1"))%>% 
  left_join(dat.D,by=c("X1"="X1"))%>% 
  left_join(dat.E,by=c("X1"="X1")) %>% 
  ggpairs(data=.,columns = 2:6,
          diag = list(continuous="densityDiag"))+
  theme_bw()+
  theme(panel.grid = element_blank())

dat.A %>% 
  left_join(dat.B,by=c("X1"="X1"))%>% 
  left_join(dat.C,by=c("X1"="X1"))%>% 
  left_join(dat.D,by=c("X1"="X1"))%>% 
  left_join(dat.E,by=c("X1"="X1")) -> pheno.dat
dat.A %>% 
  left_join(dat.B,by=c("X1"="X1"))%>% 
  left_join(dat.C,by=c("X1"="X1"))%>% 
  left_join(dat.D,by=c("X1"="X1"))%>% 
  left_join(dat.E,by=c("X1"="X1")) %>% 
  pull(X1) -> pheno.sample.list

read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/merged92.vg.filter.recode.vcf",
         comment = "##") %>% 
  mutate(across(contains("_"),function(x){str_sub(x,1,3)})) %>% 
  select(-(3:9)) %>% 
  mutate(varSite=paste(`#CHROM`,POS,sep="_")) %>% 
  select(-c("#CHROM","POS")) %>% 
  select(c("varSite",pheno.sample.list)) -> vcf.dat
vcf.dat  
pheno.dat %>% colnames()
pheno.dat
## 单果重
sf_list<-list()
for(i in 1:nrow(vcf.dat)){
  vcf.dat[i,-1] %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("V2") %>% 
    left_join(pheno.dat,by=c("V2"="X1")) %>% 
    #filter_at(vars(contains("_")),any_vars(. != "./.")) %>% 
    filter(V1 != "./.") %>% 
    mutate(V1=case_when(
      V1 == "0/0" ~ "0",
      TRUE ~ "1"
    )) -> temp.dat
  if(length(temp.dat %>% pull(V1) %>% unique()) > 1){
    temp.dat %>% 
      wilcox.test(sf~V1,data = .) %>% 
      .$p.value -> pvalue
  }
    sf_list[[i]]<- data.frame(var_Site=vcf.dat %>% pull(varSite) %>% .[i],
                              pval=pvalue)
}
sf_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  write_csv("D:/000博士毕业论文/表型数据/结构变异GWAS单果重.csv")

## 百粒重
hgw_list<-list()
for(i in 1:nrow(vcf.dat)){
  vcf.dat[i,-1] %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("V2") %>% 
    left_join(pheno.dat,by=c("V2"="X1")) %>% 
    #filter_at(vars(contains("_")),any_vars(. != "./.")) %>% 
    filter(V1 != "./.") %>% 
    mutate(V1=case_when(
      V1 == "0/0" ~ "0",
      TRUE ~ "1"
    )) -> temp.dat
  if(length(temp.dat %>% pull(V1) %>% unique()) > 1){
    temp.dat %>% 
      wilcox.test(hgw~V1,data = .) %>% 
      .$p.value -> pvalue
  }
  hgw_list[[i]]<- data.frame(var_Site=vcf.dat %>% pull(varSite) %>% .[i],
                            pval=pvalue)
}
hgw_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  write_csv("D:/000博士毕业论文/表型数据/结构变异GWAS百粒重.csv")

## 果皮厚度
fst_list<-list()
for(i in 1:nrow(vcf.dat)){
  vcf.dat[i,-1] %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("V2") %>% 
    left_join(pheno.dat,by=c("V2"="X1")) %>% 
    #filter_at(vars(contains("_")),any_vars(. != "./.")) %>% 
    filter(V1 != "./.") %>% 
    mutate(V1=case_when(
      V1 == "0/0" ~ "0",
      TRUE ~ "1"
    )) -> temp.dat
  if(length(temp.dat %>% pull(V1) %>% unique()) > 1){
    temp.dat %>% 
      wilcox.test(fst~V1,data = .) %>% 
      .$p.value -> pvalue
  }
  fst_list[[i]]<- data.frame(var_Site=vcf.dat %>% pull(varSite) %>% .[i],
                             pval=pvalue)
}
fst_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  write_csv("D:/000博士毕业论文/表型数据/结构变异GWAS果皮厚度.csv")

## 籽粒硬度
gh_list<-list()
for(i in 1:nrow(vcf.dat)){
  vcf.dat[i,-1] %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("V2") %>% 
    left_join(pheno.dat,by=c("V2"="X1")) %>% 
    #filter_at(vars(contains("_")),any_vars(. != "./.")) %>% 
    filter(V1 != "./.") %>% 
    mutate(V1=case_when(
      V1 == "0/0" ~ "0",
      TRUE ~ "1"
    )) -> temp.dat
  if(length(temp.dat %>% pull(V1) %>% unique()) > 1){
    temp.dat %>% 
      wilcox.test(gh~V1,data = .) %>% 
      .$p.value -> pvalue
  }
  gh_list[[i]]<- data.frame(var_Site=vcf.dat %>% pull(varSite) %>% .[i],
                             pval=pvalue)
}
gh_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  write_csv("D:/000博士毕业论文/表型数据/结构变异GWAS籽粒硬度.csv")


## 可滴定酸

ta_list<-list()
for(i in 1:nrow(vcf.dat)){
  vcf.dat[i,-1] %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("V2") %>% 
    left_join(pheno.dat,by=c("V2"="X1")) %>% 
    #filter_at(vars(contains("_")),any_vars(. != "./.")) %>% 
    filter(V1 != "./.") %>% 
    mutate(V1=case_when(
      V1 == "0/0" ~ "0",
      TRUE ~ "1"
    )) -> temp.dat
  if(length(temp.dat %>% pull(V1) %>% unique()) > 1){
    temp.dat %>% 
      wilcox.test(ta~V1,data = .) %>% 
      .$p.value -> pvalue
  }
  ta_list[[i]]<- data.frame(var_Site=vcf.dat %>% pull(varSite) %>% .[i],
                             pval=pvalue)
}
ta_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  write_csv("D:/000博士毕业论文/表型数据/结构变异GWAS可滴定酸.csv")


sf_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>%
  filter(padj<0.05) %>% 
  dim()

sf_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>%
  filter(padj<0.05) %>%
  pull(var_Site) %>% 
  write_lines("D:/000博士毕业论文/表型数据/SV_GWAS_sf_varSite.txt")

hgw_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  filter(padj<0.05) %>% 
  dim()

hgw_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  filter(padj<0.05) %>% 
  pull(var_Site) %>% 
  write_lines("D:/000博士毕业论文/表型数据/SV_GWAS_hgw_varSite.txt")

fst_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  filter(padj<0.05) %>% 
  dim()

gh_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>%
  filter(padj<0.05) %>% 
  dim()

ta_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  filter(padj<0.05) %>% 
  dim()

ta_list %>% bind_rows() %>% 
  mutate(padj=p.adjust(pval,method = "bonferroni")) %>% 
  filter(padj<0.05) %>% 
  pull(var_Site) %>% 
  write_lines("D:/000博士毕业论文/表型数据/SV_GWAS_ta_varSite.txt")



read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/单果重相关SV基因/all.anno.variant_function",
         col_names = FALSE) %>% 
  filter(X1=="upstream"|X1=="downstream"|X1=="exonic"|X1=="upstream;downstream") %>% 
  pull(X2) %>% 
  str_split(pattern = ",|;") %>% 
  unlist() %>% 
  str_replace(pattern = "\\(dist=[0-9]+\\)","") %>% 
  unique() %>% 
  paste("mRNA1",sep = ".") %>% 
  as.data.frame() %>% 
  rename("gene_name"=".") %>% 
  left_join(read_tsv("D:/Jupyter/panPome/ys.emapper.annotations",
                     comment = "##") %>% 
              select(1,8),
            by=c("gene_name"="#query")) %>% 
  write_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/单果重相关SV基因/单果重基因功能注释.csv")


read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/all.anno.variant_function",
         col_names = FALSE) %>% 
  filter(X1=="upstream"|X1=="downstream"|X1=="exonic"|X1=="upstream;downstream") %>% 
  pull(X2) %>% 
  str_split(pattern = ",|;") %>% 
  unlist() %>% 
  str_replace(pattern = "\\(dist=[0-9]+\\)","") %>% 
  unique() %>% 
  paste("mRNA1",sep = ".") %>% 
  as.data.frame() %>% 
  rename("gene_name"=".") %>% 
  left_join(read_tsv("D:/Jupyter/panPome/ys.emapper.annotations",
                     comment = "##") %>% 
              select(1,8),
            by=c("gene_name"="#query")) %>% 
  write_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/可滴定酸基因功能注释.csv")

pheno.dat
### 
vcf.dat %>% 
  filter(varSite=="chr6_4216577") %>% 
  select(-varSite) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("V2") %>% 
  left_join(pheno.dat,by=c("V2"="X1")) %>% 
  ggplot(aes(x=V1,y=ta))+
  geom_boxplot()

### 7月
read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/all.anno.variant_function",
         col_names = FALSE) %>% 
  filter(X1=="upstream"|X1=="downstream"|X1=="exonic"|X1=="upstream;downstream") %>% 
  pull(X2) %>% 
  str_split(pattern = ",|;") %>% 
  unlist() %>% 
  str_replace(pattern = "\\(dist=[0-9]+\\)","") %>% 
  unique() %>% 
  paste("mRNA1",sep = ".") %>% 
  as.data.frame() %>% 
  rename("gene_name"=".") %>% 
  left_join(read_tsv("D:/Jupyter/panPome/ys.emapper.annotations",
                     comment = "##") %>% 
              select(1,8),
            by=c("gene_name"="#query")) %>% 
  mutate(gene_name=str_replace(gene_name,pattern = ".mRNA1","")) %>% 
  left_join(read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/00.Jul.25/srv/3drnaseq/X2024.03.19.02.42.47.j986/result/Significant DE genes list and statistics.csv"),
            by=c("gene_name"="target")) %>% 
  na.omit() %>% 
  DT::datatable()

### 8月
read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/all.anno.variant_function",
         col_names = FALSE) %>% 
  filter(X1=="upstream"|X1=="downstream"|X1=="exonic"|X1=="upstream;downstream") %>% 
  pull(X2) %>% 
  str_split(pattern = ",|;") %>% 
  unlist() %>% 
  str_replace(pattern = "\\(dist=[0-9]+\\)","") %>% 
  unique() %>% 
  paste("mRNA1",sep = ".") %>% 
  as.data.frame() %>% 
  rename("gene_name"=".") %>% 
  left_join(read_tsv("D:/Jupyter/panPome/ys.emapper.annotations",
                     comment = "##") %>% 
              select(1,8),
            by=c("gene_name"="#query")) %>% 
  mutate(gene_name=str_replace(gene_name,pattern = ".mRNA1","")) %>% 
  left_join(read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/00.Aug.26/srv/3drnaseq/X2024.03.19.02.58.52.j881/result/Significant DE genes list and statistics.csv"),
            by=c("gene_name"="target")) %>% 
  na.omit() %>% 
  DT::datatable()


### 9月
read_tsv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/all.anno.variant_function",
         col_names = FALSE) %>% 
  filter(X1=="upstream"|X1=="downstream"|X1=="exonic"|X1=="upstream;downstream") %>% 
  pull(X2) %>% 
  str_split(pattern = ",|;") %>% 
  unlist() %>% 
  str_replace(pattern = "\\(dist=[0-9]+\\)","") %>% 
  unique() %>% 
  paste("mRNA1",sep = ".") %>% 
  as.data.frame() %>% 
  rename("gene_name"=".") %>% 
  left_join(read_tsv("D:/Jupyter/panPome/ys.emapper.annotations",
                     comment = "##") %>% 
              select(1,8),
            by=c("gene_name"="#query")) %>% 
  mutate(gene_name=str_replace(gene_name,pattern = ".mRNA1","")) %>% 
  left_join(read_csv("D:/Jupyter/panPome/Figures/结构变异图形泛基因组/可滴定酸相关SV基因/00.Sep.24/srv/3drnaseq/X2024.03.19.08.40.57.j983/result/Significant DE genes list and statistics.csv"),
            by=c("gene_name"="target")) %>% 
  DT::datatable()
