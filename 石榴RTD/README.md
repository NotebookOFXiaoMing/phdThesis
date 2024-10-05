### 参考链接

https://github.com/wyguo/genotype_specific_RTD/tree/main

https://github.com/anonconda/RTDmaker/tree/main

https://github.com/maxecoulter/BaRT-2/tree/master

https://github.com/cropgeeks/barleyPantranscriptome/tree/main

### 参考基因组构建索引

```
time STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ref.index --genomeFastaFiles ref/pome.ref.nonref.fna --outFileNamePrefix pome --limitGenomeGenerateRAM 240000000000
```

### star第一轮比对

这里如果是压缩文件会报错 FATAL ERROR in input reads: unknown file format: the read ID should start with @ or >

https://github.com/alexdobin/STAR/issues/379 需要先把文件解压

```
snakemake -s starPassOne.smk --cores 128 -p
cat 01.alignment/*/*SJ.out.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > first.pass.SJ.out.tab
```

### 第二轮参考基因组构建索引

```
mkdir ref.index.second.pass
time STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ref.index.second.pass --genomeFastaFiles ref/pome.ref.nonref.fna --outFileNamePrefix pome --limitGenomeGenerateRAM 240000000000 --sjdbFileChrStartEnd first.pass.SJ.out.tab
```

### 第二轮比对

```
snakemake -s starPassTwo.smk --cores 128 -p
```

### stringtie组装转录本

```
snakemake -s stringtieSTAR.smk --cores 128 -p
snakemake -s cpGTFSJ.smk --cores 128 -p
```

### 运行RTDmaker

```
## 运行了8个多小时
python ~/biotools/RTDmaker-main/RTDmaker.py ShortReads --assemblies 04.assemblies --references ref.gtf --SJ-data 04.SJdata --genome ref/pome.ref.nonref.fa --fastq clean.fastq --outpath pomeRTD --outname pome --prefix pomeRTD --SJ-reads 2 1 --tpm 0.1 1 --fragment-len 0.7 --antisense-len 0.5 --add intronic --keep intermediary --ram 8
```
## 更改tpm 参数 设置为0.1 1 转录本15万多 改成 1 1 最小表达量是1 在 一个样本中
## 备份之前的结果
```
cp -r pomeRTD pomeRTD.backup/
```


## clean.fastq 的数据为了节省存储空间给删掉了，新建一个clean.fastq文件夹
## 把一个bioproject的转录组数据复制到这个文件夹下

```
mkdir clean.fastq
cp ../../rnaseq/PRJNA901890/*.gz clean.fastq/
```

## 运行

```
python ~/biotools/RTDmaker-main/RTDmaker.py ShortReads --assemblies 04.assemblies --references ref.gtf --SJ-data 04.SJdata --genome ref/pome.ref.nonref.fa --fastq clean.fastq --outpath pomeRTD --outname pome --prefix pomeRTD --SJ-reads 2 1 --tpm 1 1 --fragment-len 0.7 --antisense-len 0.5 --add intronic --keep intermediary --ram 8
```

## 提示

```
 File "/data/myan/raw_data/pome/sour.pome/20231015.reanalysis/24.pomeRTD/pomeRTD/pome_RTDmaker_output/intermediary/pome_assemblies1_non_redundant.gtf" already exist (924.6MB)
Wed Aug 14 22:13:49 2024 Keeping current file
```
## 可能就不会根据表达量过滤了
## pome_assemblies1_non_redundant.gtf 把这个文件删了试试
## 这样应该是不行

## RTDmaker 的github主页说检测到intermediary 就不会从头运行，测试用一个bioproject的数据能不能行

## 用 tpm 0.1 1参数构建的参考转录本去计算表达量，获取表达量矩阵以后再进行过滤

## 跳转到下面salmon步骤


### 运行transuite

```
python ~/biotools/TranSuite-main/transuite.py Auto --gtf pomeRTD/pome_RTDmaker_output/pome.gtf --fasta pomeRTD/pome_RTDmaker_output/pome.fa --outpath transuite.output --outname pome
```

### 获取具有蛋白编码能力的基因的最长蛋白
```
python getLongestPEPfromPomeRTD.py transuite.output/pome_TranSuite_output/pome_transfeat/pome_transfeat.csv pomeRTDcoding.pep
```
### eggnog 蛋白功能注释

```
conda activate eggnogmapper
export EGGNOG_DATA_DIR=/home/myan/my_data/database/eggnog/
emapper.py -i ../pomeRTDcoding.pep -o pomeRTD --cpu 32 -m diamond --excel
```

### PH基因blast

```
makeblastdb -in ../pomeRTDcoding.pep -dbtype prot -title pomeRTD -parse_seqids -out pomeRTD
```


### salmon 转录本量化

```
salmon index -t ../pomeRTD/pome_RTDmaker_output/pome_padded.fa -k 31 -i pomeRTD_index -p 8

## Jul25
snakemake -s salmon_quant.smk --cores 128 -p
```

### 有多个文件夹，每个文件夹下有对应的转录组数据,这种
### 怎么一次性用snakemake去匹配 暂时没有搞明白
### 我这里的处理方式是 固定文件夹去匹配，然后手动修改文件夹

## R tximport 包提取表达量

## 参考 https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Import_transcript-level_estimates

```
library(tidyverse)
library(magrittr)
library(tximport)

rtracklayer::import("../pomeRTD.backup/pomeRTD/pome_RTDmaker_output//pome_padded.gtf") -> gtf.dat
gtf.dat%>%as.data.frame()%>%filter(type=="transcript")%>%select(transcript_id,gene_id)%>%set_colnames(c("TXNAME","GENEID")) -> tx2gene

## tx2gene 两列 第一列 TXNAME 第二列 GENEID

## 转录本水平的表达量
file<-"salmon.quant.output/PRJNA360679/SRR5279386/quant.sf"
names(file)<-"SRR5279386"
txi.salmon<-tximport(file, type = "salmon",txOut = TRUE)
txi.salmon$abundance%>%head()

## 基因水平的表达量
txi.salmon<-tximport(file, type = "salmon", tx2gene = tx2gene)

files<-list.files("salmon.quant.output",pattern = "quant.sf",recursive = TRUE,full.names = TRUE)
names(files)<-str_extract(files,pattern = "SRR[0-9]+")

txi.salmon.transcript<-tximport(files, type = "salmon",txOut = TRUE)
txi.salmon.gene<-tximport(files, type = "salmon",tx2gene = tx2gene)

save(txi.salmon.gene,txi.salmon.transcript,file = "txi.salmon.Rdata")

txi.salmon.transcript$abundance %>% as.data.frame() -> transcript.tpm
## 将tpm值大于0.5定义为表达
transcript.tpm[transcript.tpm > 0.5] <- 1
transcript.tpm[transcript.tpm <= 0.5] <- 0
transcript.tpm[rowSums(transcript.tpm) >= 1,]%>%dim()
## 总共是168112个转录本 0.5过滤后 146816
## RTDmaker 默认参数是 1 1 最少在1个样本中的表达量大于1

transcript.tpm[rowSums(transcript.tpm) >= 1,]%>%rownames_to_column("TXNAME")%>%left_join(tx2gene,by=c("TXNAME"="TXNAME"))%>%pull(GENEID)%>%unique()%>%length()
## 基因剩 54924

keep.transcript<-transcript.tpm[rowSums(transcript.tpm) >= 1,]%>%rownames_to_column("TXNAME")%>%select(TXNAME)
gtf.dat%>%as.data.frame()%>%inner_join(keep.transcript,by=c("transcript_id"="TXNAME"))%>%rtracklayer::export("pome_padded_filtered.gtf")

rtracklayer::import("../pomeRTD.backup/pomeRTD/pome_RTDmaker_output//pome.gtf")%>%as.data.frame()%>%inner_join(keep.transcript,by=c("transcript_id"="TXNAME"))%>%rtr
    acklayer::export("pome.filtered.gtf")
```

## 提取转录本序列

```
gffread pome_padded_filtered.gtf -g ../ref/pome.ref.nonref.fa -w pomeRTD_padded_filtered.trans.fa
gffread pome.filtered.gtf -g ../ref/pome.ref.nonref.fa -w pomeRTD.filtered.trans.fa

## 文件所在路径 05.salmon.quant
```

## 重新运行 transuite

```
conda activate rnaseq
time python ~/biotools/TranSuite-main/transuite.py Auto --gtf pome.filtered.gtf --fasta pomeRTD.filtered.trans.fa --outpath transuite.output --outname pome
## 15min
```

## 平均外显子数量

```
rtracklayer::import("pome.filtered.gtf")%>%as.data.frame()%>%filter(type=="exon")%>%group_by(transcript_id,gene_id)%>%summarise(exon_num=n())%>%filter(exon_num>1)%>%pull(exon_num)%>%summary()
```

## 盐胁迫转录组数据处理

root 1 2 鉴定到5个差异可变剪切的基因，使用eggmapper做功能注释

```
seqkit grep -r -f gene.list ../05.salmon.quant/transuite.output/pome_TranSuite_output/pome_transfeat/pome_transfeat_pep.fasta -o candidate.peps
conda activate eggnogmapper
export EGGNOG_DATA_DIR=/home/myan/my_data/database/eggnog/
emapper.py -i candidate.peps -o candidate --cpu 32 -m diamond --excel
```

## 核心转录本 总转录本数量 拟合

Pan-GP https://pangp.zhaopage.com/download.html
```
load("txi.salmon.Rdata")
library(tidyverse)
txi.salmon.gene$abundance%>%as.data.frame() -> tpm.gene
tpm.gene[tpm.gene > 0.5] <- 1
tpm.gene[tpm.gene <= 0.5] <- 0
tpm.gene[rowSums(tpm.gene) >= 1,]%>%write_delim(col_names = FALSE,file = "tpm.gene.01.matrix",delim = "")

```

## 数据路径

```
## 20231015.reanalysis/24.pomeRTD/05.salmon.quant
## 核心基因 转录本
load("txi.salmon.Rdata")

rtracklayer::import("pome.filtered.gtf")%>%as.data.frame()%>%filter(type=="transcript")%>%select(transcript_id) -> contained.trans
rtracklayer::import("pome.filtered.gtf")%>%as.data.frame()%>%filter(type=="transcript")%>%select(gene_id)%>%distinct()%>%dim()

rtracklayer::import("pome.filtered.gtf")%>%as.data.frame()%>%filter(type=="transcript")%>%select(gene_id)%>%distinct() -> contained.genes
txi.salmon.transcript$abundance%>%as.data.frame()%>%rownames_to_column("trans_id")%>%inner_join(contained.trans,by=c("trans_id"="transcript_id"))%>%write_tsv("pomeR
    TD.transcript.274.TPM.txt")

txi.salmon.gene$abundance%>%as.data.frame()%>%rownames_to_column("gene_id")%>%inner_join(contained.genes,by=c("gene_id"="gene_id"))%>%write_tsv("pomeRTD.gene.274.TP
    M.txt")


grep ",Coding," transuite.output/pome_TranSuite_output/pome_transfeat/pome_transfeat.csv | awk -F, '{print $1}' | sort -u | uniq | wc -l

grep ",Coding," transuite.output/pome_TranSuite_output/pome_transfeat/pome_transfeat.csv | awk -F, '{print $1}' | sort -u > pomeRTD.coding.gene.ids


## 基因共表达网络

data.table::fread("phdthesis/chapter6/data/pomeRTD.gene.274.TPM.txt.gz") %>% 
  inner_join(read_tsv("phdthesis/chapter6/data/pomeRTD.coding.gene.ids",col_names = FALSE),
             by=c("gene_id"="X1")) -> pome.RTD.coding.gene.exp


myfun<-function(x){
  length(which(x>0.5))
}

pome.RTD.coding.gene.exp %>% 
  #column_to_rownames("gene_id") %>% 
  dplyr::rowwise() %>% 
  mutate(group=myfun(c_across(starts_with("SRR")))) %>%
  filter(group>=55) %>% 
  dplyr::select(-group) -> pome.RTD.coding.gene.exp.filter

pome.RTD.coding.gene.exp.filter %>% column_to_rownames("gene_id") %>%
  .[1:5,1:5]

pcc <- cor(pome.RTD.coding.gene.exp.filter %>% 
             column_to_rownames("gene_id") %>% 
             t(), method = "pearson")
pcc[1:5,1:5]
write.table(pcc, file = "phdthesis/chapter6/data/pomeRTD.coding.gene.pcc.matrix", quote = F)
write.table(pcc, file = "phdthesis/chapter6/data/pomeRTD.coding.gene.pcc.matrix", quote = F)

python ~/my_data/raw_data/practice/GCN/computing_HRR_matrix_TOP420.py -p pomeRTD.coding.gene.pcc.matrix -o HRR.matrix -t 8
time python ~/my_data/raw_data/practice/GCN/top1_co_occurrence_matrix_version2_TOP420_removing_ties.py -p HRR.matrix -c non_agg_filtered_net_Cyto.csv -e non_agg_filtered_net_EGAD.csv
## 121 min

time python ~/my_data/raw_data/practice/GCN/top1_co_occurrence_matrix_version2_TOP420_keeping_ties.py -p HRR.matrix -c non_agg_full_net_Cyto.csv -e non_agg_full_net_EGAD.csv


## emapper
python ../selectMostLengthPEP.py ../transuite.output/pome_TranSuite_output/pome_transfeat/pome_transfeat_coding_pep.fa pomeRTD.coding.pep
conda activate eggnogmapper
export EGGNOG_DATA_DIR=/home/myan/my_data/database/eggnog/
emapper.py -i pomeRTD.coding.pep -o pomeRTD.coding --cpu 24 -m diamond --excel


read_tsv("../emapper.output//pomeRTD.coding.emapper.annotations",comment = "##")%>%select(`#query`,`GOs`)%>%filter(GOs!="-")%>%nest(.by="#query")%>%mutate(new_col=m
    ap(data,myfun))%>%unnest(new_col)%>%select(-data)%>%pivot_longer(!`#query`)%>%select(1,3)%>%magrittr::set_colnames(c("gene_id","GO"))%>%write_tsv("pomeRTD.coding.ge
    nes.GO.table")

net.dat<-data.table::fread("non_agg_filtered_net_EGAD.csv")%>%column_to_rownames("V1")
net.dat[net.dat>1] = 1

read_tsv("../emapper.output//pomeRTD.coding.emapper.annotations",comment = "##")%>%select(`#query`,`GOs`)%>%filter(GOs!="-")%>%nest(.by="#query")%>%mutate(new_col=m
    ap(data,myfun))%>%unnest(new_col)%>%select(-data)%>%pivot_longer(!`#query`)%>%select(1,3)%>%magrittr::set_colnames(c("gene_id","GO")) -> go.table

gene_list<-go.table%>%pull(gene_id)%>%unique()
go_term<-go.table%>%pull(GO)%>%unique()

anno <- make_annotations(go.table,gene_list,go_term)
save(anno,file = "anno.Rdata")

GO_groups_voted <- run_GBA(network = net.dat%>%as.matrix(), anno)


read_tsv("go.table.txt")%>%mutate(a=1)%>%distinct()%>%pivot_wider(names_from = "GO",values_from = "a")%>%naniar::replace_na_with(0)%>%column_to_rownames("gene_id") 
    -> anno

GO_groups_voted <- run_GBA(network = net.dat%>%as.matrix(), anno)

GO_groups_voted[[1]]%>%as.data.frame()%>%rownames_to_column()%>%write_tsv("go.auc.txt")

read_tsv("../emapper.output/pomeRTD.coding.emapper.annotations",comment = "##")%>%select(`#query`,KEGG_ko)%>%filter(KEGG_ko!="-")%>%mutate(value=1)%>%pivot_wider(na
    mes_from = "KEGG_ko",values_from = "value")%>%naniar::replace_na_with(0)%>%column_to_rownames("#query") -> kegg.anno
kegg_groups_voted <- run_GBA(network = net.dat%>%as.matrix(), kegg.anno)

kegg_groups_voted[[1]]%>%as.data.frame()%>%rownames_to_column()%>%write_tsv("kegg.auc.txt")

data.table::fread("non_agg_filtered_net_Cyto.csv")%>%magrittr::set_colnames(c("from","to","weight")) -> net.dat
nodes<-c(net.dat %>% pull(from),net.dat %>% pull(to)) %>% unique()
net_pc<-graph_from_data_frame(
  d=net.dat,vertices=nodes,
  directed=TRUE)
mean(degree(net_pc))


makeblastdb -in pomeRTD.coding.pep -dbtype prot -parse_seqids -out pomeRTD.coding.pep.db
seqkit grep -r -p "ys008G000883" ../../../11.orthofinder/all.genomes/ys.pep.fa -o ys008G000883.pep

net.dat%>%filter(from=="pomeRTD_chr8G014220"|to=="pomeRTD_chr8G014220")%>%write_tsv("pomeRTD_chr8G014220.gene.net.data")
```
