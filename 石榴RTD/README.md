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
library(ridyverse)
library(magrittr)
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

## 总共是168112个转录本 0.5过滤后 146818
## RTDmaker 默认参数是 1 1 最少在1个样本中的表达量大于1
```
