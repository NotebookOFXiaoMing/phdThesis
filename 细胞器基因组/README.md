## 叶绿体基因组组装

```
get_organelle_from_reads.py -1 ../../../upload2ncbi/YS_clean_1.fq.gz -2 ../../../upload2ncbi/YS_clean_2.fq.gz -t 16 -o getorganelle.output -F embplant_pt -R 10
get_organelle_from_assembly.py -F embplant_pt -g ys.cp.graph.gfa -o ysCP -t 2
#这里如果报错 BLAST Database error: Error: Not a valid version 4 database
#使用命令 get_organelle_config.py --add embplant_pt 重新构建一些blast库
```

### 叶绿体基因组组装和注释

http://47.96.249.172:16019/analyzer/annotate

http://47.96.249.172:16019/analyzer/view
170654022278923

## 叶绿体基因组核苷酸多样性

```
time ~/anaconda3/envs/OrganelleGenome/bin/halign -t 2 downLoad.fa
# downLoad.fa.aligned

i = 1
step = []

while i < 161224:
    step.append(i)
    i += 1000

for i in step:
    fw = open("D:/Jupyter/panPome/细胞器基因组/叶绿体基因组滑窗/" + "window_" + str(i) + ".fa",'w' )
    for rec in SeqIO.parse("D:/Jupyter/panPome/细胞器基因组/downLoad.fa.aligned",'fasta'):
        fw.write(">%s\n%s\n"%(rec.id,str(rec.seq[i-1:i+2000])))
    fw.close()

library(ape)
library(tidyverse)
myfun<-function(x){
  tmp.dna<-read.dna(x,format="fasta") 
  
  df<-data.frame(sample_name = str_extract(x,pattern = "[0-9]+") %>% as.numeric(),
                 nuc_div=pegas::nuc.div(tmp.dna))
  return(df)
}
list.files("D:/Jupyter/panPome/细胞器基因组/叶绿体基因组滑窗/",
           pattern = "*.fa",
           full.names = TRUE) %>% 
  map(myfun) %>% 
  bind_rows() %>% 
  write_tsv("D:/Jupyter/panPome/细胞器基因组/cp.pi")
```


### 叶绿体基因组circos图

http://47.96.249.172:16085/cpgview/drawmap

170884098974144

### 线粒体基因组组装

```
### 矫正ont数据
conda activate OrganelleGenome
Ratatosk correct -v -c 24 -s YS_clean_1.fq.gz YS_clean_2.fq.gz -l ys.nanopore.fq.gz -o ratatosk_correct_long_reads
~/biotools/PMAT-1.2.2/bin/PMAT autoMito -i ../../upload2ncbi/ratatosk_correct_long_reads.fastq -o ys.pome.mito -st hifi -g 340m -cpu 24 -m
```
