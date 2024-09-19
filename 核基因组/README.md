### 三代测序数据统计

```
NanoPlot -t 24 --N50 --dpi 300 --fastq ys.nanopore.fq.gz --title ys -o ys_NanoPlot.output
```

### 桃金娘目的基因组Helixer注释

```
sbatch run_helixer_Pome.slurm
```

### 基因组gap位置

```
python findGapPosition.py ys.final.fna gap.pos
python ../06.repeatModulerRepeatMasker/findGapPosition.py ys.genome.fasta gap.pos.before
```

### 按步长计算GC含量

```
python calculateGC.py ../../08.proteinCodingGenes/ys.final.masked.fna ../02.mapping.ratio/illumina.coverage gc.content
```

### 蛋白编码基因busco
```
conda activate busco5.4.6
busco -i ../ys.pep.fa -l ~/my_data/database/embryophyta_odb10 -o output -m prot -c 24 --offline
```

### 蛋白编码基因比对NR库
```
diamond blastp -d ~/my_data/database/diamond.nr/nr -q ../ys.pep.fa -o output.txt
```

### TE密度

```
grep -v "Simple_repeat" genome.nextpolish.upperbase.fasta.out | awk '{print $5"\t"$6"\t"$7}' > TE.100k
cat ../../08.proteinCodingGenes/ys.final.masked.fna.fai | awk '{print $1"\t"$2}' > chr.len
bedtools makewindows -g chr.len -w 100000 -s 50000 > chr.100k.window
bedtools coverage -a chr.100k.window -b TE.100k -counts > TE.100k.counts
```

### TRF 1-6 bp 串联重复序列

```
trf ../08.proteinCodingGenes/ys.final.masked.fna 2 7 7 80 10 5 500 -h -d -l 6 -ngs > ys.pome.trf

seqkit grep -r ../08.proteinCodingGenes/ys.final.masked.fna -p chr1 -o chr1.fna
trf chr1.fna 2 7 7 80 10 5 500 -h -d -l 6 -ngs > chr1.trf
python fix_trf_output.py chr1.trf chr1.trf.fix # 多条序列会报错

snakemake -s TRF.smk --cores 128 -p

```

### trash 结果统计

```
library(tidyverse)
list.files("./",pattern = "Summary.of.repetitive.regions.chr*",recursive = TRUE)%>%map(read.csv,row.name=1)%>%bind_rows()%>%write_csv("trash.output.summary.csv")
list.files("./",pattern = "Summary.of.repetitive.regions.sca*",recursive = TRUE)%>%map(read.csv,row.name=1)%>%bind_rows()%>%write_csv("nonchr.trash.output.summary.csv")

bedtools coverage -a ../06.repeatModulerRepeatMaskerafterNextpolish/ys.repeatmasker/chr.100k.window -b trash.TR.bed -counts > trash.TR.counts
```

### EDTA

```
time EDTA.pl --genome ../29.centromics/ys.Chr.fna --species others --step all --cds ../08.proteinCodingGenes/05.evm/ys/ys.cds.fa --sensitive 1 --anno 1 -t 96

python ~/biotools/quarTeT-main/quartet.py CentroMiner -i ../../../29.centromics/ys.Chr.fna --TE ../../../31.edta/ys.Chr.fna.mod.EDTA.anno/ys.Chr.fna.mod.EDTA.TEanno.gff3
python ~/biotools/quarTeT-main/quartet.py CentroMiner -i ../../../../29.centromics/ys.Chr.fna --TE ../../../../31.edta/ys.Chr.fna.mod.EDTA.anno/ys.Chr.fna.mod.EDTA.TEanno.gff3 -n 30 -m 2000 -t 24
grep -v "@" quarTeT.best.candidate | grep -v "#" > centrometer.candidate
```


### 可视化stainedGlass结果

```
library(data.table)
library(tidyverse)
dat<-fread("results//pome.2000.10000.bed.gz")
dat%>%filter(`#query_name`=="chr1"&reference_name=="chr1")%>%filter(query_start>=31000000&query_start<=33000000)%>%filter(reference_start>=31000000&reference_start<=33000000)%>%write_tsv("chr1.candidate.cnetrometer.stainedGlass")

filterStainedGlass<-function(chr.name,chr.start,chr.end,output.prefix){
  dat%>%
    filter(`#query_name`==chr.name&reference_name==chr.name)%>%
    filter(query_start>=chr.start&query_start<=chr.end)%>%
    filter(reference_start>=chr.start&reference_start<=chr.end)%>%
    write_tsv(paste0(output.prefix,".candidate.cnetrometer.stainedGlass"))
}
filterStainedGlass(chr.name = "chr1",chr.start = 29200000,chr.end = 29300000,output.prefix = "chr1")
filterStainedGlass(chr.name = "chr2",chr.start = 29200000,chr.end = 29300000,output.prefix = "chr2")
filterStainedGlass(chr.name = "chr3",chr.start = 21000000,chr.end = 23000000,output.prefix = "chr3")
filterStainedGlass(chr.name = "chr4",chr.start = 19600000,chr.end = 20000000,output.prefix = "chr4")
filterStainedGlass(chr.name = "chr5",chr.start = 18600000,chr.end = 19500000,output.prefix = "chr5")
filterStainedGlass(chr.name = "chr6",chr.start = 11500000,chr.end = 12500000,output.prefix = "chr6")
filterStainedGlass(chr.name = "chr7",chr.start = 12600000,chr.end = 15600000,output.prefix = "chr7")
filterStainedGlass(chr.name = "chr8",chr.start = 15400000,chr.end = 17300000,output.prefix = "chr8")
filterStainedGlass(chr.name = "chr8",chr.start = 15400000,chr.end = 17300000,output.prefix = "chr8")
```

### 可能得着丝粒区域和串联重复区取交集

```

cat ../chr1/all.repeats.from.chr1.fna.csv | grep -v "start" | awk -v FS="," '{print "chr1\t"$1"\t"$2"\t"$3"\t"$4}' > chr1.repeats.bed
bedtools intersect -a chr1.bed -b chr1.repeats.bed -wa -wb > chr1.candidate.centrometer.repeat

cat ../chr2/all.repeats.from.chr2.fna.csv | grep -v "start" | awk -v FS="," '{print "chr2\t"$1"\t"$2"\t"$3"\t"$4}' > chr2.repeats.bed
bedtools intersect -a chr2.bed -b chr2.repeats.bed -wa -wb > chr2.candidate.centrometer.repeat

cat ../chr3/all.repeats.from.chr3.fna.csv | grep -v "start" | awk -v FS="," '{print "chr3\t"$1"\t"$2"\t"$3"\t"$4}' > chr3.repeats.bed
bedtools intersect -a chr3.bed -b chr3.repeats.bed -wa -wb > chr3.candidate.centrometer.repeat

cat ../chr4/all.repeats.from.chr4.fna.csv | grep -v "start" | awk -v FS="," '{print "chr4\t"$1"\t"$2"\t"$3"\t"$4}' > chr4.repeats.bed
bedtools intersect -a chr4.bed -b chr4.repeats.bed -wa -wb > chr4.candidate.centrometer.repeatv
cat ../chr5/all.repeats.from.chr5.fna.csv | grep -v "start" | awk -v FS="," '{print "chr5\t"$1"\t"$2"\t"$3"\t"$4}' > chr5.repeats.bed
bedtools intersect -a chr5.bed -b chr5.repeats.bed -wa -wb > chr5.candidate.centrometer.repeat

cat ../chr6/all.repeats.from.chr6.fna.csv | grep -v "start" | awk -v FS="," '{print "chr6\t"$1"\t"$2"\t"$3"\t"$4}' > chr6.repeats.bed
bedtools intersect -a chr6.bed -b chr6.repeats.bed -wa -wb > chr6.candidate.centrometer.repeat

cat ../chr7/all.repeats.from.chr7.fna.csv | grep -v "start" | awk -v FS="," '{print "chr7\t"$1"\t"$2"\t"$3"\t"$4}' > chr7.repeats.bed
bedtools intersect -a chr7.bed -b chr7.repeats.bed -wa -wb > chr7.candidate.centrometer.repeat

cat ../chr8/all.repeats.from.chr8.fna.csv | grep -v "start" | awk -v FS="," '{print "chr8\t"$1"\t"$2"\t"$3"\t"$4}' > chr8.repeats.bed
bedtools intersect -a chr8.bed -b chr8.repeats.bed -wa -wb > chr8.candidate.centrometer.repeat



```

### 着丝粒区域的TE
```
## 路径 20231015.reanalysis/22.trash/01.bedtools.intersect

cat ../../06.repeatModulerRepeatMaskerafterNextpolish/ys.repeatmasker/genome.nextpolish.upperbase.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' | grep "chr" > chr.TE.bed
bedtools intersect -a chr.candidate.centrometer.bed -b chr.TE.bed -wa -wb > candidate.centrometer.TE.type

bedtools intersect -a chr.non.candidate.centrometer.bed -b ../../06.repeatModulerRepeatMaskerafterNextpolish/TE.count.100k.window -wa -wb > non.candidate.centrometer.TE
bedtools intersect -a chr.candidate.centrometer.bed -b ../../06.repeatModulerRepeatMaskerafterNextpolish/TE.count.100k.window -wa -wb > candidate.centrometer.TE

cp chr.non.candidate.centrometer.bed chr.non.candidate.centrometer.bed01 ## 看看着丝粒上下游1M
vim chr.non.candidate.centrometer.bed01

bedtools intersect -a chr.non.candidate.centrometer.bed01 -b ../../06.repeatModulerRepeatMaskerafterNextpolish/TE.count.100k.window -wa -wb > non.candidate.centrometer.TE01

## 基因数量
bedtools intersect -a chr.candidate.centrometer.bed -b ../../08.proteinCodingGenes/05.evm/ys/gene.count.100k.window -wa -wb > candidate.centrometer.Gene
bedtools intersect -a chr.non.candidate.centrometer.bed01 -b ../../08.proteinCodingGenes/05.evm/ys/gene.count.100k.window -wa -wb > non.candidate.centrometer.Gene01

grep "gene" ../../08.proteinCodingGenes/05.evm/ys/ys.rename.gff3 | grep "chr" | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk 'gsub("ID=","")' > gene.bed
bedtools intersect -a chr.candidate.centrometer.bed -b gene.bed -wa -wb | wc -l #140

bedtools intersect -a chr.candidate.centrometer.bed -b gene.bed -wa -wb | awk '{print $7}' > chr.candidate.centrometer.140genes

## GO富集分析
library(clusterProfiler)

GOinfo<-read_tsv("phdthesis/chapter4/data/go.tb")
GOinfo

GOannotation <- read_tsv("phdthesis/chapter4/data/ref.output/GOannotation.tsv") %>% 
  filter(GO!="-")
GOannotation

centrometer.gene.list<-read_lines("phdthesis/chapter2/data/chr.candidate.centrometer.140genes")%>% 
  paste0(".mRNA1")
centrometer.gene.list 

enricher(centrometer.gene.list,
         TERM2GENE = GOannotation %>% 
           filter(level=="MF") %>% 
           select(2,1),
         TERM2NAME = GOinfo[1:2],
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05)

enricher(centrometer.gene.list,
         TERM2GENE = GOannotation %>% 
           filter(level=="BP") %>% 
           select(2,1),
         TERM2NAME = GOinfo[1:2],
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) %>% 
  dotplot()+
  theme(panel.grid = element_blank())+
  labs(title = "Biological Processes") -> p10

enricher(centrometer.gene.list,
         TERM2GENE = GOannotation %>% 
           filter(level=="CC") %>% 
           select(2,1),
         TERM2NAME = GOinfo[1:2],
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) %>% 
  dotplot()+
  theme(panel.grid = element_blank())+
  labs(title = "Cell Components") -> p11

```


## 泛基因家族
```
library(tidyverse)
list.files("emapper",pattern = "*annotations$",full.names = TRUE)%>%map(read_tsv,comment="##")%>%bind_rows()%>%write_tsv("emapper.all.six.pome.tsv")

## 基因长度
library(rtracklayer)
myGeneLen<-function(x){
  import(x) %>% 
    as.data.frame() %>% 
    filter(type=="mRNA") %>% 
    select(width,ID) %>% 
    return()
}

list.files("all.genomes",pattern = "*.gff3",full.names = TRUE)%>%map(myGeneLen)%>%bind_rows()%>%write_tsv("allGenesLength.txt")

myExonNum<-function(x){
  import(x) %>% 
    as.data.frame() %>% 
    filter(type=="exon") %>% 
    mutate(ID=str_replace(ID,pattern = "exon[0-9]+","mRNA1")) %>% 
    select(ID) %>%
    group_by(ID) %>% 
    summarise(count=n()) %>% 
    return()
}

list.files("all.genomes",pattern = "*.gff3",full.names = TRUE)%>%map(myExonNum)%>%bind_rows()%>%write_tsv("allGenesExonNum.txt")

read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/01.coreGene.txt",col_names=FALSE)  %>% 
  mutate(group="Core") %>% 
  bind_rows(
    read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/02.dispensableGene.txt",col_names=FALSE)  %>% 
      mutate(group="Dispensable") 
  ) %>% 
  bind_rows(
    read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/03.privateGene.txt",col_names=FALSE)  %>% 
      mutate(group="Private") 
  ) %>% 
  inner_join(read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/allGenesExonNum.txt",col_names=TRUE),
             by=c("X1"="ID")) %>% 
  ggplot(aes(x=group,y=log2(count)))+
  geom_boxplot()


read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/01.coreGene.txt",col_names=FALSE)  %>% 
  mutate(group="Core") %>% 
  bind_rows(
    read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/02.dispensableGene.txt",col_names=FALSE)  %>% 
      mutate(group="Dispensable") 
  ) %>% 
  bind_rows(
    read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/03.privateGene.txt",col_names=FALSE)  %>% 
      mutate(group="Private") 
  ) %>% 
  inner_join(read_tsv("D:/Jupyter/panPome/Figures/第一章_核基因组/泛基因家族/allGenesLength.txt",col_names=TRUE),
             by=c("X1"="ID")) %>% 
  ggplot(aes(x=group,y=log2(width)))+
  geom_boxplot()

## 提取单拷贝cds 计算kaks 和核苷酸多样性
python extractCDSfromtotal.py CoreSingleCopyFamilyIDs.txt
conda activate doubletrouble
Rscript calculateKaKsUsingCDS.R core.kaks

~/anaconda3/envs/syri/bin/snakemake -s calculateKaKsUsingCDS.smk --cores 128 -p -k
python ../01.Core.CDS/extractCDSfromtotal.py DispensableSingleCopyFamilyIDs.txt

~/anaconda3/envs/syri/bin/snakemake -s calculateKaKsUsingCDS.smk --cores 128 -p -k

cat 01.Core.CDS.kaks/*.kaks | grep -v "dup" > core.kaks
cat 02.Dispensable.CDS.kaks/*.kaks | grep -v "dup" > dispensable.kaks


snakemake -s mafft.smk --cores 128 -p -k

## 计算核苷酸多样性
library(tidyverse)
library(pegas)
myfun<-function(x){return(nuc.div(read.dna(x,format = "fasta")))}
list.files("01.Core.CDS.aln/",pattern = "*.fasta",full.names = TRUE)%>%map(myfun)%>%unlist()%>%write_lines("core.nucldiv")
list.files("02.Dispensable.CDS.aln/",pattern = "*.fasta",full.names = TRUE)%>%map(myfun)%>%unlist()%>%write_lines("dispensable.nucldiv")

## GO富集

cat emapper/*.annotations | grep -v "#" | awk -v FS="\t" '{print $10"\t"$1}' | grep -v "-" > Term2Gene.temp
python getTerm2Gene.py Term2Gene.temp Term2Gene.txt
```

## 泛基因家族中的基因表达

```
## 只用 ys 参考基因组中的基因
## 路径 pan.raw.fq/11.SNP.smallInDel/01.YS.Genome/00.rnaseq/
myfun<-function(x){read_tsv(x)%>%mutate(sampleid=str_extract(x,pattern="SRR[0-9]+"))}

list.files(".",pattern = "*_abund.tsv",recursive = TRUE,full.names = TRUE)%>%map(myfun)%>%bind_rows()%>%select(`Gene ID`,sampleid,TPM) -> exp.dat
exp.dat%>%pivot_wider(names_from = "sampleid",values_from = "TPM")%>%write_tsv("ys.ref.exp.matrix.txt")
```

## 桃金娘目系统发育树

```
### 使用Helixer对基因组进行注释
sbatch run_helixer_Pome.slurm
## 根据gff提取pep
snakemake -s gffread.smk --cores 128 -p
## 统计pep数量
snakemake -s statPEPnum.smk --cores 128 -p

time orthofinder -f all.peps/ -M msa -S diamond -T fasttree -a 120 -t 120

mkdir iqtree.output

cp ../all.peps/OrthoFinder/Results_Apr01/MultipleSequenceAlignments/SpeciesTreeAlignment.fa ./
clipkit SpeciesTreeAlignment.fa

time iqtree2 -s SpeciesTreeAlignment.fa.clipkit -B 1000 -o Vitis_vinifera -T AUTO
```

## 基因家族收缩扩张

```
time orthofinder -f allpeps -M msa -S diamond -T fasttree -a 120 -t 120

read_tsv("../11.orthofinder/emapper/ys.emapper.annotations",comment = "##")%>%
select(`#query`,BRITE)%>%
filter(BRITE!="-")%>%select(2,1)%>%write_tsv("keggTerm2Gene.temp",col_names = FALSE)
python ../11.orthofinder/getTerm2Gene.py keggTerm2Gene.temp keggTerm2Gene.txt
```
