## 代码

### megahit 组装
```
## fastp对原始数据进行过滤
snakemake -s fastp.smk --cores 52 -p
## megahit组装得到contigs
sbatch run_megahit.slurm

## megahit组装得到的contig id重新命名，在id里添加了样本名字
snakemake -s renameMegahitOutout.smk --cores 52 -p

conda activate quast 
snakemake -s quast.smk --cores 52 -p
```

### 提取没有比对上的序列

```
snakemake -s extractUnalignContigs.smk --cores 32 -p
cat 04.unalignContigs/*.fa > allUnalignContigs.fa
```

### 去除可能的线粒体和叶绿体数据
```
mkdir 05.rmCpMitolikeContig
cd 05.rmCpMitolikeContig
samtools faidx ../allUnalignContigs.fa
awk 'OFS="\t" {print $1,0,$2}' ../allUnalignContigs.fa.fai > allUnalignContigs.fa.bed
minimap2 -x asm5 -t 24 ../../sour.pome/20231015.reanalysis/02.pome_CP_Mito/pomeMito.fa ../allUnalignContigs.fa -o mito.paf
bedtools coverage -a allUnalignContigs.fa.bed -b mito.paf.bed > mito.paf.bed.cov
awk '{ if ($7>.5) print $1} ' mito.paf.bed.cov > potentialMito.ids

minimap2 -x asm5 -t 24 ../../sour.pome/20231015.reanalysis/02.pome_CP_Mito/pomeCP.fa ../allUnalignContigs.fa -o cp.paf
awk 'OFS="\t" {print $1,$3,$4}' cp.paf > cp.paf.bed
bedtools coverage -a allUnalignContigs.fa.bed -b cp.paf.bed > cp.paf.bed.cov
awk '{ if ($7>.5) print $1} ' cp.paf.bed.cov > potentialCp.ids

cat potentialCp.ids potentialMito.ids > potentialMitoCp.ids
cd ../
seqkit grep -r -v -f 05.rmCpMitolikeContig/potentialMitoCp.ids allUnalignContigs.fa -o allUnalignExcludeMiotCpLikeContigs.fa
```

### 聚类去冗余 参考https://github.com/SJTU-CGM/TGSRICEPAN/blob/master/Genomes2Pan.sh
```
 ~/biotools/gclust/script/sortgenome.pl --genomes-file ../allUnalignExcludeMiotCpLikeContigs.fa --sortedgenomes-file unalign_no2_merge.sorted.fa
time ~/biotools/gclust/gclust -loadall -minlen 20 -both -nuc -threads 32 -ext 1 -sparse 2 -memiden 90 unalign_no2_merge.sorted.fa > unalign.gclust.clu
# 9 min
python cluster2fasta.py unalign_no2_merge.sorted.fa unalign.gclust.clu unalign.gclust.fa

~/biotools/EUPAN-v0.44/bin/eupan rmRedundant blastCluster -c 0.9 -t 32 unalign.gclust.fa selfblast /home/myan/biotools/EUPAN-v0.44/tools/ncbi-blast-2.2.28+/bin
seqkit stats selfblast/non-redundant.fa

## 121,240  100,064,379
```

### NT数据库进行比对
```
python splitSeqs.py ../05.gclust/selfblast/non-redundant.fa 5000
mkdir split.seqs
mv seq* split.seqs
snakemake -s blastnNT.smk --cores 96 -p
snakemake -s getNonPlantTaxid.smk --cores 32 -p

cat nonPlantTaxid/*.taxid | sort -u | wc -l # 53744
cat nonPlantTaxid/*.taxid | sort -u > nonPlantTaxid.ids
cd ../
seqkit grep -r -v -f 06.blastn2NT/nonPlantTaxid.ids 05.gclust/selfblast/non-redundant.fa -o final.NonRefSeq.fa
seqkit stats final.NonRefSeq.fa

# 66,486  49,338,095 bp

python renameNRSseqID.py final.NonRefSeq.fa final.NonRefSeq.rename.fa
```

### 使用酸石榴基因组构建的重复序列库进行重复序列注释

```
conda activate edta03
RepeatMasker -e rmblast -pa 32 -qq -xsmall \
-lib /data/myan/raw_data/pome/sour.pome/20231015.reanalysis/06.repeatModulerRepeatMaskerafterNextpolish/00.ys.rep.lib/ys-families.fa \
/data/myan/raw_data/pome/pan.raw.fq/final.NonRefSeq.rename.fa -dir pomeNRS.repeatmasker
```

### 注释蛋白编码基因
```
conda activate galba
time perl ~/biotools/GALBA-main/scripts/galba.pl --species=pomeNRS --genome=../07.repeatMasker/pomeNRS.repeatmasker/final.NonRefSeq.rename.fa.masked --prot_seq=all.six.Pome.peps --threads 32

conda activate rnaseq
snakemake -s genomeAnnotationStep01.smk --configfiles=step01.yaml --cores 64 -p

conda activate braker2
snakemake -s genomeAnnotationStep02.smk --configfiles=step02.yaml --cores 1 -p
snakemake -s genomeAnnotationStep03.smk --configfile=step03.yaml --cores 48 -p

snakemake -s genomeAnnotationStep03_04.smk --configfile=step04.yaml --cores 2 -p
conda activate EVM
snakemake -s genomeAnnotationStep04.smk --configfile=step04.yaml --cores 32 -p
```

### 蛋白编码基因功能注释
```
conda activate eggnogmapper
export EGGNOG_DATA_DIR=/home/myan/my_data/database/eggnog/
time emapper.py -i ../05.evm/NonRefSeq/NonRefSeq.pep.fa -o NonRefSeq --cpu 24 -m diamond
conda activate interproscan01
interproscan.sh -i ../05.evm/NonRefSeq/NonRefSeq.pep.fa -f tsv --goterms -dp -cpu 24 -o NonRefSeq.pep.fa.interproscan
```


### 检测SNP和InDel

```
# 单个基因组
bwa index ../../../sour.pome/20231015.reanalysis/08.proteinCodingGenes/ys.final.masked.fna
time snakemake -s bwaSamtoolsBcftools.smk --configfiles=config_SNP.yaml --cores 128 -p

## 叶城酸石榴 加 非参考序列
cat ../../../sour.pome/20231015.reanalysis/08.proteinCodingGenes/ys.final.masked.fna ../../07.repeatMasker/pomeNRS.repeatmasker/final.NonRefSeq.rename.fa.masked > ys.final.plus.nonRefSeq.masked.fna
cat ../../../sour.pome/20231015.reanalysis/08.proteinCodingGenes/05.evm/ys/ys.rename.gff3 ../../10.evm.proteinCodingGenes/05.evm/NonRefSeq/NonRefSeq.rename.gff3 > ys.final.plus.nonRefSeq.masked.gff3

## 单个基因组加上非参考序列

snakemake -s bwaSamtoolsBcftools.smk --configfiles=config_SNP.yaml --cores 32 -p
```

### 公共转录组数据

```
# 单个参考基因组
conda activate rnaseq
snakemake -s hisat2Stringtie.smk --configfiles=PRJNA360679.yaml --cores 64 -p
snakemake -s hisat2Stringtie.smk --configfiles=PRJNA361285.yaml --cores 64 -p

```

### 计算基因的覆盖度 生成基因存在缺失矩阵

```
## gff格式转换为gtf格式
gffread ys.final.plus.nonRefSeq.masked.gff3 -T ys.final.plus.nonRefSeq.masked.tmp.gtf
## gffread转换得到的gtf文件第三列没有gene 而是transcript,应该是有参数没有设置对，暂时找不到对应参数了 下面脚本改一下
python changeGTF.py ys.final.plus.nonRefSeq.masked.tmp.gtf ys.final.plus.nonRefSeq.masked.gtf
#~/biotools/software.package/EUPAN-v0.44/bin/ccov abcd.gtf 02.eupan.input.bam/Ch_BHYSZ.sorted/Ch_BHYSZ.sorted.bam > output.dir/data/Ch_BHYSZ.sorted/Ch_BHYSZ.sorted.sta
snakemake -s geneANDcdsCOV.smk --cores 24 -p
#~/biotools/software.package/EUPAN-v0.44/bin/eupan geneCov -t 16 04.practice/ output.dir01 ys.final.plus.nonRefSeq.masked.gtf

Rscript summary_geneCov.R
Rscript summary_cdsCov.R

library(tidyverse)
gene.df<-read_tsv("summary_gene.cov")
cds.df<-read_tsv("summary_cds.cov")
identical(colnames(gene.df),colnames(cds.df))
identical(gene.df$`#name`,cds.df$`#name`)

~/biotools/software.package/EUPAN-v0.44/bin/eupan geneExist summary_gene.cov summary_cds.cov 0.8 0.95 > genePAV.matrix
```

### coinfinder

```
### iqtree 进化树

coinfinder -i coinfinder.geneinfo -p genePAV.matrix.phy.treefile -a

```

### RGA 抗病基因类似物
```
conda activate interproscan01
module load perl5/5.18.2-thread-multi

```

## 根据fastp的输出统计数据量

```
## 文件路径 /home/myan/my_data/raw_data/pome/pan.raw.fq
## 总共96份材料 92份RY的数据 + Tunisia + Nana + PT + YT 这四个样本的测序深度比较高
## 01.fastp.report.summary 92个样本
## 直接用92个数据统计
snakemake -s fastp_report_stat.smk --cores 32 -p
list.files("01.fastp.report.summary",pattern = "*.csv",full.names = TRUE)%>%map(read_csv)%>%bind_rows() -> dat
dat%>%mutate(total=r1+r2)%>%mutate(depth=total/340000000)%>%write_csv("sequencingDataStat92.csv")
```

## 统计megahit组装的结果 contig数量 N50 组装的基因组大小 BUSCO值等

```
## contig数量 基因组大小 N50
snakemake -s megahit_assemble_stat.smk --cores 32 -p
list.files("02.megahit.stats",pattern = "*.txt",full.names = TRUE,recursive = FALSE)%>%map(read_tsv)%>%bind_rows() -> dat

dat%>%filter(str_detect(file,"_"))%>%arrange(desc(num_seqs))
dat%>%filter(str_detect(file,"_"))%>%arrange(num_seqs)
```
## 统计quast 没有比对上的序列数量

```
list.files("03.quast",pattern = "*.unaligned.info",recursive = TRUE,full.names = TRUE)%>%map(read_tsv)%>%bind_rows() -> dat

## 完全没有比对的序列
dat%>%filter(str_starts(Contig,pattern = "Ch|Ti|Gl"))%>%filter(Unaligned_type=="full")%>%mutate(group=str_sub(Contig,1,2))%>%group_by(group)%>%summarise(count=n())

dat%>%filter(str_starts(Contig,pattern = "Ch|Ti|Gl"))%>%filter(Unaligned_type=="full")%>%mutate(group=str_sub(Contig,1,2))%>%group_by(group)%>%summarise(count=n())%>%pull(count)%>%sum()
dat%>%filter(str_starts(Contig,pattern = "Ch|Ti|Gl"))%>%filter(Unaligned_type=="full")%>%pull(Total_length)%>%sum()

## 局部没有比对的序列
dat%>%filter(str_starts(Contig,pattern = "Ch|Ti|Gl"))%>%filter(Unaligned_type!="full")%>%filter(Unaligned_length/Total_length>=0.5)%>%mutate(group=str_extract(Contig,pattern =".*_k"))%>%select(-5)%>%group_by(group)%>%summarise(count=n())%>%pull(count)%>%sum()

```
## 用基因组序列去跑BUSCO

```
## 测试
time busco -i 02.megahit.output.rename/Ch_BHYSZ.contigs.fa -l /home/myan/my_data/database/eukaryota_odb10/ -o 02.megahit.busco.output -m genome --force --offline --augustus --cpu 32
## 如果用太多核心可能会有报错 12个核心还好
## 5min运行完
## eukaryota_odb10 数据库是255
conda activate busco5.4.6
snakemake -s busco.smk --cores 128 -p

## eukaryota_odb10/ 数据库 
time busco -i 02.megahit.output.rename/Ch_BHYSZ.contigs.fa -l /home/myan/my_data/database/embryophyta_odb10 -o 02.megahit.busco.output -m genome --force --offline --augustus --cpu 12
## 34 min
snakemake -s busco.smk --cores 128 -p
```

## 非参考序列中蛋白编码基因长度 外显子数量

```
rtracklayer::import("10.evm.proteinCodingGenes/05.evm/NonRefSeq/NonRefSeq.rename.gff3")%>%as.data.frame()%>%filter(type=="gene")%>%pull(width)%>%mean()
rtracklayer::import("10.evm.proteinCodingGenes/05.evm/NonRefSeq/NonRefSeq.rename.gff3")%>%as.data.frame()%>%filter(type=="exon")%>%mutate(ID=str_replace(ID,pattern = ".exon[0-9]+",""))%>%group_by(ID)%>%summarise(count=n())%>%pull(count)%>%mean()

rtracklayer::import("10.evm.proteinCodingGenes/05.evm/NonRefSeq/NonRefSeq.rename.gff3")%>%as.data.frame()%>%filter(type=="gene")%>%select(width)%>%mutate(group="Non
    ref") -> nonref.gene.len

rtracklayer::import("/data/myan/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/ys.rename.gff3")%>%as.data.frame()%>%filter(type=="gene")%>%s
    elect(width)%>%mutate(group="Ref") -> ref.gene.len

t.test(width~group,data=bind_rows(ref.gene.len,nonref.gene.len))
bind_rows(ref.gene.len,nonref.gene.len)%>%write_tsv("gene.len.txt")

rtracklayer::import("/data/myan/raw_data/pome/sour.pome/20231015.reanalysis/11.orthofinder/all.genomes/ys.rename.gff3")%>%as.data.frame()%>%filter(type=="exon")%>%m
    utate(ID=str_replace(ID,pattern = ".exon[0-9]+",""))%>%group_by(ID)%>%summarise(count=n())%>%mutate(group="Ref") -> ref.exon.number

rtracklayer::import("10.evm.proteinCodingGenes/05.evm/NonRefSeq/NonRefSeq.rename.gff3")%>%as.data.frame()%>%filter(type=="exon")%>%mutate(ID=str_replace(ID,pattern 
    = ".exon[0-9]+",""))%>%group_by(ID)%>%summarise(count=n())%>%mutate(group="Ref") -> nonref.exon.number

bind_rows(ref.exon.number,nonref.exon.number)%>%write_tsv("exon.number.txt")
```

## 基因表达量

```
## 路径 pan.raw.fq/11.SNP.smallInDel/02.YS.NoRef.Genome
library(tidyverse)
myfun<-function(x){
  read_tsv(x) %>% 
    select(`Gene ID`,TPM) %>% 
    mutate(study=str_extract(x,pattern = "PRJNA[0-9]+"),
           sampleid=str_extract(x,pattern = "SRR[0-9]+")) -> dat
  return(dat)
}
list.files("00.rnaseq",pattern = "gene_abund.tsv",full.names = TRUE,recursive = TRUE)%>%map(myfun)%>%bind_rows() -> tpm.longer

tpm.longer%>%pull(sampleid)%>%unique()%>%length #335
tpm.longer%>%pull(study)%>%unique()%>%length # 28

tpm.longer%>%select(-3)%>%pivot_wider(names_from = "sampleid",values_from = "TPM")%>%column_to_rownames("Gene ID") -> tpm.wider.raw

tpm.wider<-tpm.wider.raw
tpm.wider[tpm.wider <= 0.5] <- 0

rowSums(tpm.wider)%>%as.data.frame()%>%rownames_to_column("geneid")%>%rename("tpm"=".")%>%filter(str_starts(geneid,"Non"))%>%filter(tpm>0)%>%dim
rowSums(tpm.wider)%>%as.data.frame()%>%rownames_to_column("geneid")%>%rename("tpm"=".")%>%filter(!str_starts(geneid,"Non"))%>%filter(tpm>0)%>%dim

save(tpm.longer,tpm.wider,tpm.wider.raw,file="tpm335samples.Rdata")

```

## 转录组数据的比对率

```
## 路径 pan.raw.fq/11.SNP.smallInDel/02.YS.NoRef.Genome

myfun<-function(x){
  align_rate<-read_lines(x)[15]%>%str_extract(pattern = "[0-9]+.[0-9]+")%>%as.numeric()
  sampleid<-str_extract(x,pattern = "SRR[0-9]+")
  
  return(data.frame(x=sampleid,y=align_rate))
}

list.files("00.rnaseq",pattern = ".txt",recursive = TRUE,full.names = TRUE)%>%map(myfun)%>%bind_rows()%>%mutate(group="refplusnonref") -> rpn.dat
list.files("../01.YS.Genome/00.rnaseq",pattern = ".txt",recursive = TRUE,full.names = TRUE)%>%map(myfun)%>%bind_rows()%>%mutate(group="ref") -> r.dat

bind_rows(r.dat,rpn.dat)%>%write_tsv("rnaseq_align_rate.txt")

```

## 非参考序列中注释到PFAM结构域基因的比例

```
read_tsv("../../10.evm.proteinCodingGenes/06.emapper/NonRefSeq.emapper.annotations",comment = "##")%>%select(`#query`,PFAMs)%>%filter(PFAMs!="-")%>%dim()
read_tsv("../../../sour.pome/20231015.reanalysis/09.function.annotation/01.emapper/ys.emapper.annotations",comment = "##")%>%select(`#query`,PFAMs)%>%filter(PFAMs!=
    "-")%>%dim()
```

## 非参考序列中注释的蛋白的GO富集分析

```
## 参考微信公众号推文 使用clusterProfiler对非模式植物进行注释
python parse_go_obofile.py -i go-basic.obo -o go.tb
## 下面的代码需要联网
python parse_eggNOG.py -i ../10.evm.proteinCodingGenes/06.emapper/NonRefSeq.emapper.annotations -g go.tb -O ath,osa -o nonref.output
python parse_eggNOG.py -i ../../sour.pome/20231015.reanalysis/09.function.annotation/01.emapper/ys.emapper.annotations -g go.tb -O ath,osa -o ref.output

library(clusterProfiler)

GOinfo<-read_tsv("phdthesis/chapter4/data/go.tb")
GOinfo

list.files("phdthesis",pattern = "GOanno*",recursive = TRUE,
           full.names = TRUE) %>% 
  map(read_tsv) %>% 
  bind_rows() %>% 
  filter(GO!="-") -> GOannotation

GOannotation %>% 
  filter(str_starts(gene,"NonRef")) %>% 
  pull(gene) %>% unique() -> nonrefgenelist

enricher(nonrefgenelist,
         TERM2GENE = GOannotation %>% 
           filter(level=="MF") %>% 
           select(2,1),
         TERM2NAME = GOinfo[1:2],
         pvalueCutoff = 0.5,
         qvalueCutoff = 0.5) -> nonref.go.enricher.MF

dotplot(nonref.go.enricher.MF)

enricher(nonrefgenelist,
         TERM2GENE = GOannotation %>% 
           filter(level=="BP") %>% 
           select(2,1),
         TERM2NAME = GOinfo[1:2],
         pvalueCutoff = 0.5,
         qvalueCutoff = 0.5) -> nonref.go.enricher.BP
dotplot(nonref.go.enricher.BP)

enricher(nonrefgenelist,
         TERM2GENE = GOannotation %>% 
           filter(level=="CC") %>% 
           select(2,1),
         TERM2NAME = GOinfo[1:2],
         pvalueCutoff = 0.5,
         qvalueCutoff = 0.5) -> nonref.go.enricher.CC
dotplot(nonref.go.enricher.CC)

```

## 使用resistify注释抗病基因

```
## https://github.com/SwiftSeal/resistify
conda activate resistify
resistify ../10.evm.proteinCodingGenes/05.evm/NonRefSeq/NonRefSeq.pep.fa nonref.resistify.output
## 运行11min
resistify ../08.proteinCodingGenes/05.evm/ys/ys.pep.fa ref.resistify.output
## 运行了2个多小时
```

## 基因存在缺失矩阵构建进化树

```
read_tsv("phdthesis/chapter4/data/genePAV92.matrix") %>% 
  column_to_rownames("geneid") %>%
  t() %>% 
  as.data.frame() %>% 
  unite("newcol",everything(),sep="") %>% 
  rownames_to_column() %>% 
  magrittr::set_colnames(c("92","31998")) %>% 
  write_tsv("phdthesis/chapter4/data/genePAV92.iqtree.phy")

iqtree2 -s genePAV92.iqtree.phy -bb 1000

```
