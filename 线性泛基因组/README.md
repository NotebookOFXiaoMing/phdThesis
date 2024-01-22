## 代码

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
```

### NT数据库进行比对
```
python splitSeqs.py ../clusterRes_rep_seq.fasta 5000
```
