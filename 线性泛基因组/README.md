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
time emapper.py -i ../05.evm/NonRefSeq/NonRefSeq.pep.fa -o NonRefSeq --cpu 24 -m diamond
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
