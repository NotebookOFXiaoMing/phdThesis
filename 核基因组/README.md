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
```
