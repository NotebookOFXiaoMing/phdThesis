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
