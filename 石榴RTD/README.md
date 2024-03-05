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

### 运行transuite

```
python ~/biotools/TranSuite-main/transuite.py Auto --gtf pomeRTD/pome_RTDmaker_output/pome.gtf --fasta pomeRTD/pome_RTDmaker_output/pome.fa --outpath transuite.output --outname pome

```
