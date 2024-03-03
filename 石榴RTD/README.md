### 参考链接

https://github.com/wyguo/genotype_specific_RTD/tree/main

https://github.com/anonconda/RTDmaker/tree/main

### 参考基因组构建索引

```
time STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ref.index --genomeFastaFiles ref/pome.ref.nonref.fna --outFileNamePrefix pome --limitGenomeGenerateRAM 240000000000
```

### star第一轮比对

这里如果是压缩文件会报错 FATAL ERROR in input reads: unknown file format: the read ID should start with @ or >

https://github.com/alexdobin/STAR/issues/379 需要先把文件解压

```
snakemake -s starPassOne.smk --cores 128 -p
```
