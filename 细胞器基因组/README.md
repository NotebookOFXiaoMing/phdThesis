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
halign -t 2 downLoad.fa
```

### 叶绿体基因组circos图

http://47.96.249.172:16085/cpgview/drawmap

170884098974144

