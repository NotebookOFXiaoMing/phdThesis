```
minigraph -cxggs --inv no -t 24 ../11.orthofinder/all.genomes/ys.Genome.fa.masked ../11.orthofinder/all.genomes/tunisia.Genome.fa.masked ../11.orthofinder/all.genomes/bhagwa.Genome.fa.masked ../11.orthofinder/all.genomes/dabenzi.Genome.fa.masked ../11.orthofinder/all.genomes/taishanhong.Genome.fa.masked ../11.orthofinder/all.genomes/azerbaijan.Genome.fa.masked > pome.gfa

snakemake -s minigraph.smk --cores 24 -p
Rscript mergeCov.R
```
