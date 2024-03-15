library(tidyverse)
dat<-jsonlite::fromJSON(snakemake@input[[1]])
read1<-dat$read1_after_filtering$total_bases
read2<-dat$read2_after_filtering$total_bases



data.frame(sample_id=snakemake@params[[1]],
           r1=read1,
           r2=read2) %>% 
  write_csv(snakemake@output[[1]])