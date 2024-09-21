library(tidyverse)

myfun<-function(x){
  return(read_tsv(x) %>% 
    select(`Gene ID`,`Coverage`) %>% 
    mutate(sampleName=str_extract(x,pattern = "SRR[0-9]+")))
}
list.files("03.expression",
           pattern = "gene_abund.tsv",
           recursive = TRUE,
           full.names = TRUE) %>% 
  map(myfun) %>% 
  bind_rows() %>% 
  pivot_wider(names_from = "sampleName",values_from = "Coverage") %>% 
  write_tsv("coverage.tsv")
  