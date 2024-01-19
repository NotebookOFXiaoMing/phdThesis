library(tidyverse)

read_tsv("02.cov/ysCov.tsv") %>% 
  left_join(read_tsv("02.cov/tunisiaCov.tsv"),
            by=c("nodeid"="nodeid")) %>% 
  left_join(read_tsv("02.cov/bhagwaCov.tsv"),
            by=c("nodeid"="nodeid")) %>% 
  left_join(read_tsv("02.cov/taishanhongCov.tsv"),
            by=c("nodeid"="nodeid")) %>% 
  left_join(read_tsv("02.cov/dabenziCov.tsv"),
            by=c("nodeid"="nodeid")) %>% 
  left_join(read_tsv("02.cov/azerbaijanCov.tsv"),
            by=c("nodeid"="nodeid")) %>% 
  select(-c(3:5)) -> datmat

datmat %>% write_tsv("datmat.tsv")
