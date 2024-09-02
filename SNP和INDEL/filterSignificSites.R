args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
read_tsv(args[1],col_names = FALSE) %>% 
  filter(X4<=1.04e-5) %>% 
  pull(X1) %>% 
  write_lines(args[2])