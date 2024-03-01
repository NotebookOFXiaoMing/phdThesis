library(readxl)
library(tidyverse)

df<-read_tsv("D:/000博士毕业论文/表型数据/表型样本.txt",col_names = FALSE)

df
read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="单果重") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/01.单果重.tsv")



read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="横纵径") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/02.横纵经.tsv")


read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="果皮厚度") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/03.果皮厚度mm.tsv")


read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="百粒重") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/04.百粒重.tsv")

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="籽粒硬度") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/05.籽粒硬度.tsv")

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="可溶性固形物") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/06.可溶性固形物.tsv")


read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="可溶性糖") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/07.可溶性糖.tsv")

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="总酚") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/08.总酚.tsv")
