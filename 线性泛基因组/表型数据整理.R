library(readxl)
library(tidyverse)

df<-read_tsv("D:/000博士毕业论文/表型数据/表型样本.txt",col_names = FALSE)

df

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="单果重") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/01.单果重.xlsx")


read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="单果重") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  write_tsv("D:/000博士毕业论文/表型数据/01.单果重.tsv")



read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="横纵径") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/02.横纵经.xlsx")


read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="果皮厚度") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/03.果皮厚度mm.xlsx")


read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="百粒重") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/04.百粒重.xlsx")

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="籽粒硬度") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/05.籽粒硬度.xlsx")

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="可溶性固形物") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/06.可溶性固形物.xlsx")


read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="可溶性糖") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/07.可溶性糖.xlsx")

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="总酚") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/08.总酚.xlsx")

read_excel("D:/000博士毕业论文/原始数据汇总2020-董建梅.xlsx",
           sheet="可滴定酸") %>% 
  right_join(df,by=c("品种"="X2")) %>% 
  rio::export("D:/000博士毕业论文/表型数据/09.可滴定酸.xlsx")


library(smplot2)
readxl::read_excel("D:/000博士毕业论文/表型数据/09.可滴定酸.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`可滴定酸（%）`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)

readxl::read_excel("D:/000博士毕业论文/表型数据/08.总酚.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`含量 mg/100g`,fill=X2))+
  #geom_boxplot()
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL,y="总酚 (mg/100g)")

readxl::read_excel("D:/000博士毕业论文/表型数据/07.可溶性糖.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`可溶性糖含量`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)


readxl::read_excel("D:/000博士毕业论文/表型数据/06.可溶性固形物.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`可溶性固形物`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)


readxl::read_excel("D:/000博士毕业论文/表型数据/05.籽粒硬度.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`籽粒硬度`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)

readxl::read_excel("D:/000博士毕业论文/表型数据/04.百粒重.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`百粒重`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)


readxl::read_excel("D:/000博士毕业论文/表型数据/03.果皮厚度mm.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`果皮厚mm`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)



readxl::read_excel("D:/000博士毕业论文/表型数据/02.横纵经.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`纵径mm`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)

readxl::read_excel("D:/000博士毕业论文/表型数据/02.横纵经.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`横径mm`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)

readxl::read_excel("D:/000博士毕业论文/表型数据/01.单果重.xlsx") %>% 
  mutate(X2=str_extract(X1,pattern = "Ch|Ti")) %>% 
  ggplot(aes(x=X2,y=`重量`,fill=X2))+
  #geom_boxplot()+
  sm_raincloud()+
  scale_x_discrete(labels=c("农家品系","西藏石榴"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(x=NULL)
