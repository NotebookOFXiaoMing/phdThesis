read_tsv("D:/000博士毕业论文/表型数据/表型样本.txt",col_names = FALSE)

read_tsv("D:/Jupyter/panPome/genePAV.matrix") %>% 
  select("geneID",read_tsv("D:/000博士毕业论文/表型数据/表型样本01.txt",col_names = FALSE) %>% 
           pull(X1)) %>% 
  #column_to_rownames("geneID") %>% 
  rowwise(geneID) %>% 
  mutate(total_gene = sum(c_across(where(is.numeric)))) %>% 
  filter(total_gene<=24&total_gene>=3) %>% 
  ungroup() %>% 
  select(-total_gene) -> dat.pav

dat.pav %>% dim()
read_excel("D:/000博士毕业论文/表型数据/01.单果重.xlsx") %>% 
  group_by(X1) %>% 
  summarise(mean_value=mean(重量)) %>% 
  ungroup() -> dat.pheno.01
dat.pheno.01
new.dat.pav<-dat.pav %>% 
  column_to_rownames("geneID")

rownames(new.dat.pav)[1]
new.dat.pav[1,] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  right_join(dat.pheno.01,by=c("rowname"="X1")) %>% 
  wilcox.test(mean_value ~ get("ys001G000065"),data=.) %>% 
  .$p.value

pheno.01.pvalue<-data.frame(gene_name=character(),
                            pvalue=numeric())

pheno.01.pvalue
for(i in 1:nrow(dat.pav)){
  
  pheno.01.pvalue<-add_row(pheno.01.pvalue,
                           gene_name=rownames(new.dat.pav)[i],
                           pvalue=new.dat.pav[i,] %>% 
                             t() %>% 
                             as.data.frame() %>% 
                             rownames_to_column() %>% 
                             right_join(dat.pheno.01,by=c("rowname"="X1")) %>% 
                             wilcox.test(mean_value ~ get(rownames(new.dat.pav)[i]),data=.) %>% 
                             .$p.value)
}

pheno.01.pvalue %>% 
  mutate(fdr=p.adjust(pvalue,method = "fdr")) %>% 
  ggplot(aes(x="A",y=-log10(fdr)))+
  geom_boxplot()

p.adjust(pheno.01.pvalue %>% pull(pvalue),
         method = "fdr")
ggplot(data=pheno.01.pvalue,aes(x="A",y=-log10(pvalue)))+
  geom_boxplot()
