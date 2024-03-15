### 统计二代测序数据的数据量
```
snakemake -s fastp_report_stat.smk --cores 2 -p

library(tidyverse)
list.files("01.fastp.report.summary",pattern = "*.csv",full.names = TRUE)%>%map(read_csv)%>%bind_rows()%>%write_csv("total_base.csv")

```
