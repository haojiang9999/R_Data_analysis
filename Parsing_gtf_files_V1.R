##import gtf files
FANTOM5.gtf <- rtracklayer::import("/data8t_4/JH/scRNA_seq/genome/FANTOM5/lv3_robust/FANTOM_CAT.lv3_robust.gtf")
##could be transformed to data.frame 
gtf_df=as.data.frame(FANTOM5.gtf)
head(gtf_df)
## parse the big table
# convert data.frame to tbl format
library(tidyverse)
gtf_df_tbl<-gtf_df %>% 
            rownames_to_column() %>% 
            as_tibble()
# subset genes not transcripts
gtf_df_tbl_genes<-gtf_df_tbl %>%
                  filter(type == "gene") 
# Remove the NA columns in the table
gtf_df_tbl_genes<-gtf_df_tbl_genes%>%
                  select(-c(9,10,21:25))
# Summary info
table(gtf_df_tbl_genes[, 10])
# save data
saveRDS(gtf_df_tbl_genes, file = "Genes_FANTOM5_CAT.lv3_robust_gtf_tbl.rds")
saveRDS(gtf_df_tbl, file = "FANTOM5_CAT.lv3_robust_gtf_tbl.rds")
write.csv(gtf_df_tbl_genes, file = "Genes_FANTOM5_CAT.lv3_robust_gtf_tbl.csv")
