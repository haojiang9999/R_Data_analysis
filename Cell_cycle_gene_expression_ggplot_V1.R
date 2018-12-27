# convert data.frame to tbl format
library(tidyverse)
df_tbl<-sc_type[, rownames(subset(pData(hsmm3), cell_type == cellType))] %>% 
  rownames_to_column() %>% 
  as_tibble()
#########################
# add gene_type info
# convert data.frame to tbl
colnames(gene.type)[1]<-"rowname" 
gene.type_tbl<- gene.type %>%
  as_tibble()
# Merge the gene_typy to df_tbl
df_tbl<- left_join(df_tbl, gene.type_tbl, by = "rowname") # the by value must be the same
###########################
# Subset the data.frame select genes expression level >1 and expressed at one cell at least
markerGenes<- "POU5F1"
df_tbl.test<- df_tbl%>%
  filter(rowname == markerGenes)
## prepare for ggplot
# changes a wide data format into a long data format fot ggplot2
# just use key as each sample!!
df_tbl.test.gathered <- df_tbl.test %>%
  gather(colnames(df_tbl.test)[2:length(colnames(df_tbl.test))],
         key =  "samplename",
         value = "counts")
# check out
summary(df_tbl.test.gathered)
# add cell cycle order
# convert pData to tbl
pData(sub_type)<-cbind(pData(sub_type),ordIndex)
pData_tbl<- pData(sub_type) %>%
  as_tibble()
colnames(pData_tbl)[1]<-"samplename" 
# Merge the psuedotime to df_tbl
df_tbl.test.gathered<- left_join(df_tbl.test.gathered, pData_tbl, by = "samplename") # the by value must be the same

# expression level across samples log2 transformed counts
ggplot(df_tbl.test.gathered, aes(x=ordIndex, y=counts)) + 
 #geom_step()
geom_point()

# expression level across samples log2 transformed counts
ggplot(df_tbl.test.gathered, aes(x=ordIndex, y=log2(counts))) + 
  #geom_step()
  geom_point()
