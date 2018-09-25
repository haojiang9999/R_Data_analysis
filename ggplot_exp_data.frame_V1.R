# Load the data.frame
df<-readRDS('hg19_SMART_heart_sn.rds')
# convert data.frame factor to numeric
indx <- sapply(df, is.factor)
df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
# convert data.frame to tbl format
library(tidyverse)
df_tbl<-df %>% 
        rownames_to_column() %>% 
        as_tibble()
# add sum of counts and gene expression across cells 
df_tbl<-df_tbl %>%
        mutate(sum = rowSums(.[2:6])) %>%         # Sum of every row counts
        mutate(cell_exp = rowSums(!.[2:6] == 0))  # how many cells expression this gene by find non 0 vlues 
# add gene_type info
# convert data.frame to tbl
colnames(gene.type)[1]<-"rowname" 
gene.type_tbl<- gene.type %>%
                as_tibble()
# Merge the gene_typy to df_tbl
df_tbl<- left_join(df_tbl, gene.type_tbl, by = "rowname") # the by value must be the same
# Subset the data.frame select genes expression level >1 and expressed at one cell at least
df_tbl.test<- df_tbl%>%
              filter(sum > 0 & cell_exp >=1)
## prepare for ggplot
# changes a wide data format into a long data format fot ggplot2
# just use key as each sample!!
df_tbl.test.gathered <- df_tbl.test %>%
                        gather(colnames(df_tbl.test)[2:6],
                        key =  "samplename",
                        value = "counts")
# check out
summary(df_tbl.test.gathered)
# expression level across samples log2 transformed counts
ggplot(df_tbl.test.gathered, aes(x=samplename, y=log2(counts))) + 
  geom_boxplot() + 
  coord_flip()
# 
ggplot(df_tbl.test.gathered, aes(x=log2(counts), group=samplename)) + 
  geom_line(stat="density")
ggplot(df_tbl.test.gathered, aes(samplename, fill = gene_type.x)) + geom_bar(counts > 0)
table(df_tbl.test.gathered$gene_type.x)

df_tbl.test.gathered %>%
 filter(log2(counts) >1 )  %>%
  ggplot(aes(samplename, fill = gene_type.x, order = c("lincRNA"))) + geom_bar()+ scale_fill_brewer(palette = "Paired")
RColorBrewer::display.brewer.all()
