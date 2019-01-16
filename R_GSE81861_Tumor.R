tumor.all.counts <- cbind(geneName, tumor.all.counts)
head(NM.all.counts)
# convert data.frame to tbl format
library(tidyverse)
tumor.all_tbl<-tumor.all.counts %>% 
  #  rownames_to_column() %>% 
  as_tibble()
# reform the tbl
tumor.all_tbl <- tumor.all_tbl %>% 
  select(-X)
tumor.all_tbl[,376:377]
# add sum of counts and gene expression across cells 
tumor.all_tbl<-tumor.all_tbl %>%
  mutate(sum = rowSums(.[2:376])) %>% # Sum of every row counts
  filter(sum > 0)

# add gene_type info
# convert data.frame to tbl
tumor.gene.sum <- table(gencode.v19.gene.type[tumor.all_tbl$geneName, ]$gene_type)
tumor.gene.sum <- as.data.frame(tumor.gene.sum)
tumor.lncRNA <-  tumor.gene.sum[tumor.gene.sum$Var1 %in% c("antisense",
                                                  "lincRNA",
                                                  "processed_transcript",
                                                  "pseudogene"), ]
sum(tumor.lncRNA$Freq)
lncRNA <- data.frame("lncRNA",sum(tumor.lncRNA$Freq))
names(lncRNA) <- c("Var1", "Freq")
tumor.gene.sum.sub <-rbind(tumor.gene.sum[tumor.gene.sum$Var1 == "protein_coding", ],
                        lncRNA)
# all genes expression
library(ggplot2)
# Barplot
bp<- ggplot(tumor.gene.sum.sub, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
bp
# Create a blank theme 
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
pie <- bp + coord_polar("y", start=0) + blank_theme +
  theme(axis.text.x=element_blank()) +
  ggtitle("Tumor expressed gene summary") +
  geom_text(aes(label = Freq), position = "identity",vjust = 10)

png(filename = "Tumor expressed gene summary.png",
    width = 1200, height = 703)
pie
dev.off()
# lncRNA expressionl
bp<- ggplot(tumor.lncRNA, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0) + blank_theme +
  theme(axis.text.x=element_blank()) +
  ggtitle("Tumor expressed lncRNA summary") +
  geom_text(aes(label = Freq), hjust = -0.5)

png(filename = "Tumor expressed lncRNA summary.png",
    width = 1200, height = 703)
pie
dev.off()
## Tumor tissues expression median ##
for (i in 2:268) {
  med.lnc <- NM.all_tbl[,c(1,i)]
  med.lnc <- filter(med.lnc, med.lnc[2] >0)
  med.lnc <- filter(med.lnc,geneName %in% NM.lncRNA.name$gene_id)
  med.lnc<-as.data.frame(med.lnc)
  med.lnc[,2]
  #  median(med.lnc[,2])
  
  med.mRNA <- NM.all_tbl[,c(1,i)]
  med.mRNA  <- filter(med.mRNA, med.mRNA[2] >0)
  med.mRNA  <- filter(med.mRNA, geneName %in% NM.mRNA.name$gene_id)
  med.mRNA<-as.data.frame(med.mRNA)
  #  median(med.mRNA[,2])
  #  print(i)
  NM.scMed.ratio[i-1]<-median(med.lnc[,2])/median(med.mRNA[,2])
  #  print(x)
}
