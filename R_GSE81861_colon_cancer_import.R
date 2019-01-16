# R step_1 read expression table from GEO download data
# 1) read from txt files
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE81861_human_colorectal_tumors/"
NM.all.counts <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_NM_all_cells_COUNT.csv.gz")), header = T)
tumor.all.counts <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_tumor_all_cells_COUNT.csv.gz")), header = T)
# 2) gene annotation
head(NM.all.counts[, 1])
sub(".*_", "", head(NM.all.counts[, 1]))
sub(".*_", "", head(tumor.all.counts[, 1]))
head(rownames(NM.all.counts))
head(rownames(tumor.all.counts))
# cell information
colnames(NM.all.counts)
colnames(tumor.all.counts)
# gene annotation
gencode.v19.gene.type
geneName <- sub(".*_", "", NM.all.counts[, 1])
match(geneName, gencode.v19.gene.type[, 2])
table(match(geneName, gencode.v19.gene.type[, 2]))
table(geneName %in% gencode.v19.gene.type[, 2]) # how many genes have annotation all
table(gencode.v19.gene.type$gene_type)
# log2 transform
NM.all.counts.log2 <- log2(NM.all.counts[, 2:267] +1 )
head(NM.all.counts)
NM.all.counts.log2 <- cbind(geneName, NM.all.counts.log2)

head(NM.all.counts.log2)
length(unique(geneName))
# convert data.frame to tbl format
library(tidyverse)
NM.all_tbl<-NM.all.counts.log2 %>% 
#  rownames_to_column() %>% 
  as_tibble()

# reform the tbl
NM.all_tbl <- NM.all_tbl %>% 
  select(-X)

# add sum of counts and gene expression across cells 
NM.all_tbl<-NM.all_tbl %>%
  mutate(sum = rowSums(.[2:267])) %>% # Sum of every row counts
  filter(sum > 0)

# add gene_type info
# convert data.frame to tbl
NM.gene.sum <- table(gencode.v19.gene.type[NM.all_tbl$geneName, ]$gene_type)
NM.gene.sum <- as.data.frame(NM.gene.sum)
NM.lncRNA <-  NM.gene.sum[NM.gene.sum$Var1 %in% c("antisense",
                                                  "lincRNA",
                                                  "processed_transcript",
                                                  "pseudogene"), ]
sum(NM.lncRNA$Freq)
lncRNA <- data.frame("lncRNA",sum(NM.lncRNA$Freq))
names(lncRNA) <- c("Var1", "Freq")
NM.gene.sum.sub <-rbind(NM.gene.sum[NM.gene.sum$Var1 == "protein_coding", ],
                        lncRNA)
library(ggplot2)
# Barplot
bp<- ggplot(NM.gene.sum.sub, aes(x="", y=Freq, fill=Var1))+
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
  ggtitle("Normal expressed gene summary") +
  geom_text(aes(label = Freq), position = "identity",vjust = 10)

png(filename = "Normal expressed gene summary.png",
     width = 1200, height = 703)
pie
dev.off()
# lncRNA expressionl
library(ggrepel)
bp<- ggplot(NM.lncRNA, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0) + blank_theme +
  theme(axis.text.x=element_blank()) +
  ggtitle("Normal expressed lncRNA summary") +
  geom_text(aes(label = Freq), hjust = -0.5)

png(filename = "Normal expressed lncRNA summary.png",
    width = 1200, height = 703)
pie
dev.off()
# mRNA and lncRNA median compare

the_median <- apply(log2(NM.all_tbl[,2:267]+1), 1, median)
NM.all_tbl$Median <- the_median
NM.all_tbl[,267:269]
summary(the_median)
# lncRNA genes
NM.lncRNA.name <- gencode.v19.gene.type[gencode.v19.gene.type$gene_type %in% c("antisense",
                                                             "lincRNA",
                                                             "processed_transcript",
                                                             "pseudogene"), ]
# mRNA genes
NM.mRNA.name<- gencode.v19.gene.type[gencode.v19.gene.type$gene_type %in% "protein_coding", ]

# bulk-like lncRNA expression median
the_median <- NM.all_tbl %>% 
              filter(geneName %in% NM.lncRNA.name$gene_id) %>% 
              select(sum) 
median(the_median$sum)
# bulk-like mRNA median expression
NM.mRNA.name<- gencode.v19.gene.type[gencode.v19.gene.type$gene_type %in% "protein_coding", ]
the_median <- NM.all_tbl %>% 
              filter(geneName %in% NM.mRNA.name$gene_id) %>% 
              select(sum) 
median(the_median$sum)
### for single cell lncRNA median vs mRNA median
med.lnc <- NM.all_tbl[,1:2] %>%
            filter(NM.all_tbl[2] >0) %>%
            filter(geneName %in% NM.lncRNA.name$gene_id)
med.lnc<-as.data.frame(med.lnc)
med.lnc[,2]
median(med.lnc[,2])

med.mRNA <- NM.all_tbl[,1:2] %>%
  filter(NM.all_tbl[2] >0) %>%
  filter(geneName %in% NM.mRNA.name$gene_id)
med.mRNA<-as.data.frame(med.mRNA)
median(med.mRNA[,2])

median(med.lnc[,2])/median(med.mRNA[,2])

for (i in 2:268) {
  med.lnc <- NM.all_tbl[,c(1,i)]
  med.lnc <- filter(med.lnc, med.lnc[2] >0)
  med.lnc <- filter(med.lnc,geneName %in% NM.lncRNA.name$gene_id)
  med.lnc<-as.data.frame(med.lnc)
  med.lnc[,2]
  median(med.lnc[,2])
  
  med.mRNA <- NM.all_tbl[,c(1,i)]
  med.mRNA  <- filter(med.mRNA, med.mRNA[2] >0)
  med.mRNA  <- filter(med.mRNA, geneName %in% NM.mRNA.name$gene_id)
  med.mRNA<-as.data.frame(med.mRNA)
  median(med.mRNA[,2])
  print(i)
  x<-median(med.lnc[,2])/median(med.mRNA[,2])
  print(x)
}


