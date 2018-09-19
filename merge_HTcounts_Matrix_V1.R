# not tested
files <- list.files(path="/data8t_4/JH/scRNA_seq/E_MTAB_3929/FANTOM_gene_name_ssno")
genes <- read.table(files[1], header=FALSE, sep="\t")[,1]     # gene names
df    <- do.call(cbind,lapply(files,function(fn)read.table(fn,header=FALSE, sep="\t")[,2]))
df    <- cbind(genes,df)
sample.name<-gsub("readCounts_","",files)
sample.name<-gsub(".txt","",sample.name)
colnames(df)<-c("genes",sample.name)
head(sample.name)
dim(df)
class(df)
saveRDS(df,file = "FANTOM_E_MTAB_3929_exp.rds")
