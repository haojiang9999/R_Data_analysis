# not tested
files <- list.files(path="/data8t_4/JH/scRNA_seq/E_MTAB_3929/FANTOM_gene_name_ssno")
genes <- read.table(files[1], header=FALSE, sep="\t")[,1]     # gene names
#Using multiple core to speed up
library(parallel)
#find cores number
detectCores(logical = F)
#set core numbers
mc <- getOption("mc.cores", 16)
#Using the mclapply
df    <- do.call(cbind,mclapply(files,function(fn)read.table(fn,header=FALSE, sep="\t")[,2],mc.cores = mc))
#Stop
stopCluster(mc)
#add gene names
row.names(df)<- genes
#add colnames
sample.name<-gsub("readCounts_","",files)
sample.name<-gsub(".txt","",sample.name)
colnames(df)<-sample.name
head(sample.name)
dim(df)
class(df)
#Save the data
saveRDS(df,file = "FANTOM_E_MTAB_3929_exp.rds")

