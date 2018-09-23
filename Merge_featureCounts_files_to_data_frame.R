##Merge featureCoubts output
files <- list.files(path="/data8t_4/JH/scRNA_seq/Heart_project/ANBJ2296_PM-BJ2296-43_AHNJK2CCXY_2018-09-19/clean_map/",
                    pattern="counts.txt$")
genes <- read.table(files[1], header=FALSE, sep="\t")[,1]     # gene names
#Using multiple core to speed up
library(parallel)
#find cores number
detectCores(logical = F)
#set core numbers
mc <- getOption("mc.cores", 20)
#Using the mclapply
##cbind() will change the data some time so cbind.data.frame maybe help
df<- do.call(cbind.data.frame,mclapply(files,function(fn)read.table(fn,header=F, sep="\t")[,7],mc.cores = mc))
#Stop
stopCluster(mc)
#add gene names
row.names(df)<- genes
# Assign headers based on existing row in dataframe
colnames(df)<-as.character(unlist(df[1,])) #The key here is to unlist the row first.
df<- df[-1,]
# Save the files
saveRDS(df, file = "hg19_SMART_heart_sn.rds")

