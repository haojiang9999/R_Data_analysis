##import gtf files
FANTOM.gtf <- rtracklayer::import("FANTOM_CAT.lv3_robust.gtf")
##could be transformed to data.frame 
#gtf_df=as.data.frame(gtf)
x<-FANTOM.gtf$gene_name
lenght(unique(x))
length(unique(x))
length(intersect(genes,unique(x)))
table(is.na(x))
table(FANTOM.gtf$coding_status)
table(FANTOM.gtf$geneSuperClass)
##classified