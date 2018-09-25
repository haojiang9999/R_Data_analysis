# import gtf file
gencode.v19.annotation.gtf<- rtracklayer::import("gencode.v19.annotation.gtf")
# select genes want to be annotated by gene_name
gene.type<-gencode.v19.annotation.gtf[gene_name %in% gencode.v19.annotation.gtf$gene_name, c("gene_name","gene_type","gene_status")]
# get a data.frame of annotation feilds
gene.type<-as.data.frame(gene.type)
# select columns want to keep
gene.type<-gene.type[, 6:7]
# remove the duplicate rows
gene.type<-gene.type[!duplicated(gene.type), ]
# Find out gene names that have different annotation
gene.type[duplicated(gene.type[, 1]), ]
# Check out!!
gene.type[gene.type$gene_name == "KMT2B", ]
# Remove duplicated gene name rows RANDOMLY
gene.type<-gene.type[!duplicated(gene.type[, 1]), ]
# Save the data
saveRDS(gene.type, file = "gene_type_hg19_by_symbol.rds")
