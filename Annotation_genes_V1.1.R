library(mygene)
#load data
gene_name<-readRDS('gene.name.rds')
## Using package 'mygene' annotaed genes.
gene<-queryMany(gene_name,scopes=c("symbol","alias","hgnc"),return.as="DataFrame",fields= "type_of_gene",species="human")

#Threre were two types of missing
###find genes were not find in mygene
#if notfound is TRUE
gene[!is.na(gene$notfound),]$query
##Find genes were missed in mygene database(I don`t known why, but it happens!!! )
gene_left<-gene_name[!gene_name %in% gene$query]

##According to the genes were lefted find the annotation files
#load HGNC files downloda from https://www.genenames.org/cgi-bin/statistics
hgnc<-read.delim2('hgnc_complete_set.txt')
class(hgnc)
##Still some genes can not be annotated!!
gene_left[!gene_left %in% hgnc$symbol]

##The rest of the genes annotated by HGNC
gene_left_hgnc<-hgnc[hgnc$symbol %in% intersect(gene_left,hgnc$symbol),]
table(gene_left_hgnc$locus_group)

##Extract gene names belong to "lncRNA,ncRNA,pseudogene"
table(gene$type_of_gene)
table(gene_left_hgnc$locus_group)