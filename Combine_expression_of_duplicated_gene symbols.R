#### combine duplicated genes` expression by gene symbol
# find the duplicated gene symbols
dupGenes <- TCGA.COAD.TPM.dataset$geneAnno$gene[duplicated(TCGA.COAD.TPM.dataset$geneAnno$gene)]
dupGenes <- unique(as.character(dupGenes))
# add all gene symbol to the expression matrix
COAD.tpm <- TCGA.COAD.TPM.dataset$COAD.tpm
COAD.tpm.df <- as.data.frame(COAD.tpm)
COAD.tpm.df$geneSymbol <- as.character(TCGA.COAD.TPM.dataset$geneAnno$gene)
head(COAD.tpm.df)
# i=dupGenes[4]
mergeGene <- list()
for (i in 1:length(dupGenes)) {
  dupPos <- COAD.tpm.df$geneSymbol == dupGenes[i]
  dupTable <- COAD.tpm.df[dupPos,]
  mergeGene[[i]] <- apply(dupTable[1:length(dupTable)-1], 2, sum)
  #apply(mergeGene, FUN = rbind)
}
names(mergeGene) <- dupGenes
mergeGene.mx<- do.call("rbind", mergeGene)
# remove duplicated gene symbols from data
COAD.tpm.df.dupRemo <- COAD.tpm.df[!(COAD.tpm.df$geneSymbol %in% dupGenes),]
rownames(COAD.tpm.df.dupRemo) <- COAD.tpm.df.dupRemo$geneSymbol
# remove geneSymbol column
COAD.tpm.df.dupRemo <- subset(COAD.tpm.df.dupRemo, select=-c(geneSymbol))
# add duplicated genes` expression
COAD.tpm.symbol <- rbind(COAD.tpm.df.dupRemo, mergeGene.mx)
head(COAD.tpm.symbol)













