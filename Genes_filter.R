# using TPM as input
# select genes that sum expression larger than 3
TPM.fil <- df.TPM[rowSums(df.TPM) > 3,]
# find most varince genes 
CV <- function(x){
  (sd(x)/mean(x))*100
}
genesSelected <- apply(TPM.fil, 1, CV)
# Top 8000 most variable genes
genesSelected <- tail(sort(genesSelected), 8000)
TPM.fil <- TPM.fil[names(genesSelected), ]
