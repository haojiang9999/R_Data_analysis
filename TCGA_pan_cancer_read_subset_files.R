#### TCGA PAN-Cancer data log2(TPM + 0.001)
### file path
filePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/"
# read TCGA TPM data that transformed by log2(x+0.001)
TCGA.tpm.file <- paste0(filePath,"gene_expression_RNAseq/tcga_RSEM_gene_tpm/tcga_RSEM_gene_tpm.gz")
TCGA.tpm.log2 <- read.delim2(gzfile(TCGA.tpm.file))
# convert factors to numeric and add rownames
ENSG <- TCGA.tpm.log2$sample
TCGA.tpm.log2 <- apply(TCGA.tpm.log2[2:length(TCGA.tpm.log2)],2, 
                  function(x){
                              as.numeric(as.character(x))
                              })
rownames(TCGA.tpm.log2) <- ENSG
TCGA.tpm.log2[1:10,1:10]
### convert the data from log2(TPM+0.001) to real TPM
power.jh <- function(x, base = 2){
  `^`(base,x)
}
TCGA.tpm <- apply(TCGA.tpm.log2, 2, power.jh)
TCGA.tpm[1:10,1:10]
sum(TCGA.tpm[,3])
TCGA.tpm <- TCGA.tpm - 0.001
TCGA.tpm[TCGA.tpm<0] <- 0
saveRDS(TCGA.tpm, file = "TCGA_TPM.rds")
# Check the TPM sum for each sample
#TCGA.tpm[1:10,1:10]
### read clinical metadata
TCGA.sampelType.file <- paste0(filePath,"phenotype/TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
TCGA.sampelType <- read.table(file = TCGA.sampelType.file , sep = '\t', header = TRUE)
TCGA.survival.file <- paste0(filePath,"phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz")
TCGA.survival <- read.table(file = TCGA.survival.file , sep = '\t', header = TRUE)
# sample types
head(TCGA.sampelType)
table(TCGA.sampelType$sample_type)
z <- table(TCGA.sampelType$X_primary_disease)
length(z)
# suvival data looks contain  more info 
head(TCGA.survival)
table(TCGA.survival$cancer.type.abbreviation)
TCGA.survival$histological_type
TCGA.survival$
### subset expression and metadata
##33 cancer types 
#ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH 
#KIRC KIRP LAML  LGG LIHC LUAD LUSC MESO   OV PAAD PCPG 
#PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC  UCS  UVM
# cancer type
cancerType <- "COAD"
selectedSamples<- TCGA.survival[TCGA.survival$cancer.type.abbreviation == cancerType,]$sample
# convert - to .
selectedSamples <- gsub("-",".",selectedSamples) 
# cancer type have data in the TCGA.tpm data
selectedSamples <- selectedSamples[selectedSamples %in% colnames(TCGA.tpm)]
# subset expression matrix by cancer type
subTPM <- TCGA.tpm[, as.character(selectedSamples)]
subTPM[1:10,1:10]
# subset clinical metadata
rownames(TCGA.sampelType)<- gsub("-",".",TCGA.sampelType$sample) 
rownames(TCGA.survival)<- gsub("-",".",TCGA.survival$sample) 
TCGA.sampelType[selectedSamples,]
TCGA.survival[selectedSamples,]
subClinic <- cbind(TCGA.sampelType[selectedSamples,],TCGA.survival[selectedSamples,])
### gene annotation files
TCGA.geneAnno.file <- paste0(filePath,"gene_expression_RNAseq/tcga_RSEM_gene_tpm/gencode.v23.annotation.gene.probemap")
TCGA.geneAnno <- read.delim2(gzfile(TCGA.geneAnno.file))
head(TCGA.geneAnno)
length(TCGA.geneAnno$id)
## reorder genes by anno file
sortByname <- match(as.character(TCGA.geneAnno$id), as.character(TCGA.tpm$sample))
head(cbind(as.character(TCGA.geneAnno$id), rownames(subTPM)[sortByname]))
subTPM.order.by.anno.gene <- subTPM[sortByname,]
### Output the datasets as a list
outPut <- list(subTPM.order.by.anno.gene, subClinic, TCGA.geneAnno)
names(outPut) <- c(paste0(cancerType,".tpm"), paste0(cancerType,".clinic"), "geneAnno")
saveRDS(outPut, file = paste0("TCGA.",cancerType,".TPM.dataset.rds"))
