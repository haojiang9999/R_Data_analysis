#### TCGA PAN-Cancer data
### file path
filePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/"
# read TCGA TPM data that transformed by log2(x+0.001)
TCGA.tpm.file <- paste0(filePath,"gene_expression_RNAseq/tcga_RSEM_gene_tpm/tcga_RSEM_gene_tpm.gz")
TCGA.tpm <- read.delim2(gzfile(TCGA.tpm.file))
# TCGA.tpm[1:10,1:10]
# Check the TPM sum for each sample
y <- as.numeric(as.character(TCGA.tpm[,4]))
sum(`^`(2,y)-0.001)
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
selectedSamples <- gsub("-",".",selectedSamples)
test <- TCGA.tpm[,as.character(selectedSamples)]
table(colnames(TCGA.tpm) %in% selectedSamples)

TCGA.tpm[, as.character(selectedSamples)]
as.character(selectedSamples)
