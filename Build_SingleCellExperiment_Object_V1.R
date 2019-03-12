# build up a SingleCellExperiment object
library(scater)
library(SingleCellExperiment)
# read the sample meta data colData from csv file
library(readr)
scHEART_metaData <- read_csv("scHEART_metaData.csv")
View(scHEART_metaData)
scHEART_metaData<-as.data.frame(scHEART_metaData) # convert meata data to data frame
rownames(scHEART_metaData)<-scHEART_metaData[, 1] # add row names to the data frame and row names was column name in the exp Matrix  
# read the rows/genes meta data
Genes_FANTOM5_CAT.lv3_robust_gtf_tbl <- readRDS("/data8t_4/JH/scRNA_seq/genome/FANTOM5/lv3_robust/R_FANTOM5_gtf_files_parse
                                                /Genes_FANTOM5_CAT.lv3_robust_gtf_tbl.rds")
View(Genes_FANTOM5_CAT.lv3_robust_gtf_tbl)
rowData<-as.data.frame(Genes_FANTOM5_CAT.lv3_robust_gtf_tbl)       # convert tbl data to data.frame
rownames(rowData)<-rowData[, "gene_id"]                            # add row names to the rowData and row name was row.name in the exp Matrix
rowData<- rowData[rownames(txi.sum$abundance), ]                   # subset the rows/genes in tne exp Matrix
rowData<-GenomicRanges::makeGRangesFromDataFrame(rowData, keep.extra.columns = T) # convert data.frame to  GRanges object
# build SingleCellExperiment object
scHeart_FAMTOM5 <- SingleCellExperiment(
  assays = list(TMP = txi.sum$abundance,             # TMP from salmon and tximport
                counts = txi.sum$counts,             # NumReads in salmon
                log2TMP = log2(txi.sum$abundance)),  # log2 transformed TMP
  colData = scHEART_metaData,
  rowData = rowData
)
assay(scHeart_FAMTOM5, "log2TMP")
