# set the folder contain quant files
dir<- "/data8t_4/JH/scRNA_seq/Heart_project/ANBJ2296_PM-BJ2296-44_AHNT5NCCXY_2018-09-26/clean_map/salmon_transcripts_quant"
# extract folders name
samples<- dir("/data8t_4/JH/scRNA_seq/Heart_project/ANBJ2296_PM-BJ2296-44_AHNT5NCCXY_2018-09-26/clean_map/salmon_transcripts_quant")
samples<-samples[1:5]
# construct each quant files PATH
files <- file.path(dir,samples, "quant.sf")
# give names to these files
names(files) <- samples
all(file.exists(files))
#need to construct a tx2gene table
# using GenomicFeatures::makeTxDbFromGFF and gencodeV29.gtf file
txdb<-GenomicFeatures::makeTxDbFromGFF(file = "/data8t_4/JH/scRNA_seq/genome/gencode/release_29/gencode.v29.annotation.gtf.gz")
# check out the features
columns(txdb)
keytypes(txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene<- select(txdb,keys = k ,"GENEID", "TXNAME")
# import data sets through tximport
library(tximport)
# the default sets were get gene-level summarization
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = T) # using ignoreAfterBar to get rid of "|"
names(txi)
# to avoid gene-level summarization by setting txOut=TRUE, giving the original transcript level estimates as a list of matrices.
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
#These matrices can then be summarized afterwards using the function summarizeToGene. 
#This then gives the identical list of matrices as using  txOut=FALSE (default) in the first tximport call.
txi.sum <- summarizeToGene(txi.tx, tx2gene, ignoreAfterBar = T)
all.equal(txi$counts, txi.sum$counts)
