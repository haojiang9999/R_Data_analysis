# Pseudotime assignment using monocle 2
# load data from reCAT
# 2.data preprocessing extract genes related to 
setwd("/data8t_4/JH/scRNA_seq/Test/cell_cycle/reCAT-master/R")
source("get_test_exp.R")
cell_cycle<-get_test_exp(sc_type)
class(cell_cycle)
pd1 <- data.frame(colnames(t(cell_cycle))) # cell names
rownames(pd1) <- pd1[, 1]
fd1 <- data.frame(rownames(t(cell_cycle))) # gene names
rownames(fd1) <- fd1[, 1]
library(monocle)
pd <- new("AnnotatedDataFrame", data = pd1)
fd <- new("AnnotatedDataFrame", data = fd1)
data_es <- newCellDataSet(t(cell_cycle),
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = tobit() #Choosing a distribution for your data
)
# Trajectory step 3: order cells along the trajectory
ordering_genes <- rownames(t(cell_cycle))
hsmm1 <- setOrderingFilter(data_es, ordering_genes)
hsmm2 <- reduceDimension(hsmm1, method = "DDRTree", max_components = 2)
hsmm3 <- orderCells(hsmm2, reverse = F)
plot_cell_trajectory(hsmm3)
pData(hsmm3)
plot_cell_trajectory(hsmm3, color_by = "Pseudotime")
# add cell type to data
cell_type <- stringr::str_extract(rownames(pData(hsmm3)), "^[^_]+(?=_)")
pData(hsmm3) <-cbind(pData(hsmm3),cell_type)

hsmm3<- orderCells(hsmm3)
plot_cell_trajectory(hsmm3,
                     color_by = "cell_type")
# plot score on monocle psuedotime
score_result <- get_score(exprs(hsmm3))
print(plot_bayes(score_result$bayes_score, order(pData(hsmm3)$Pseudotime)))
print(plot_mean(score_result$mean_score, order(pData(hsmm3)$Pseudotime)))
# plot score by Bayes value
score_result <- get_score(exprs(hsmm3))
print(plot_bayes(score_result$bayes_score, order(score_result$bayes_score$G2M.score)))
print(plot_mean(score_result$mean_score, order(score_result$bayes_score$G2M.score)))

############################################################################################
#monocle2 order cells
filePath = "/data8t_4/JH/scRNA_seq/GEO/GSE75748_Hum_embryonic_stem_cell/R_GSE75748_Hum_embryonic_stem_cell/"
# plot different cell types
table(pData(hsmm3)$cell_type)  #how many cell in each cell types
cellType <- "H1"
for (cellType in as.character(unique(pData(hsmm3)$cell_type))){
sub_type<-hsmm3[, rownames(subset(pData(hsmm3), cell_type == cellType))]
sub_type <- reduceDimension(sub_type, method = "DDRTree", max_components = 2)
sub_type <- orderCells(sub_type, reverse = F)
png(file=paste0(filePath,cellType,"_monocle_trajectory.png"), width = 1260, height = 667)
print(plot_cell_trajectory(sub_type))
print(title(main = paste(cellType,"_monocle_trajectory.png")))
dev.off()
png(file=paste0(filePath,cellType,"_monocle_trajectory_psuedotime.png"), width = 1260, height = 667)
print(plot_cell_trajectory(sub_type, color_by = "Pseudotime"))
dev.off()
source("get_score.R")
score_result <- get_score(exprs(sub_type))
#head(score_result)
#score_result$bayes_score
#score_result$mean_score
# 5.plot1
source("plot.R")
png(file=paste0(filePath,cellType,"_monocle_Bayes.png"), width = 1260, height = 667)
print(plot_bayes(score_result$bayes_score, order(pData(sub_type)$Pseudotime)))
dev.off()
png(file=paste0(filePath, cellType,"_monocle_Mean.png"), width = 1260, height = 667)
print(plot_mean(score_result$mean_score, order(pData(sub_type)$Pseudotime)))
dev.off()
}

#############################################################################################
# reCAT order cells
#reCAT
filePath = "/data8t_4/JH/scRNA_seq/GEO/GSE75748_Hum_embryonic_stem_cell/R_GSE75748_Hum_embryonic_stem_cell/"
cellType <- "H1"
for (cellType in as.character(unique(pData(hsmm3)$cell_type))){
  print(cellType)
  sub_type<-hsmm3[, rownames(subset(pData(hsmm3), cell_type == cellType))]
  source("get_test_exp.R")
  test<-exprs(sub_type)
  test_exp <- get_test_exp(test)
  test<-test_exp[rowSums(test_exp) > 0, colSums(test_exp) > 1] # make sure rowSums and colSums >0 but I do not know why colSums > 1 it works
  test<-as.matrix(test)                                        # make sure it`s a matrix
  table(colSums(test_exp) > 0)
  source("get_ordIndex.R")
  ordIndex <- get_ordIndex(test , 20)
  source("get_score.R")
  score_result <- get_score(t(test))
  #head(score_result)
  #score_result$bayes_score
  #score_result$mean_score
  print(paste0(cellType,"_plotting"))
  png(file=paste0(filePath,cellType,"_reCAT_Bayes.png"), width = 1260, height = 667)
  # have to add print()
  print(plot_bayes(score_result$bayes_score, ordIndex))
  dev.off()
  png(file=paste0(filePath, cellType,"_reCAT_Mean.png"), width = 1260, height = 667)
  print(plot_mean(score_result$mean_score, ordIndex))
  dev.off()
  print(paste0(cellType,"_finished"))
}
#########################Order cell by Bayes score############
for (cellType in as.character(unique(pData(hsmm3)$cell_type))){
  print(cellType)
  sub_type<-hsmm3[, rownames(subset(pData(hsmm3), cell_type == cellType))]
  source("get_test_exp.R")
  test<-exprs(sub_type)
  test_exp <- get_test_exp(test)
  test<-test_exp[rowSums(test_exp) > 0, colSums(test_exp) > 1] # make sure rowSums and colSums >0 but I do not know why colSums > 1 it works
  test<-as.matrix(test)                                        # make sure it`s a matrix
  table(colSums(test_exp) > 0)
  source("get_score.R")
  score_result <- get_score(t(test))
  ordIndex <- order(score_result$bayes_score$G2M.score)
  #head(score_result)
  #score_result$bayes_score
  #score_result$mean_score
  print(paste0(cellType,"_plotting"))
  png(file=paste0(filePath,cellType,"_order_by_Bayes_Bayes.png"), width = 1260, height = 667)
  # have to add print()
  print(plot_bayes(score_result$bayes_score, ordIndex))
  dev.off()
  png(file=paste0(filePath, cellType,"_order_by_Bayes_Mean.png"), width = 1260, height = 667)
  print(plot_mean(score_result$mean_score, ordIndex))
  dev.off()
  print(paste0(cellType,"_finished"))
}
### Back to reCAT
ordIndex <- pData(hsmm3)$Pseudotime
# 4.get bayes-score and mean-score
source("get_score.R")
score_result <- get_score(t(cell_cycle))
head(score_result)
score_result$bayes_score
score_result$mean_score
# 5.plot1
source("plot.R")
plot_bayes(score_result$bayes_score, order(pData(hsmm3)$Pseudotime))
plot_mean(score_result$mean_score, order(pData(hsmm3)$Pseudotime))
# 6.HMM
source("get_hmm.R")
#myord = c(4:1, 295:5)

hmm_result <- get_hmm_order(bayes_score = score_result$bayes_score, 
                            mean_score = score_result$mean_score, 
                            ordIndex = ordIndex, cls_num = 4, 
                            #                            myord = myord, 
                            rdata = rdata
)
# 7.choose the start
source("get_start.R")
start = get_start(bayes_score = score_result$bayes_score, 
                  mean_score = score_result$mean_score, 
                  ordIndex = ordIndex, cls_num = 3, rdata = rdata, nthread = 3)

#8.plot2
source("plot.R")
load("../data/ola_mES_2i_hmm.RData")
plot_bayes(score_result$bayes_score, ordIndex, cls_result = hmm_result, cls_ord = myord, colorbar = 1)
plot_mean(score_result$mean_score, ordIndex, hmm_result, hmm_order, 1)
# 9.cluster
source("get_clsuter_result.R")
load("../data/Flo_test_exp.RData")
load("../data/bestEnsembleComplexTSP 10 - 216 Flo .RData")
cls_result = get_cluster_result(test_exp = test_exp, ensembleResultLst = ensembleResultLst, 
                                resultLst = resultLst, cls_num = 20)
cls_result$cls
table(cls_result$cls)
