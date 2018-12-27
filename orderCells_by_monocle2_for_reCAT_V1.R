# Pseudotime assignment using monocle 2
# load data from reCAT
# 2.data preprocessing extract genes related to 
setwd("/data8t_4/JH/scRNA_seq/Test/cell_cycle/reCAT-master/R")
source("get_test_exp.R")
GSE70630_test<-get_test_exp(GSE70630_test)
class(GSE70630_test)
pd1 <- data.frame(colnames(t(GSE70630_test))) # cell names
rownames(pd1) <- pd1[, 1]
fd1 <- data.frame(rownames(t(GSE70630_test))) # gene names
rownames(fd1) <- fd1[, 1]
pd <- new("AnnotatedDataFrame", data = pd1)
fd <- new("AnnotatedDataFrame", data = fd1)
library(monocle)
data_es <- newCellDataSet(t(GSE70630_test),
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = tobit() #Choosing a distribution for your data
                      )
# Trajectory step 3: order cells along the trajectory
ordering_genes <- rownames(t(GSE70630_test))
hsmm1 <- setOrderingFilter(data_es, ordering_genes)
hsmm2 <- reduceDimension(hsmm1, method = "DDRTree", max_components = 2)
hsmm3 <- orderCells(hsmm2, reverse = F)
plot_cell_trajectory(hsmm3)
pData(hsmm3)
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")

hsmm3<- orderCells(hsmm3)
plot_cell_trajectory(hsmm3, color_by = "Pseudotime")
# I get the pseudotime Yeah!
head(pData(hsmm3))
summary(pData(hsmm3)$Pseudotime)
length(pData(hsmm3)$Pseudotime)
table(duplicated(pData(hsmm3)$Pseudotime))

### Back to reCAT
ordIndex <- pData(hsmm3)$Pseudotime
# 4.get bayes-score and mean-score
source("get_score.R")
score_result <- get_score(t(GSE70630_test))
head(score_result)
score_result$bayes_score
score_result$mean_score
# 5.plot1
source("plot.R")
plot_bayes(score_result$bayes_score, order(pData(hsmm3)$Pseudotime))
plot_mean(score_result$mean_score, order(pData(hsmm3)$Pseudotime))
# 6.HMM
source("get_hmm.R")
myord = c(4:1, 295:5)
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
