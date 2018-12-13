#### reCAT #####


# 2.data preprocessing
setwd("/data8t_4/JH/scRNA_seq/Test/cell_cycle/reCAT-master/R")
source("get_test_exp.R")
load("../data/ola_mES_2i.RData")
cycle.genes <- get_test_exp(test_exp)
cycle.genes<-row.names(cycle.genes)
load("../data/ola_mES_2i.RData")
head(test_exp)
length(test_exp)
test_exp<-test_exp[,cycle.genes]
# 3.get order (please see * for additional parameters)
source("get_ordIndex.R")
ordIndex <- get_ordIndex(test_exp, 20)

# 4.get bayes-score and mean-score
source("get_score.R")
score_result <- get_score(t(test_exp))
head(score_result)
score_result$bayes_score
score_result$mean_score
# 5.plot1
source("plot.R")
plot_bayes(score_result$bayes_score, ordIndex)
plot_mean(score_result$mean_score, ordIndex)
# 6.HMM
source("get_hmm.R")
load("../data/ola_mES_2i_ordIndex.RData")
load("../data/ola_mES_2i_region.RData")
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
