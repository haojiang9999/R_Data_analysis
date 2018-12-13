# read test example
GSE70630_expr<-read.table(gzfile("/data8t_4/JH/scRNA_seq/GEO/GSE70630/GSE70630_processed/GSE70630_OG_processed_data_v2.txt.gz"), header=T)
# too large, so subset first 500 cells
GSE70630_test<-GSE70630_expr[, 1:1000]
GSE70630_test<-get_test_exp(GSE70630_test)
# 3.get order (please see * for additional parameters)
source("get_ordIndex.R")
GSE70630_ordIndex <- get_ordIndex(GSE70630_test, 20)
# 4.get bayes-score and mean-score
source("get_score.R")
GSE70630_score_result <- get_score(t(GSE70630_test))
head(GSE70630_score_result)
GSE70630_score_result$bayes_score
GSE70630_score_result$mean_score
# 5.plot1
source("plot.R")
plot_bayes(GSE70630_score_result$bayes_score, GSE70630_ordIndex)
plot_mean(GSE70630_score_result$mean_score, GSE70630_ordIndex)
