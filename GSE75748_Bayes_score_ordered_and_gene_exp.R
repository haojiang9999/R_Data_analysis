# get cell cycle Bayes-score from reCAT and ordered by scores 
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
 
}