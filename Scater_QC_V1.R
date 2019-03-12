# Quality control with scater
library(scater)
example_sce <- calculateQCMetrics(scHeart_FAMTOM5, 
                                  )
colnames(colData(example_sce))
colnames(rowData(example_sce))
plotQC(example_sce, type = "highest-expression")
plotQC(example_sce, type = "exprs-freq-vs-mean")
plotColData(example_sce, x = "total_features_by_counts",
            y = "pct_counts_feature_control", colour = "Mutation_Status") +
  theme(legend.position = "top") +
  stat_smooth(method = "lm", se = FALSE, size = 2, fullrange = TRUE)
plotScater(example_sce, block1 = "Mutation_Status", block2 = "Treatment",
           colour_by = "Cell_Cycle", nfeatures = 300, exprs_values = "counts")
example_sce2 <- example_sce
example_sce2$plate_position <- paste0(
  rep(LETTERS[1:5], each = 8), 
  rep(formatC(1:8, width = 2, flag = "0"), 5)
)
plotPlatePosition(example_sce2, colour_by = "Gene_0001",
                  by_exprs_values = "counts") 
plotRowData(example_sce, x = "n_cells_by_counts", y = "mean_counts")
# Filtering the SingleCellExperiment
# By cells
keep.total <- example_sce$total_counts > 1e5
keep.n <- example_sce$total_features_by_counts > 500
filtered <- example_sce[,keep.total & keep.n]
dim(filtered)
example_sce <- runPCA(example_sce, use_coldata = TRUE,
                      detect_outliers = TRUE)
plotReducedDim(example_sce, use_dimred="PCA_coldata")
# By features
keep_feature <- nexprs(example_sce, byrow=TRUE) >= 4
example_sce <- example_sce[keep_feature,]
dim(example_sce)

example_sce <- normalize(example_sce)
plotQC(example_sce, type = "expl")
plotQC(example_sce, type = "expl",
       variables = c("total_features_by_counts", "total_counts",
                     "Mutation_Status", "Treatment", "Cell_Cycle"))
plotQC(example_sce, type = "expl", method = "pairs", theme_size = 6)


