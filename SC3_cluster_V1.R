library(SingleCellExperiment)
library(SC3)
library(scater)
## SC3 cluster
test1 <- runPCA(test1, ncomponents=3)
plotReducedDim(test1, "PCA", ncomponents=3, colour_by="disease_type",
               shape_by="patient")
reducedDim(test1)
plotPCA(test1, ncomponents = 3, colour_by = "disease_type",
        shape_by = "patient")
test1.sce <- sc3(test1, ks = 2:4, biology = TRUE)
# compare Normal nuclus vs hole cell
colData(example_sce)
test2<- example_sce[, example_sce$disease_type == "Normal"]
dim(test2)
plotPCA(test2, colour_by = "sample_type", shape_by = "patient")
