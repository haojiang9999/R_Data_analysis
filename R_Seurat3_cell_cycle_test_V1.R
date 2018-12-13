############ Seruat 3.0 ############

library(Seurat)

# Read in the expression matrix The first row is a header row, the first
# column is rownames
exp.mat <- read.table(file = "/data8t_4/JH/scRNA_seq/Test/cell_cycle/nestorawa_forcellcycle_expressionMatrix.txt", 
                      header = TRUE, as.is = TRUE, row.names = 1)

# Also read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "/data8t_4/JH/scRNA_seq/Test/cell_cycle/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(counts = exp.mat,project = "marrow")
marrow <- NormalizeData(object = marrow)
marrow <- FindVariableFeatures(object = marrow, do.plot = FALSE, display.progress = FALSE)
marrow <- ScaleData(object = marrow, display.progress = FALSE)
# run a PCA on our object, using the variable genes we found in FindVariableFeatures above
marrow <- RunPCA(object = marrow, pc.genes = marrow@var.genes, pcs.print = 1:4, 
                 genes.print = 10)
DimHeatmap(object = marrow, dims = 5, cells = 500, balanced = TRUE)
## Assign Cell-Cycle Scores
marrow <- CellCycleScoring(object = marrow, s.features = s.genes, g2m.features = g2m.genes, 
                           set.ident = TRUE)
# view cell cycle scores and phase assignments
head(x = marrow@meta.data)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
          ncol = 2)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells
# separate entirely by phase
marrow <- RunPCA(object = marrow, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
DimPlot(object = marrow)
## Regress out cell cycle scores during data scaling
marrow <- ScaleData(object = marrow, vars.to.regress = c("S.Score", "G2M.Score"), 
                    display.progress = FALSE)
# Now, a PCA on the variable genes no longer returns components associated
# with cell cycle
marrow <- RunPCA(object = marrow, pc.genes = marrow@var.genes, genes.print = 10)
# When running a PCA on only cell cycle genes, cells no longer separate by
# cell-cycle phase
marrow <- RunPCA(object = marrow, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
DimPlot(object = marrow)
## Alternate Workflow
marrow@meta.data$CC.Difference <- marrow@meta.data$S.Score - marrow@meta.data$G2M.Score
marrow <- ScaleData(object = marrow, vars.to.regress = "CC.Difference", display.progress = FALSE)

# cell cycle effects strongly mitigated in PCA
marrow <- RunPCA(object = marrow, pc.genes = marrow@var.genes, genes.print = 10)
# when running a PCA on cell cycle genes, actively proliferating cells
# remain distinct from G1 cells however, within actively proliferating
# cells, G2M and S phase cells group together
marrow <- RunPCA(object = marrow, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
DimPlot(object = marrow)
