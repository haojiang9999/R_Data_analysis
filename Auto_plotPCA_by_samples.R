#read files
CGfiltered.rds <- "/stor/jianghao/Myjobs/MethyKit_read_samples/"
MethyKit.objDBs <- list.files(CGfiltered.rds, pattern = "filtered.1Kb.rds")

#i=9
#read rds files for plot 
for (i in 2:length(MethyKit.objDBs)) {
  library(methylKit)
  objDB <- readRDS(MethyKit.objDBs[i])
  meth.CRC.1Kb <- unite(objDB)
  objDB.250Kb <- tileMethylCounts(objDB, win.size=250000 , step.size = 250000, mc.cores = 8)
  meth.CRC.250Kb <- unite(objDB.250Kb)
  ########################### 1Kb bins PCA analysis #################################
  # 1.mehtylKit using prcomp to do PCA and can return prcomp object
  pca <- PCASamples(meth.CRC.1Kb, obj.return = T)
  #2. PCA analysis using prcomp
  #plot(pca$x[,1], pca$x[,2])
  #3. Using ggplot2 to revise this plot:
  #First, a new dataframe should be created, with the information of sample-group.
  pca_out <- as.data.frame(pca$x)
  pca_out$group <- sapply( strsplit(as.character(row.names(pca_out)), "_"), "[[", 2 )
  #head(pca_out)
  #Second, some prepartions for ggplot2:
  library(ggplot2)
  #library(grid)
  #library(gridExtra)
  #3.4 The percentage:
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste( colnames(pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  p.1Kb<-ggplot(pca_out,aes(x=PC1,y=PC2,color=group ))
  p.1Kb<-p.1Kb+geom_point() + xlab(percentage[1]) + ylab(percentage[2]) +
    labs(title = MethyKit.objDBs[i],
         subtitle = "1Kb windows")
  ########################### 250Kb bins PCA analysis #################################
  pca <- PCASamples(meth.CRC.250Kb, obj.return = T)
  pca_out <- as.data.frame(pca$x)
  pca_out$group <- sapply( strsplit(as.character(row.names(pca_out)), "_"), "[[", 2 )
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste( colnames(pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  p.250Kb<-ggplot(pca_out,aes(x=PC1,y=PC2,color=group ))
  p.250Kb<-p.250Kb+geom_point() + xlab(percentage[1]) + ylab(percentage[2]) +
    labs(title = MethyKit.objDBs[i],
         subtitle = "250Kb windows")
  png(filename = paste0("PCA_plot_of_",MethyKit.objDBs[i],".png"),
      width = 500, height = 1000)
  # a function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  multiplot(p.1Kb, p.250Kb)
  dev.off()
}

  
