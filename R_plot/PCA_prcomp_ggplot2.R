MethyKit.objDB.CRC10.127cells.1Kb
meth.CRC10.1Kb <- unite(MethyKit.objDB.CRC10.127cells.1Kb)
PCASamples(meth.CRC10.1Kb, obj.return = FALSE)

# 1.mehtylKit using prcomp to do PCA and can return prcomp object
pca <- PCASamples(meth.CRC10.1Kb, obj.return = T)
head(pca)
pca$rotation # feartures contribution to each PCs
pca$x        # samples position on that direction PCs
pca$sdev     # sd of samples on the PCs direction
#2. PCA analysis using prcomp
plot(pca$x[,1], pca$x[,2])
#3. Using ggplot2 to revise this plot:
#First, a new dataframe should be created, with the information of sample-group.
pca_out <- as.data.frame(pca$x)
pca_out$group <- sapply( strsplit(as.character(row.names(pca_out)), "_"), "[[", 2 )
head(pca_out)
#Second, some prepartions for ggplot2:
library(ggplot2)
#library(grid)
#library(gridExtra)
p<-ggplot(pca_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()
p
#3.2 plot with a theme
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))
p<-ggplot(pca_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()+theme
p
#3.3 Put the words on the figure:
p<-ggplot(pca_out,aes(x=PC1,y=PC2,color=group, label=row.names(pca_out) ))
p<-p+geom_point()+ geom_text(size=3)+theme
p
#3.4 The percentage:
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- paste( colnames(pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(pca_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2])
p
#3.5 Change the order for Sample group:
#reorder the group names in the right hand by convert to factors
pca_out$group <- factor(pca_out$group, levels = c("NC", "LN1", "LN2","LN3","PT1",
                                                  "PT2","PT3","PT4"))

p<-ggplot(pca_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2]) 
#+ scale_color_manual(values=c("#FFFF00", "#00FFFF", "#FF00FF"))
p
#3.6 Save in PDF file or some formats
pdf(  "file_out.pdf",width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

#3.7 Plot features that contribute to the classification
# something wrong with this resault

pca_out_r <- as.data.frame(pca$rotation)
pca_out_r$feature <- row.names(pca_out_r)

pca_out_r

p<-ggplot(pca_out_r,aes(x=PC1,y=PC2,label=feature,color=feature ))
p<-p+geom_point()+theme + geom_text(size=3)
p


