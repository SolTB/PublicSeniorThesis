### Thesis: RNA-Seq Basic Statistics
## Author: Sol Taylor-Brill
## Version: 10/02/19

#read RNA Input file
RNA<-read.csv(file = "VehicleCBDStd.csv", 
              header=TRUE)

head(RNA)
x<-length(RNA$X)

hist(RNA$lfcSE)
hist(RNA$log2FoldChange)
plot(RNA$log2FoldChange,RNA$pvalue, main="pvalue vs. log(foldchange)", xlab="log(foldchange)", ylab="p values")
#plot(RNA$X, RNA$log2FoldChange, main="Scatterplot of FoldChanges",
     #xlab="Genes ", ylab="log(FoldChange)", pch=19)
#n<-dim(RNA)[1]
#xbar<-mean(RNA$baseMean)
#xbar
#sd<-sd(RNA$baseMean)
#sd
