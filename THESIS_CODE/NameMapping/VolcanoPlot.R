res <- read.csv("p<0.001fgseaCSV.csv", header=TRUE)
head(res)

# Make a basic volcano plot
with(res, plot(ES, -log10(pval), pch=20, main="Gene Set Enrichment Pathways", xlim=c(-1,1), ylim = c(0,4)))

# Add colored points: red if pval<=0.1, orange if ES>0.5 green if both)
with(subset(res, pval<=0.1 ), points(ES, -log10(pval), pch=20, col="red"))
with(subset(res, abs(ES)>0.5), points(ES, -log10(pval), pch=20, col="orange"))
with(subset(res, pval<0.1 & abs(ES)>0.5), points(ES, -log10(pval), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<=0.1 & abs(ES)>0.1), textxy(ES, -log10(pval), labs=pathway, cex=.8))