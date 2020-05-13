res <- read.csv("NameMappedDE.csv", header=TRUE)
head(res)

# Make a basic volcano plot
with(res, plot(Log2FC, -log10(Pval), pch=20, main="Volcano plot", xlim=c(-8,8), ylim = c(0,150)))

# Add colored points: red if padj<=0.001, orange of log2FC>2, green if both)
with(subset(res, P_adj<=.001 ), points(Log2FC, -log10(Pval), pch=20, col="red"))
with(subset(res, abs(Log2FC)>2), points(Log2FC, -log10(Pval), pch=20, col="orange"))
with(subset(res, P_adj<.001 & abs(Log2FC)>2), points(Log2FC, -log10(Pval), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, P_adj<.0000000000001 & abs(Log2FC)>5), textxy(Log2FC, -log10(Pval), labs=Name, cex=.8))