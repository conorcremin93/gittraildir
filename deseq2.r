################


res1 <- results(dds, contrast=c("condition", "Wildtype", "Mock"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH")
res2 <- results(dds, contrast=c("condition", "delNS", "Mock"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH")
res3 <- results(dds, contrast=c("condition", "3M", "Mock"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH")
res4 <- results(dds, contrast=c("condition", "Wildtype", "delNS"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH")
res5 <- results(dds, contrast=c("condition", "Wildtype", "3M"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH")
res6<- results(dds, contrast=c("condition", "3M", "delNS"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH")

resLFC1 <- lfcShrink(dds, contrast=c("condition", "Wildtype", "Mock"), res=res1)
resLFC2 <- lfcShrink(dds, contrast=c("condition", "delNS", "Mock"), res=res2)
resLFC3 <- lfcShrink(dds, contrast=c("condition", "3M", "Mock"), res=res3)
resLFC4 <- lfcShrink(dds, contrast=c("condition", "Wildtype", "delNS"), res=res4)
resLFC5 <- lfcShrink(dds, contrast=c("condition", "Wildtype", "3M"), res=res5)
resLFC6 <- lfcShrink(dds, contrast=c("condition", "3M", "delNS"), res=res6)

resSig1 <- na.omit(resLFC1, cols = c("log2FoldChange", "padj"))
resSig2 <- na.omit(resLFC2, cols = c("log2FoldChange", "padj"))
resSig3 <- na.omit(resLFC3, cols = c("log2FoldChange", "padj"))
resSig4 <- na.omit(resLFC4, cols = c("log2FoldChange", "padj"))
resSig5 <- na.omit(resLFC5, cols = c("log2FoldChange", "padj"))
resSig6 <- na.omit(resLFC6, cols = c("log2FoldChange", "padj"))

Sig1 <- resSig1[(resSig1[,2] > 2) & resSig1[,6] < 0.05,]
Sig2 <- resSig2[(resSig2[,2] > 2) & resSig2[,6] < 0.05,]
Sig3 <- resSig3[(resSig3[,2] > 2) & resSig3[,6] < 0.05,]
Sig4 <- resSig4[(resSig4[,2] > 2) & resSig4[,6] < 0.05,]
Sig5 <- resSig5[(resSig5[,2] > 2) & resSig5[,6] < 0.05,]
Sig6 <- resSig6[(resSig6[,2] > 2) & resSig6[,6] < 0.05,]

reorder1 <- Sig1[order(Sig1[,2]),]
reorder2 <- Sig2[order(Sig2[,2]),]
reorder3 <- Sig3[order(Sig3[,2]),]
reorder4 <- Sig4[order(Sig4[,2]),]
reorder5 <- Sig5[order(Sig5[,2]),]
reorder6 <- Sig6[order(Sig6[,2]),]

write.table(as.data.frame(reorder1), quote = FALSE, sep="\t", file="Functional_Analysis/lncRNA_wildtype_vs_Mock.txt")
write.table(as.data.frame(reorder2), quote = FALSE, sep="\t", file="Functional_Analysis/lncRNA_delNS_vs_Mock.txt")
write.table(as.data.frame(reorder3), quote = FALSE, sep="\t", file="Functional_Analysis/lncRNA_3M_vs_Mock.txt")
write.table(as.data.frame(reorder4), quote = FALSE, sep="\t", file="Functional_Analysis/lncRNA_wildtype_vs_delNS.txt")
write.table(as.data.frame(reorder5), quote = FALSE, sep="\t", file="Functional_Analysis/lncRNA_wildtype_vs_3M.txt")
write.table(as.data.frame(reorder6), quote = FALSE, sep="\t", file="Functional_Analysis/lncRNA_3M_vs_delNS.txt")



d <- plotCounts(dds, gene=which.min(res1$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(dds, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
plotPCA(rld, intgroup = "condition")

###Distribution Plots
par(mfrow=c(3,1))
#Hist of raw counts
hist(as.matrix(countData), col="blue", border="white", breaks=100) 
#
hist(as.matrix(countData), col="blue", border="white",
     breaks=20000, xlim=c(0,2000), main="Counts per gene",
     xlab="Counts (truncated axis)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
#
epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(countData + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

par(mfrow=c(1,1))
##Box plots of non-normalized log2(counts+1) per sample.

boxplot(log2(countData + epsilon), col=colData$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(counts +1)")
###

# We will require one function from the affy package
a

###

###
###