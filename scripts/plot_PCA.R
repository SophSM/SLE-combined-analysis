# PCA
#####
library(DESeq2)
library(ggplot2)
load('/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData')
load("/mnt/Citosina/amedina/ssalazar/meta/combined/vsd2.RData")
outdir = '/mnt/Citosina/amedina/ssalazar/meta/combined/figures/'
######
mat <- as.matrix(assay(vsd2))
pc <- prcomp(t(mat))
dtp <- data.frame('DISEASE' = all_data$DISEASE, pc$x[,1:2])

pcaData <- plotPCA(vsd2, 'DISEASE', returnData = TRUE)

pcaData <- tibble::rownames_to_column(pcaData, "study")

percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=group, group = group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  stat_ellipse() + theme_classic() +
  labs(fill = 'Samples')

save(pca_plot, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/pca_object.RData")

ggsave(paste0(outdir,"PCA-ellipse.png"), width = 3000, height = 3000,
       units = 'px', dpi = 300, bg = "white", plot = pca_plot)
dev.off()
