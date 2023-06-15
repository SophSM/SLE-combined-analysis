# Barplots for all DEGs
#####
library(tidyverse)
library(ggplot2)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"
DGE.list <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/DGElist_withNames_noncoding.csv")
DGE.list <- DGE.list[,-1]
#####


upreg <- na.omit(subset(DGE.list, padj < 0.05 & (log2FoldChange > 1)))
downreg <- na.omit(subset(DGE.list, padj < 0.05 & (log2FoldChange < -1 )))

levels(as.factor(upreg$transcript_biotype))

count.types1 <- as.data.frame(table(upreg$transcript_biotype))
keep1 <- c("lncRNA","miRNA","protein_coding","protein_coding_CDS_not_defined","snoRNA","snRNA")
df1.small <-count.types1[count.types1$Var1 %in% keep1,]


levels(as.factor(downreg$transcript_biotype))
count.types2 <- as.data.frame(table(downreg$transcript_biotype))
keep2 <- c("lncRNA","protein_coding","protein_coding_CDS_not_defined")
df2.small <-count.types2[count.types2$Var1 %in% keep2,]


## barplots

p <- ggplot(data = df1.small, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', fill = '#c21f35') + 
  geom_text(aes(label = Freq, x = Var1, y = Freq, hjust = -0.5, size = 2)) +
  coord_flip() +
  xlab(NULL) + 
  theme_classic() +
  theme(legend.position="none",
        axis.ticks.y=element_blank())

outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"
save(p, file=paste0(outdir,"barplot_up_nc.RData"))
# ggsave(paste0(outdir, "biotypes_barplot.png"), plot = p, dpi = 300)

p2 <- ggplot(data = df2.small, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', fill = '#195ae6') + 
  geom_text(aes(label = Freq, x = Var1, y = Freq, hjust = -0.5, size = 2)) +
  coord_flip() +
  xlab(NULL) +
  theme_classic() +
  theme(legend.position="none",
        axis.ticks.y=element_blank())
# ggsave(paste0(outdir, "biotypes_barplot_2.png"), plot = p2, dpi = 300)

save(p2, file=paste0(outdir,"barplot_down_nc.RData"))
