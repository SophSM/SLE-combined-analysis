# Barplots for all DEGs
#####
library(tidyverse)
library(ggplot2)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"
DGE.list <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/DGElist_withNames_noncoding.csv")
DGE.list <- DGE.list[,-1]
#####


diffexpressed1 <- subset(DGE.list, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1 ))
diffexpressed2 <- subset(DGE.list, padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2 ))

levels(as.factor(diffexpressed1$transcript_biotype))
count.types1 <- as.data.frame(table(diffexpressed1$transcript_biotype))
keep1 <- c("lncRNA","miRNA","protein_coding","protein_coding_CDS_not_defined","snoRNA","snRNA")
df1.small <-count.types1[count.types1$Var1 %in% keep1,]

levels(as.factor(diffexpressed2$transcript_biotype))
count.types2 <- as.data.frame(table(diffexpressed2$transcript_biotype))
keep2 <- c("lncRNA","protein_coding","protein_coding_CDS_not_defined")
df2.small <-count.types2[count.types2$Var1 %in% keep2,]


## barplots

p <- ggplot(data = df1.small, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', fill = '#5d996c') + 
  geom_text(aes(label = Freq, vjust = -0.3, size = 3.5)) +
  coord_flip() +
  labs(title="Log2FoldChange < or > 1",x = NULL, y = "Frequency") +
  theme_classic() +
  theme(legend.position="none")

outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"
ggsave(paste0(outdir, "biotypes_barplot.png"), plot = p, dpi = 300)

p2 <- ggplot(data = df2.small, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', fill = '#a1506c') + 
  geom_text(aes(label = Freq, vjust = -0.3, size = 3.5)) +
  coord_flip() +
  labs(title="Log2FoldChange < or > 2",x = NULL, y = "Frequency") +
  theme_classic() +
  theme(legend.position="none")

ggsave(paste0(outdir, "biotypes_barplot_2.png"), plot = p2, dpi = 300)