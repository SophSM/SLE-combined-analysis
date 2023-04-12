# Volcano plot

#####

library(ggplot2)

load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"

#####

upGenes <- subset(df_names, padj < 0.05 & log2FoldChange >= 1)
veryUP <- subset(df_names, padj < 0.05 & log2FoldChange >= 2)

downGenes <- subset(df_names, padj < 0.05 & log2FoldChange <= -0.5)
veryDOWN <- subset(df_names, padj < 0.05 & log2FoldChange <= -2)

data <- df_names %>% 
  mutate(Expression = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))
head(data) %>% 
  knitr::kable()

volcanoplot <- ggplot(data, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-adj")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

data %>% 
  count(Expression) %>% 
  knitr::kable()


top <- 10 # no. of highlighted genes in plot
top_genes <- bind_rows(data %>% 
                         filter(Expression == 'Up-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top),
                       data %>% 
                         filter(Expression == 'Down-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top)
)
top_genes %>% 
  knitr::kable()

volcanoplot_names <-  volcanoplot +
  ggrepel::geom_label_repel(data = top_genes,
                            mapping = aes(log2FoldChange, -log(padj,10), label = gene_name),
                            size = 2)

ggsave(paste0(outdir,"volcanoPlotWithTopGenes.png"),
       plot = p3, dpi = 300)
dev.off()