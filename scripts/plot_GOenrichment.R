
# GO enrichment
#####

library(gprofiler2)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)

df_names <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/DGElist_withNames.csv")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"

#####

data <- df_names %>% 
  mutate(Expression = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))

up <- data %>%filter(Expression == 'Up-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))

down <- data %>%filter(Expression == 'Down-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))

# query

multi_gp <- gost(list("Upregulated" = up$gene_name, "Downregulated" = down$gene_name), 
                 correction_method = "fdr", multi_query = F, ordered_query = T, 
                 organism = 'hsapiens')
gost_query <- as.data.frame(multi_gp$result)
colnames(gost_query)
class(multi_gp)
save(multi_gp, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/GOST_result.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/GOST_result.RData")
# manhattan plot

gostp1 <- gostplot(multi_gp, interactive = FALSE)
ggsave(paste0(outdir,"manhattanGO.png"),
       plot = gostp1, dpi = 300)
dev.off()

#####

# barplots

bar_data <- data.frame("term" = as.factor(gost_query$term_name), "condition" = gost_query$query, "count" = gost_query$term_size, "p.adjust" = gost_query$p_value, 'category' = as.factor(gost_query$source))

top_terms <- head(bar_data[order(bar_data$p.adjust),],40)
top_terms <- subset(top_terms, p.adjust < 1e-20)
write.csv(top_terms, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/Gprofiler_TopTerms.csv")

bar_data_up <- subset(bar_data, condition == 'Upregulated')
bar_data_up <- head(bar_data_up[order(bar_data_up$p.adjust),],15)

bar_data_down <- subset(bar_data, condition == 'Downregulated')
bar_data_down <- head(bar_data_down[order(bar_data_down$p.adjust),],15)

bar_data_reduced <- rbind(bar_data_up, bar_data_down)

bar_data_ordered <- bar_data_reduced[order(bar_data_reduced$p.adjust),] # order by count
bar_data_ordered<- bar_data_ordered[order(bar_data_ordered$category),] # order by category
bar_data_ordered$num <- seq(1:nrow(bar_data_ordered)) # num category for plot

g <- ggplot(bar_data_ordered, aes(count, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = count),
    color = "black",
    hjust = -0.1,
    size = 4,
    position = position_dodge(0.9)
  ) +
  labs(x = "Gene counts" , y = NULL) +
  scale_fill_manual(name='Category', labels = c('Biological Process', 'Cellular Component',
                                                'Molecular Function', 'REAC'), values = c('#3C6997', '#DD7230','#B4DC7F','#25ced1')) +
  theme(
    legend.position = "right",
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.background = element_blank()
  )

ggsave(paste0(outdir,"barplotGO.png"),
       plot = g, dpi = 300, width = 1000, height = 800, units = 'px')
dev.off()



## barplot only upregulated

bar_data_up <- subset(bar_data, condition == 'Upregulated')
bar_data_up <-head(bar_data_up[order(bar_data_up$p.adjust),],40) # order by pvalue
bar_data_up_ordered <- bar_data_up[order(bar_data_up$p.adjust),] # order by pvalue
bar_data_up_ordered<- bar_data_up_ordered[order(bar_data_up_ordered$category),] # order by category
bar_data_up_ordered$p.val <- round(-log10(bar_data_up_ordered$p.adjust), 2)
bar_data_up_ordered$num <- seq(1:nrow(bar_data_up_ordered)) # num category for plot

g.up <- ggplot(bar_data_up_ordered, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = p.val),
    color = "black",
    hjust = 0,
    size = 2.2,
    position = position_dodge(0)
  ) +
  labs(x = "-log10(p-value)" , y = NULL) +
  scale_fill_manual(name='Category', labels = c('Biological Process', 'REAC', 'TF'), values = c('#3C6997', '#DD7230','#B4DC7F','#25ced1')) +
  theme(
    legend.position = "right",
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.background = element_blank() 
  ) + theme_classic()

ggsave(paste0(outdir,"barplotUP_GO.png"),
       plot = g.up + theme_classic(), dpi = 300, width = 20, height = 10)
dev.off()


save(g.up, file = paste0(outdir, "barplotGO_UP.RData"))


## barplot only downregulated

bar_data_down <- subset(bar_data, condition == 'Downregulated')
bar_data_down <-head(bar_data_down[order(bar_data_down$p.adjust),],40) # order by pvalue
bar_data_down_ordered <- bar_data_down[order(bar_data_down$p.adjust),] # order by pvalue
bar_data_down_ordered<- bar_data_down_ordered[order(bar_data_down_ordered$category),] # order by category
bar_data_down_ordered$p.val <- round(-log10(bar_data_down_ordered$p.adjust), 2)
bar_data_down_ordered$num <- seq(1:nrow(bar_data_down_ordered)) # num category for plot

g.down <- ggplot(bar_data_down_ordered, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = p.val),
    color = "black",
    hjust = 0,
    size = 2.2,
    position = position_dodge(0)
  ) +
  labs(x = "-log10(p-value)" , y = NULL) +
  scale_fill_manual(name='Category', labels = c('CORUM', 'Biological Process', 'Cellular Component', 'Molecular Function'), values = c('#e810cb', '#3C6997','#d9c621','#b30039')) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )+ theme_classic()

ggsave(paste0(outdir,"barplotDOWN_GO.png"),
  plot = g.down + theme_classic(), dpi = 300, width = 20, height = 10)
# dev.off()

save(g.down, file = paste0(outdir, "barplotGO_DOWN.RData"))


#####