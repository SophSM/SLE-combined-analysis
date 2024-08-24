
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

## reduce terms with rvgo
library(rrvgo)

upTerms <- gost_query %>% filter(query == 'Upregulated' & source == "GO:BP")
downTerms <- gost_query %>% filter(query == 'Downregulated' & source == "GO:BP")


simMat_up <- calculateSimMatrix(upTerms$term_id, # vector GO terms
                             orgdb="org.Hs.eg.db",
                             ont="BP", 
                             method="Rel")
scores_up <- setNames(-log10(upTerms$p_value), upTerms$term_id)
reducedTerms_up <- reduceSimMatrix(simMat_up,
                                   scores_up,
                                threshold= 0.7,
                                orgdb="org.Hs.eg.db")

simMat_down <- calculateSimMatrix(downTerms$term_id, # vector GO terms
                                orgdb="org.Hs.eg.db",
                                ont="BP", 
                                method="Rel")
scores_down <- setNames(-log10(downTerms$p_value), downTerms$term_id)
reducedTerms_down <- reduceSimMatrix(simMat_down,
                                     scores_down,
                                   threshold= 0.7,
                                   orgdb="org.Hs.eg.db")

reducedTerms_down <- reducedTerms_down[reducedTerms_down$size != 0,]
reducedTerms_up <- reducedTerms_up[reducedTerms_up$size != 0,]

# get most significant term per parent term
reducedTerms_downDF <- reducedTerms_down %>%
  group_by(parentTerm) %>%
  summarise(max_logpval = max(score))%>%
  arrange(max_logpval)%>%
  mutate(direction= "Downregulated")

reducedTerms_upDF <- reducedTerms_up %>%
  group_by(parentTerm) %>%
  summarise(max_logpval = max(score))%>%
  arrange(desc(max_logpval))%>%
  head(10)%>%
  arrange(max_logpval)%>%
  mutate(direction= "Upregulated")

parentTerms <- rbind(reducedTerms_downDF, reducedTerms_upDF)

parentTerms_df <- parentTerms %>%
  distinct(parentTerm, .keep_all = TRUE)
parentTerms_df$parentTerm <- factor(parentTerms_df$parentTerm, levels = parentTerms_df$parentTerm)

p_parent <- ggplot(parentTerms_df, aes(x = parentTerm,y = direction, size = max_logpval))+
  geom_point(aes(color = direction), position = position_dodge(width = 0.6), alpha = 0.7)+
  coord_flip() +
  theme_minimal()+
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.5, "lines"))+
  scale_size_continuous(range = c(3, 10)) + 
  scale_color_manual(values = c("Upregulated" = "firebrick3", "Downregulated"="dodgerblue"))+
  labs(x = "Biological process", y = "Expression", colour = "Expression", size = "-log10(p-value)")

ggsave(filename = glue::glue("{outdir}/parentTerms_dotplot.png"), plot = p_parent,
       width = 25, height = 20, units = "cm", dpi = 300)


# the scores column corresponds to the -log10(pvalue of the term)
#####
sessionInfo()

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS:   /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRblas.so
# LAPACK: /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRlapack.so

# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] forcats_0.5.1          stringr_1.4.0          dplyr_1.0.10          
# [4] purrr_0.3.4            readr_1.4.0            tidyr_1.2.1           
# [7] tibble_3.1.3           tidyverse_1.3.0        ggplot2_3.3.5         
# [10] clusterProfiler_3.18.1 DOSE_3.16.0            enrichplot_1.10.2     
# [13] gprofiler2_0.2.0      

# loaded via a namespace (and not attached):
#   [1] fs_1.5.0             lubridate_1.7.9.2    bit64_4.0.5         
# [4] RColorBrewer_1.1-2   httr_1.4.2           tools_4.0.2         
# [7] backports_1.2.1      utf8_1.2.2           R6_2.5.0            
# [10] DBI_1.1.1            lazyeval_0.2.2       BiocGenerics_0.36.1 
# [13] colorspace_2.0-2     withr_2.4.2          tidyselect_1.2.0    
# [16] gridExtra_2.3        bit_4.0.4            compiler_4.0.2      
# [19] rvest_0.3.6          cli_3.6.0            Biobase_2.50.0      
# [22] xml2_1.3.2           scatterpie_0.1.6     plotly_4.9.3        
# [25] shadowtext_0.0.8     scales_1.1.1         digest_0.6.27       
# [28] pkgconfig_2.0.3      htmltools_0.5.1.1    dbplyr_2.2.1        
# [31] fastmap_1.1.0        readxl_1.3.1         htmlwidgets_1.5.3   
# [34] rlang_1.0.6          rstudioapi_0.13      RSQLite_2.2.7       
# [37] farver_2.1.0         generics_0.1.0       jsonlite_1.7.2      
# [40] BiocParallel_1.24.1  GOSemSim_2.16.1      magrittr_2.0.1      
# [43] GO.db_3.12.1         Matrix_1.3-4         Rcpp_1.0.7          
# [46] munsell_0.5.0        S4Vectors_0.28.1     fansi_0.5.0         
# [49] viridis_0.5.1        lifecycle_1.0.3      stringi_1.6.2       
# [52] ggraph_2.0.5         MASS_7.3-53          plyr_1.8.6          
# [55] qvalue_2.22.0        grid_4.0.2           blob_1.2.1          
# [58] parallel_4.0.2       ggrepel_0.9.1        DO.db_2.9           
# [61] crayon_1.4.1         lattice_0.20-41      graphlayouts_0.7.1  
# [64] haven_2.3.1          cowplot_1.1.1        splines_4.0.2       
# [67] hms_1.1.0            pillar_1.6.2         fgsea_1.16.0        
# [70] igraph_1.2.8         reshape2_1.4.4       stats4_4.0.2        
# [73] fastmatch_1.1-0      reprex_1.0.0         glue_1.4.2          
# [76] downloader_0.4       modelr_0.1.8         data.table_1.14.0   
# [79] BiocManager_1.30.21  vctrs_0.5.1          tweenr_1.0.2        
# [82] cellranger_1.1.0     gtable_0.3.0         polyclip_1.10-0     
# [85] assertthat_0.2.1     cachem_1.0.5         ggforce_0.3.3       
# [88] broom_0.7.9          tidygraph_1.2.0      viridisLite_0.4.0   
# [91] rvcheck_0.1.8        AnnotationDbi_1.52.0 memoise_2.0.0       
# [94] IRanges_2.24.1       ellipsis_0.3.2      