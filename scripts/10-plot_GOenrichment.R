
# GO enrichment
#####

library(gprofiler2)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
DIR = "/Users/sofiasalazar/Desktop/SLE-combined-analysis"
FIGDIR = glue::glue("{DIR}/figures_updated")
DGElist = data.table::fread(glue::glue("{DIR}/out/DGElist_all.csv"))
metadata <- read.table(glue::glue("{DIR}/data/all_data.csv"), header = T, row.names = 1,
                       sep = ",")

#####

DGElist_proteincoding <- DGElist %>% 
  filter(transcript_biotype == "protein_coding") %>%
  mutate(Expression = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))

up <- DGElist_proteincoding %>%filter(Expression == 'Up-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))

down <- DGElist_proteincoding %>%filter(Expression == 'Down-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))

# query

multi_gp <- gost(list("Up-regulated" = up$gene_name, "Down-regulated" = down$gene_name), 
                 correction_method = "fdr", multi_query = F, ordered_query = T, 
                 organism = 'hsapiens')
gost_query <- as.data.frame(multi_gp$result)
colnames(gost_query)
class(multi_gp)


#####

# barplots

bar_data <- data.frame("term" = as.factor(gost_query$term_name),
                       "direction" = gost_query$query, 
                       "intersection_size" = gost_query$intersection_size, 
                       "p.adjust" = gost_query$p_value, 
                       'category' = as.factor(gost_query$source))


bar_data_reduced <- bar_data %>%
  group_by(direction, category) %>%
  slice_max(order_by = p.adjust, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  filter(category %in% c("GO:BP", "REAC"))
  
bar_data_table <-rbind( bar_data_reduced ,
                          bar_data %>%
                          filter(direction == "Down-regulated") %>%
                          arrange(p.adjust) %>% 
                          head(15)
                          )
# write.table(bar_data_table, file = glue::glue("{DIR}/out/protein_coding_go_terms.csv"),
            # sep = ",", row.names = F, quote = F)
g <- ggplot(bar_data_reduced, 
            aes(
                x = direction, 
                y = reorder(term, -p.adjust),
                size = -log10(p.adjust),
                color = category)) +
  geom_point() +
  # geom_text(
  #   aes(label = count),
  #   color = "black",
  #   hjust = -0.1,
  #   size = 4,
  #   position = position_dodge(0.9)
  # ) +
  labs(x = "" , y = "") +
  scale_color_manual(name='Category', labels = c('Biological Process', 'REAC','Cellular Component',
                                                'Molecular Function', 'REAC', 'TF'), 
                    values = c('#3C6997','#DD7230', '#B4DC7F','#25ced1', '#8237de')) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "right",
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.background = element_blank()
  )

png(filename = glue::glue("{FIGDIR}/DEG_protein_coding_GOterms.png"),
    height = 18, width = 25, units = "cm", res = 500)
print(g)
dev.off()
ggsave(paste0(outdir,"barplotGO.png"),
       plot = g, dpi = 300, width = 1000, height = 800, units = 'px')
dev.off()


bar_data_small_terms<- data.table::fread(glue::glue("{DIR}/out/protein_coding_go_terms.csv"))
bar_data_small_terms <- bar_data_small_terms %>%
  filter(exclude == "no")

mat_small_meta <- bar_data_small_terms %>%
  dplyr::select(c(small_term, direction, category, p.adjust)) %>%
  unique()

mat_small <- mat_small_meta %>%
  mutate(logp = -log10(p.adjust)) %>%
  dplyr::select(c(small_term, direction, logp)) %>%
  pivot_wider(names_from = direction, values_from = logp) %>%
  column_to_rownames("small_term") %>%
  as.matrix()

col_pval <- colorRamp2(c(0, 1, 1.4, 1.6, 1.8, 2, 2.2),
                       c("white", "#FFFFCC", "#C7E9B4",
                         "#7FCDBB", "#41B6C4", "#2C7FB8", "#253494"))

category_annot <- rowAnnotation("Category" = mat_small_meta$category,
                                col = list("Category" = c("GO:BP" = "#a61b52",
                                                          "REAC" = "#004d01",
                                                          "GO:CC" = "#25ced1",
                                                          "GO:MF" = "#DD7230")),
                                gp = gpar(col = "white"))
ht_terms <- Heatmap(mat_small, name = "-log10(p-value)", cluster_rows = F,
        cluster_columns = F, row_names_side = "left", rect_gp = gpar(col = "white"),
        col = col_pval, left_annotation = category_annot,
        row_names_max_width = grid::unit(12, "cm"))

png(filename = glue::glue("{FIGDIR}/heatmap_protein_coding_terms.png"), height = 22,
    width = 14, units = "cm", res = 500)
draw(ht_terms)
dev.off()


category_annot_h <- HeatmapAnnotation("Category" = mat_small_meta$category,
                                col = list("Category" = c("GO:BP" = "#a61b52",
                                                          "REAC" = "#004d01",
                                                          "GO:CC" = "#25ced1",
                                                          "GO:MF" = "#DD7230")),
                                gp = gpar(col = "white"))

ht_terms_horizontal <- Heatmap(t(mat_small), name = "-log10(p-value)", cluster_rows = F,
                    cluster_columns = F, row_names_side = "left", rect_gp = gpar(col = "white"),
                    col = col_pval, top_annotation = category_annot_h,
                    column_names_rot = 45,
                    # row_names_max_width = grid::unit(25, "cm"),
                    column_names_max_height = grid::unit(12, "cm"),
                    column_names_gp = gpar(fontsize = 9))
png(filename = glue::glue("{FIGDIR}/heatmap_horizontal_protein_coding_terms.png"),
    height = 7,
    width = 30, units = "cm", res = 500)
draw(ht_terms_horizontal)
dev.off()
## reduce terms with rvgo
library(rrvgo)

upTerms <- gost_query %>% filter(query == 'Up-regulated' & source == "GO:BP")
downTerms <- gost_query %>% filter(query == 'Down-regulated' & source == "GO:BP")


simMat_up <- calculateSimMatrix(upTerms$term_id, # vector GO terms
                             orgdb="org.Hs.eg.db",
                             ont="BP", 
                             method="Rel")
scores_up <- setNames(-log10(upTerms$p_value), upTerms$term_id)
reducedTerms_up <- reduceSimMatrix(simMat_up,
                                   scores_up,
                                threshold= 0.7,
                                orgdb="org.Hs.eg.db")

# simMat_down <- calculateSimMatrix(downTerms$term_id, # vector GO terms
#                                 orgdb="org.Hs.eg.db",
#                                 ont="BP", 
#                                 method="Rel")
# scores_down <- setNames(-log10(downTerms$p_value), downTerms$term_id)
# reducedTerms_down <- reduceSimMatrix(simMat_down,
#                                      scores_down,
#                                    threshold= 0.7,
#                                    orgdb="org.Hs.eg.db")
# 
# reducedTerms_down <- reducedTerms_down[reducedTerms_down$size != 0,]
# reducedTerms_up <- reducedTerms_up[reducedTerms_up$size != 0,]

# get most significant term per parent term
# reducedTerms_downDF <- reducedTerms_down %>%
#   group_by(parentTerm) %>%
#   summarise(max_logpval = max(score))%>%
#   arrange(max_logpval)%>%
#   mutate(direction= "Downregulated")

reducedTerms_upDF <- reducedTerms_up %>%
  group_by(parentTerm) %>%
  summarise(max_logpval = max(score))%>%
  arrange(desc(max_logpval))%>%
  head(20) %>%
  arrange(max_logpval)%>%
  mutate(direction= "Up-regulated")

parentTerms <- (reducedTerms_upDF)

parentTerms_df <- parentTerms %>%
  distinct(parentTerm, .keep_all = TRUE)
parentTerms_df$parentTerm <- factor(parentTerms_df$parentTerm, levels = parentTerms_df$parentTerm)
parentTerms_df <- rbind(parentTerms_df,
                        data.frame(parentTerm = downTerms$term_name,
                                   max_logpval = -log10(downTerms$p_value),
                                   direction  = "Down-regulated"))
p_parent <- ggplot(parentTerms_df, aes(x = direction, size = max_logpval, y = parentTerm, color = direction))+
  geom_point() +
  scale_color_manual(values = c( "Up-regulated" = "firebrick3", "Down-regulated" = "dodgerblue3"))+
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12, hjust = 1, vjust = 1, angle = 45 ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.5, "lines")) +
  labs(y = "Biological process", size = expression("-log"[10]*"Adj. p-value"),
       color = "Direction", x = "") 

png(filename = glue::glue("{FIGDIR}/parentTerms_dotplot.png"),
       width = 25, height = 20, units = "cm", res = 500)
print(p_parent)
dev.off()


#####
sessionInfo()

# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Chicago
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] ggplotify_0.1.3             patchwork_1.3.2             RColorBrewer_1.1-3         
# [4] biomaRt_2.64.0              DESeq2_1.48.2               SummarizedExperiment_1.38.1
# [7] Biobase_2.68.0              MatrixGenerics_1.20.0       matrixStats_1.5.0          
# [10] GenomicRanges_1.60.0        GenomeInfoDb_1.44.3         IRanges_2.42.0             
# [13] S4Vectors_0.46.0            BiocGenerics_0.54.1         generics_0.1.4             
# [16] shiny_1.13.0                rrvgo_1.20.0                ggvenn_0.1.19              
# [19] circlize_0.4.17             gprofiler2_0.2.4            ComplexHeatmap_2.25.2      
# [22] lubridate_1.9.5             forcats_1.0.1               stringr_1.6.0              
# [25] dplyr_1.2.0                 purrr_1.2.1                 readr_2.2.0                
# [28] tidyr_1.3.2                 tibble_3.3.1                tidyverse_2.0.0            
# [31] ggplot2_4.0.2              
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.5.0           later_1.4.6             bitops_1.0-9           
# [4] filelock_1.0.3          R.oo_1.27.1             lifecycle_1.0.5        
# [7] httr2_1.2.2             doParallel_1.0.17       NLP_0.3-2              
# [10] lattice_0.22-9          MASS_7.3-65             magrittr_2.0.4         
# [13] plotly_4.12.0           rmarkdown_2.30          yaml_2.3.12            
# [16] httpuv_1.6.16           otel_0.2.0              askpass_1.2.1          
# [19] reticulate_1.45.0       cowplot_1.2.0           DBI_1.2.3              
# [22] abind_1.4-8             pkgload_1.5.0           R.utils_2.13.0         
# [25] RCurl_1.98-1.17         yulab.utils_0.2.4       rappdirs_0.3.4         
# [28] GenomeInfoDbData_1.2.14 tm_0.7-16               ggrepel_0.9.6          
# [31] pheatmap_1.0.13         umap_0.2.10.0           RSpectra_0.16-2        
# [34] codetools_0.2-20        DelayedArray_0.34.1     DOSE_4.2.0             
# [37] xml2_1.5.2              tidyselect_1.2.1        shape_1.4.6.1          
# [40] UCSC.utils_1.4.0        farver_2.1.2            BiocFileCache_2.16.2   
# [43] jsonlite_2.0.0          GetoptLong_1.1.0        iterators_1.0.14       
# [46] bbmle_1.0.25.1          foreach_1.5.2           tools_4.5.0            
# [49] progress_1.2.3          Rcpp_1.1.1              glue_1.8.0             
# [52] SparseArray_1.8.1       xfun_0.56               qvalue_2.40.0          
# [55] withr_3.0.2             numDeriv_2016.8-1.1     BiocManager_1.30.27    
# [58] fastmap_1.2.0           openssl_2.3.4           digest_0.6.39          
# [61] timechange_0.4.0        R6_2.6.1                mime_0.13              
# [64] gridGraphics_0.5-1      colorspace_2.1-2        GO.db_3.21.0           
# [67] dichromat_2.0-0.1       RSQLite_2.4.6           R.methodsS3_1.8.2      
# [70] utf8_1.2.6              data.table_1.18.2.1     prettyunits_1.2.0      
# [73] httr_1.4.8              htmlwidgets_1.6.4       S4Arrays_1.8.1         
# [76] pkgconfig_2.0.3         gtable_0.3.6            blob_1.3.0             
# [79] S7_0.2.1                XVector_0.48.0          htmltools_0.5.9        
# [82] fgsea_1.34.2            clue_0.3-67             scales_1.4.0           
# [85] png_0.1-8               wordcloud_2.6           knitr_1.51             
# [88] rstudioapi_0.18.0       reshape2_1.4.5          tzdb_0.5.0             
# [91] rjson_0.2.23            coda_0.19-4.1           curl_7.0.0             
# [94] bdsmatrix_1.3-7         org.Hs.eg.db_3.21.0     cachem_1.1.0           
# [97] GlobalOptions_0.1.3     parallel_4.5.0          AnnotationDbi_1.70.0   
# [100] treemap_2.4-4           apeglm_1.30.0           pillar_1.11.1          
# [103] vctrs_0.7.1             slam_0.1-55             promises_1.5.0         
# [106] dbplyr_2.5.2            xtable_1.8-4            cluster_2.1.8.2        
# [109] evaluate_1.0.5          mvtnorm_1.3-3           cli_3.6.5              
# [112] locfit_1.5-9.12         compiler_4.5.0          rlang_1.1.7            
# [115] crayon_1.5.3            labeling_0.4.3          emdbook_1.3.14         
# [118] plyr_1.8.9              fs_1.6.6                stringi_1.8.7          
# [121] viridisLite_0.4.3       gridBase_0.4-7          BiocParallel_1.42.2    
# [124] Biostrings_2.76.0       lazyeval_0.2.2          GOSemSim_2.34.0        
# [127] Matrix_1.7-4            hms_1.1.4               bit64_4.6.0-1          
# [130] KEGGREST_1.48.1         igraph_2.2.2            memoise_2.0.1          
# [133] fastmatch_1.1-8         bit_4.6.0   