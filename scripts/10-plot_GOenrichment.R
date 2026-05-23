
# GO enrichment
#####

library(gprofiler2)
# library(enrichplot)
library(DOSE)
# library(clusterProfiler)
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
  arrange(desc(abs(log2FoldChange)))

down <- DGElist_proteincoding %>%filter(Expression == 'Down-regulated') %>% 
  arrange(desc(abs(log2FoldChange)))


# rank gene lists

go_query_lst <- list("Up-regulated" = up$gene_name,
                     "Down-regulated" = down$gene_name)
# query

multi_gp <- gost(go_query_lst, 
                 correction_method = "fdr", multi_query = T, ordered_query = T, 
                 organism = 'hsapiens', sources = c("GO:BP", "GO:MF","KEGG", "REAC"))
# gost_query <- as.data.frame(multi_gp$result)
# colnames(gost_query)
# class(multi_gp)

go_results <- multi_gp$result %>%
  dplyr::select(c(term_name, source, p_values, significant, intersection_sizes)) %>% 
  mutate(direction = list(names(go_query_lst))) %>%
  unnest(direction, p_values, significant, intersection_sizes) %>%
  filter(significant == T)
write.table(go_results, file = glue::glue("{DIR}/out/DEG_full_go_terms.csv"),
            sep = ",", row.names = F, quote = F)
#####

# barplots


bar_data_reduced <- go_results %>%
  group_by(direction, source) %>%
  slice_min(order_by = p_values, n = 10, with_ties = FALSE) %>%
  ungroup()
  
# write.table(bar_data_reduced, file = glue::glue("{DIR}/out/DEG_reduced_go_terms_2.csv"),
            # sep = ",", row.names = F, quote = F)
g <- ggplot(bar_data_reduced, 
            aes(x = direction, 
                y = reorder(term_name, -p_values),
                size = -log10(p_values),
                color = source)) +
  geom_point() +
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



bar_data_small_terms<- read.table(glue::glue("{DIR}/out/DEG_reduced_go_terms_2.csv"),
                                         sep = ",", header = T)
bar_data_small_terms <- bar_data_small_terms %>%
  na.omit()

mat_small_meta <- bar_data_small_terms %>%
  dplyr::select(c(small_term, direction, source, p_values)) %>%
  unique()

mat_small <- mat_small_meta %>%
  mutate(logp = -log10(p_values)) %>%
  dplyr::select(c(small_term, direction, logp)) %>%
  group_by(direction, small_term) %>%
  slice_max(logp, n = 1) %>%
  ungroup() %>%
  pivot_wider(names_from = direction, values_from = logp) %>%
  column_to_rownames("small_term") %>%
  as.matrix()

mat_small_meta_2 <- mat_small_meta %>%
  slice(match(rownames(mat_small) , small_term))

col_pval <- colorRamp2(c(0, 0.5, 1, 3, 6, 8, 16, 30, 40),
                       c("white", "#FFFFD9", "#C7E9B4", "#7FCDBB", "#41B6C4",
                         "#1D91C0","#225EA8", "#253494", "#081D58"))

category_annot <- rowAnnotation("Category" = mat_small_meta_2$source,
                                col = list("Category" = c("GO:BP" = "#a61b52",
                                                          "GO:MF" = "#DD7230",
                                                          "KEGG" = '#004d01',
                                                          "REAC" = '#8237de' )),
                                gp = gpar(col = "white"))
ht_terms <- Heatmap(mat_small, name = "-log10(p-value)", cluster_rows = F,
        cluster_columns = F, row_names_side = "left", rect_gp = gpar(col = "white"),
        col = col_pval, left_annotation = category_annot,
        row_names_max_width = grid::unit(12, "cm"))

png(filename = glue::glue("{FIGDIR}/heatmap_protein_coding_terms.png"), height = 22,
    width = 14, units = "cm", res = 500)
draw(ht_terms)
dev.off()

pdf(glue::glue("{FIGDIR}/heatmap_protein_coding_terms.pdf"), height = (15 / 2.54),
    width = (14/2.54))
draw(ht_terms)
dev.off()
category_annot_h <- HeatmapAnnotation("Category" = mat_small_meta_2$source,
                                col = list("Category" = c("GO:BP" = "#a61b52",
                                                          "GO:MF" = "#DD7230",
                                                          "KEGG" = '#004d01',
                                                          "REAC" = '#8237de' )),
                                gp = gpar(col = "white"))

ht_terms_horizontal <- Heatmap(t(mat_small), name = "-log10(p-value)", cluster_rows = F,
                    cluster_columns = F, row_names_side = "left", rect_gp = gpar(col = "white"),
                    col = col_pval, top_annotation = category_annot_h,
                    column_names_rot = 45,
                    na_col = "lightgray",
                    # row_names_max_width = grid::unit(25, "cm"),
                    column_names_max_height = grid::unit(15, "cm"),
                    column_names_gp = gpar(fontsize = 10))
png(filename = glue::glue("{FIGDIR}/heatmap_horizontal_protein_coding_terms.png"),
    height = 8,
    width = 30, units = "cm", res = 500)
draw(ht_terms_horizontal)
dev.off()

pdf(glue::glue("{FIGDIR}/heatmap_horizontal_protein_coding_terms.pdf"),
    height = (8 / 2.54),
    width = (30 / 2.54))
draw(ht_terms_horizontal)
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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DOSE_4.2.0            gprofiler2_0.2.4      ggpubr_0.6.2         
# [4] ggpattern_1.2.1       patchwork_1.3.2       circlize_0.4.17      
# [7] ComplexHeatmap_2.25.2 lubridate_1.9.5       forcats_1.0.1        
# [10] stringr_1.6.0         dplyr_1.2.0           purrr_1.2.1          
# [13] readr_2.2.0           tidyr_1.3.2           tibble_3.3.1         
# [16] ggplot2_4.0.2         tidyverse_2.0.0      
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-9            DBI_1.2.3               rlang_1.1.7            
# [4] magrittr_2.0.4          clue_0.3-67             GetoptLong_1.1.0       
# [7] otel_0.2.0              matrixStats_1.5.0       compiler_4.5.0         
# [10] RSQLite_2.4.6           reshape2_1.4.5          png_0.1-8              
# [13] vctrs_0.7.1             pkgconfig_2.0.3         shape_1.4.6.1          
# [16] crayon_1.5.3            fastmap_1.2.0           backports_1.5.0        
# [19] XVector_0.48.0          labeling_0.4.3          utf8_1.2.6             
# [22] tzdb_0.5.0              UCSC.utils_1.4.0        bit_4.6.0              
# [25] xfun_0.56               cachem_1.1.0            GenomeInfoDb_1.44.3    
# [28] jsonlite_2.0.0          blob_1.3.0              BiocParallel_1.42.2    
# [31] broom_1.0.12            parallel_4.5.0          cluster_2.1.8.2        
# [34] R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3     
# [37] car_3.1-5               GOSemSim_2.34.0         Rcpp_1.1.1             
# [40] iterators_1.0.14        knitr_1.51              R.utils_2.13.0         
# [43] IRanges_2.42.0          splines_4.5.0           Matrix_1.7-4           
# [46] timechange_0.4.0        tidyselect_1.2.1        qvalue_2.40.0          
# [49] rstudioapi_0.18.0       dichromat_2.0-0.1       abind_1.4-8            
# [52] doParallel_1.0.17       codetools_0.2-20        curl_7.0.0             
# [55] plyr_1.8.9              lattice_0.22-9          Biobase_2.68.0         
# [58] withr_3.0.2             KEGGREST_1.48.1         S7_0.2.1               
# [61] evaluate_1.0.5          Biostrings_2.76.0       pillar_1.11.1          
# [64] carData_3.0-6           foreach_1.5.2           stats4_4.5.0           
# [67] plotly_4.12.0           generics_0.1.4          RCurl_1.98-1.17        
# [70] S4Vectors_0.46.0        hms_1.1.4               scales_1.4.0           
# [73] glue_1.8.0              lazyeval_0.2.2          tools_4.5.0            
# [76] data.table_1.18.2.1     fgsea_1.34.2            ggsignif_0.6.4         
# [79] fs_1.6.6                fastmatch_1.1-8         cowplot_1.2.0          
# [82] AnnotationDbi_1.70.0    colorspace_2.1-2        GenomeInfoDbData_1.2.14
# [85] Formula_1.2-5           cli_3.6.5               rappdirs_0.3.4         
# [88] viridisLite_0.4.3       gtable_0.3.6            yulab.utils_0.2.4      
# [91] R.methodsS3_1.8.2       rstatix_0.7.3           digest_0.6.39          
# [94] BiocGenerics_0.54.1     rjson_0.2.23            htmlwidgets_1.6.4      
# [97] farver_2.1.2            R.oo_1.27.1             memoise_2.0.1          
# [100] htmltools_0.5.9         lifecycle_1.0.5         httr_1.4.8             
# [103] GlobalOptions_0.1.3     GO.db_3.21.0            bit64_4.6.0-1   