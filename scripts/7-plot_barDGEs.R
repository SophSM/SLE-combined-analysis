# Barplots for all DEGs
#####
library(tidyverse)
library(plyr)
library(ggplot2)
DIR = "/Users/sofiasalazar/Desktop/SLE-combined-analysis"
FIGDIR = glue::glue("{DIR}/figures_updated")

DGElist = data.table::fread(glue::glue("{DIR}/out/DGElist_all.csv"))
#####
DGElist <- DGElist %>%
  mutate(expression = case_when(
    log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
    TRUE ~ "Unchanged"))


# count.types1 <- as.data.frame(table(upreg$transcript_biotype))
keep_biotypes <- c("lncRNA","miRNA","protein_coding","protein_coding_CDS_not_defined","snoRNA","snRNA")

summary_biotypes <- DGElist %>%
  filter(transcript_biotype %in% keep_biotypes) %>%
  group_by(expression, transcript_biotype) %>%
  summarise(count = n_distinct(ID))

## barplots

p <- ggplot(data = summary_biotypes %>%
              filter(expression == "Up-regulated"),
            aes(x = count, y = transcript_biotype)) + 
  geom_bar(stat = 'identity', fill = "firebrick3") + 
  # scale_fill_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  scale_x_continuous(limits = c(0,306), expand = c(0, 0)) + 
  geom_text(aes(label = count, x = count, hjust = -0.5, size = 2)) +
  labs(x = "Number of transcripts", y = "", title = "Up-regulated transcripts") +
  theme_minimal() +
  theme(legend.position="none",
        axis.ticks.y=element_blank())

save(p, file=glue::glue("{FIGDIR}/barplot_up_biotypes.RData"))
png(filename = glue::glue("{FIGDIR}/barplot_up_biotypes.png"),
    height = 10, width = 15, units = "cm", res = 500)
print(p)
dev.off()


p2 <- ggplot(data = summary_biotypes %>%
               filter(expression == "Down-regulated"),
             aes(x = count, y = transcript_biotype)) + 
  geom_bar(stat = 'identity', fill = "dodgerblue") + 
  # scale_fill_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  scale_x_continuous(limits = c(0,35), expand = c(0, 0)) + 
  geom_text(aes(label = count, x = count, hjust = -0.5, size = 2)) +
  labs(x = "Number of transcripts", y = "", title = "Down-regulated transcripts") +
  theme_minimal() +
  theme(legend.position="none",
        axis.ticks.y=element_blank())

save(p2, file=glue::glue("{FIGDIR}/barplot_down_biotypes.RData"))
png(filename = glue::glue("{FIGDIR}/barplot_down_biotypes.png"),
    height = 10, width = 15, units = "cm", res = 500)
print(p2)
dev.off()

######
# sessionInfo()
# 
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
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] biomaRt_2.64.0              DESeq2_1.48.2              
# [3] SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [5] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# [7] GenomicRanges_1.60.0        GenomeInfoDb_1.44.3        
# [9] IRanges_2.42.0              S4Vectors_0.46.0           
# [11] BiocGenerics_0.54.1         generics_0.1.4             
# [13] shiny_1.13.0                rrvgo_1.20.0               
# [15] ggvenn_0.1.19               circlize_0.4.17            
# [17] gprofiler2_0.2.4            ComplexHeatmap_2.25.2      
# [19] lubridate_1.9.5             forcats_1.0.1              
# [21] stringr_1.6.0               dplyr_1.2.0                
# [23] purrr_1.2.1                 readr_2.2.0                
# [25] tidyr_1.3.2                 tibble_3.3.1               
# [27] tidyverse_2.0.0             ggplot2_4.0.2              
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.18.0       jsonlite_2.0.0         
# [4] shape_1.4.6.1           umap_0.2.10.0           magrittr_2.0.4         
# [7] farver_2.1.2            rmarkdown_2.30          GlobalOptions_0.1.3    
# [10] fs_1.6.6                vctrs_0.7.1             memoise_2.0.1          
# [13] askpass_1.2.1           progress_1.2.3          S4Arrays_1.8.1         
# [16] htmltools_0.5.9         curl_7.0.0              SparseArray_1.8.1      
# [19] htmlwidgets_1.6.4       httr2_1.2.2             plyr_1.8.9             
# [22] plotly_4.12.0           cachem_1.1.0            igraph_2.2.2           
# [25] mime_0.13               lifecycle_1.0.5         iterators_1.0.14       
# [28] pkgconfig_2.0.3         Matrix_1.7-4            R6_2.6.1               
# [31] fastmap_1.2.0           GenomeInfoDbData_1.2.14 clue_0.3-67            
# [34] numDeriv_2016.8-1.1     digest_0.6.39           colorspace_2.1-2       
# [37] AnnotationDbi_1.70.0    RSpectra_0.16-2         RSQLite_2.4.6          
# [40] org.Hs.eg.db_3.21.0     filelock_1.0.3          labeling_0.4.3         
# [43] timechange_0.4.0        abind_1.4-8             httr_1.4.8             
# [46] compiler_4.5.0          bit64_4.6.0-1           withr_3.0.2            
# [49] doParallel_1.0.17       S7_0.2.1                BiocParallel_1.42.2    
# [52] DBI_1.2.3               R.utils_2.13.0          MASS_7.3-65            
# [55] openssl_2.3.4           DelayedArray_0.34.1     rappdirs_0.3.4         
# [58] rjson_0.2.23            tools_4.5.0             otel_0.2.0             
# [61] httpuv_1.6.16           R.oo_1.27.1             glue_1.8.0             
# [64] GOSemSim_2.34.0         promises_1.5.0          gridBase_0.4-7         
# [67] cluster_2.1.8.2         gtable_0.3.6            tzdb_0.5.0             
# [70] R.methodsS3_1.8.2       data.table_1.18.2.1     hms_1.1.4              
# [73] xml2_1.5.2              utf8_1.2.6              XVector_0.48.0         
# [76] ggrepel_0.9.6           foreach_1.5.2           pillar_1.11.1          
# [79] emdbook_1.3.14          yulab.utils_0.2.4       later_1.4.6            
# [82] BiocFileCache_2.16.2    lattice_0.22-9          bit_4.6.0              
# [85] tidyselect_1.2.1        locfit_1.5-9.12         GO.db_3.21.0           
# [88] tm_0.7-16               Biostrings_2.76.0       knitr_1.51             
# [91] NLP_0.3-2               xfun_0.56               pheatmap_1.0.13        
# [94] stringi_1.8.7           UCSC.utils_1.4.0        lazyeval_0.2.2         
# [97] yaml_2.3.12             evaluate_1.0.5          codetools_0.2-20       
# [100] bbmle_1.0.25.1          wordcloud_2.6           BiocManager_1.30.27    
# [103] cli_3.6.5               xtable_1.8-4            reticulate_1.45.0      
# [106] treemap_2.4-4           dichromat_2.0-0.1       Rcpp_1.1.1             
# [109] dbplyr_2.5.2            coda_0.19-4.1           png_0.1-8              
# [112] bdsmatrix_1.3-7         parallel_4.5.0          blob_1.3.0             
# [115] prettyunits_1.2.0       mvtnorm_1.3-3           apeglm_1.30.0          
# [118] viridisLite_0.4.3       slam_0.1-55             scales_1.4.0           
# [121] crayon_1.5.3            GetoptLong_1.1.0        rlang_1.1.7            
# [124] KEGGREST_1.48.1  