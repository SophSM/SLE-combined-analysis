# PCA
#####
library(DESeq2)
library(ggplot2)
library(tidyverse)
DIR = "/Users/sofiasalazar/Desktop/SLE-combined-analysis"
FIGDIR = glue::glue("{DIR}/figures_updated")
# load(glue::glue("{DIR}/data/LRT-dds.RData"))
load(glue::glue("{DIR}/data/vsd2.RData"))
metadata <- read.table(glue::glue("{DIR}/data/all_data.csv"), header = T, row.names = 1,
                       sep = ",")
DGElist = data.table::fread(glue::glue("{DIR}/out/DGElist_all.csv"))
DGElist_prot <- DGElist %>%
  filter(transcript_biotype == "protein_coding")
######

mat <- as.matrix(assay(vsd2)) # col = samples, rows  = genes
rownames(mat) <- gsub("\\..*","", rownames(mat))
mat_prot <- mat[DGElist_prot$ID,]

pc <- prcomp(t(mat_prot), center = TRUE, scale. = TRUE) # col = genes , rows = samples
pca_df <- data.frame(pc$x[,1:2])
pca_df <- pca_df %>%
  rownames_to_column("sampleID") %>%
  full_join(metadata, by = c("sampleID" = "samples"))

var_explained <- pc$sdev^2 / sum(pc$sdev^2)

pc1_pct <- round(var_explained[1] * 100, 3)
pc2_pct <- round(var_explained[2] * 100, 3)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = DISEASE)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')) +
  labs(x = paste0("PC1: ",pc1_pct,"% variance"),
       y = paste0("PC2: ",pc2_pct,"% variance")) +
  stat_ellipse() +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),  
        axis.text.y = element_text(size = 25),
        legend.key.size = unit(1.5, "lines"))



save(pca_plot, file = glue::glue("{FIGDIR}/pca_object.RData"))
png(filename = glue::glue("{FIGDIR}/pca.png"), height = 20, width = 22,
    units = "cm", res = 500)
print(pca_plot)
dev.off()

ggsave(filename =glue::glue("{FIGDIR}/pca.svg"), height = 20, width = 22,
       units = "cm", dpi = 500, plot = pca_plot)

# Color by study
pca_df <- pca_df %>%
  mutate(study =  case_when(
    study == "SRP322015" ~ "GSE175839",
    study == "SRP168421" ~ "GSE122459",
    study == "SRP311059" ~ "GSE169080",
    study == "SRP296987" ~ "GSE162828",
    study == "SRP111941" ~ "GSE101437",
    study == "SRP136102" ~ "GSE112087",
    study == "SRP073191" ~ "GSE80183",
    TRUE ~"GSE72509"))




pca_plot2 <- ggplot(pca_df, aes(PC1, PC2, color=study, shape = DISEASE)) +
  geom_point(size=3) +
  labs(x = paste0("PC1: ",pc1_pct,"% variance"),
       y = paste0("PC2: ",pc2_pct,"% variance")) +
  scale_color_manual(values = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 
                                'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 
                                'GSE101437' = '#a1645c', 'GSE112087'='#599163', 
                                'GSE80183'='#755c91','GSE72509'='#e68a8a')) +
  # stat_ellipse(aes(group = DISEASE, color = DISEASE)) + 
  theme_classic() +
  labs(color = 'Study', shape = "Samples") +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),  
        axis.text.y = element_text(size = 25),
        legend.key.size = unit(1.5, "lines"))

save(pca_plot2, file = glue::glue("{FIGDIR}/pca_study_object.RData"))
png(filename = glue::glue("{FIGDIR}/pca_study.png"), height = 20, width = 22,
    units = "cm", res = 500)
print(pca_plot2)
dev.off()
#######
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