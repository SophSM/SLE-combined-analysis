# Heatmaps

#####
library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
DIR = "/Users/sofiasalazar/Desktop/SLE-combined-analysis"
FIGDIR = glue::glue("{DIR}/figures_updated")
DGElist = data.table::fread(glue::glue("{DIR}/out/DGElist_all.csv"))
metadata <- read.table(glue::glue("{DIR}/data/all_data.csv"), header = T, row.names = 1,
                       sep = ",")

# load("/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData")
load(glue::glue("{DIR}/data/vsd2.RData"))
#####

metadata <- metadata %>%
  mutate(study =  case_when(
    study == "SRP322015" ~ "GSE175839",
    study == "SRP168421" ~ "GSE122459",
    study == "SRP311059" ~ "GSE169080",
    study == "SRP296987" ~ "GSE162828",
    study == "SRP111941" ~ "GSE101437",
    study == "SRP136102" ~ "GSE112087",
    study == "SRP073191" ~ "GSE80183",
    TRUE ~"GSE72509"))



# DGE_names <- merge(df_names, DGE, by = c('log2FoldChange', 'pvalue', 'padj')) # getting gene names
# dim(DGE_names) # 18982
# DGE_names<-DGE_names[,-c(10,11,12)]

DGE_protein_coding <- DGElist %>%
  filter(transcript_biotype == "protein_coding")

norm_counts <- as.data.frame(assay(vsd2))
dim(norm_counts) # 49465
zscore <- t(scale(t(norm_counts)))
zscore <- as.data.frame(zscore) %>%
  rownames_to_column("geneID") %>% 
  mutate(geneID = str_remove(geneID, "\\..*$")) %>%
  column_to_rownames("geneID") %>% 
  as.matrix()



###

# order dataframes
DGE_protein_coding <- DGE_protein_coding %>%
  arrange(desc(log2FoldChange))
ordered_norm <- zscore[DGE_protein_coding$ID, ]

dim(ordered_norm)  # 19191 318

# get top genes rows
DGE.top <- DGElist %>%
  filter(transcript_biotype == "protein_coding") %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  mutate(direction = ifelse(log2FoldChange > 1, "Up-regulated", "Down-regulated")) %>%
  group_by(direction) %>%
  arrange(desc(log2FoldChange)) %>%
  ungroup()
  
dim(DGE.top) # 290 8



# get only rows with common names
common_names = intersect(DGE.top$ID, rownames(zscore) )


# order norm counts according to logfoldchange in DGE list

ordered_norm <- ordered_norm[common_names,]
all(rownames(ordered_norm)==DGE.top$ID)
dim(ordered_norm) # 290 318

##############

# FOR ALL GENES

l2_val <- as.matrix(DGE.top$log2FoldChange)

# color map for log fold change
col_logFC <- colorRamp2(c(-2, -1, 0, 1,2), c("#2C7BB6","#ABD9E9", "white", "#e08788","#D7191C"))
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))
draw(row_ha)
# samples ordered (no clustering)
ordered_samples <- metadata[order(metadata$DISEASE),] # reorder all_data

# order count matrix
ordered_norm <- ordered_norm[, ordered_samples$samples]
all(colnames(ordered_norm)==ordered_samples$samples)


split = data.frame(Samples = ordered_samples$DISEASE) # make block split
ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 'GSE101437' = '#a1645c', 'GSE112087'='#599163', 'GSE80183'='#755c91','GSE72509'='#e68a8a')))
l2_val <- as.matrix(DGE.top$log2FoldChange)
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))
Colsplit <- data.frame(Samples = ordered_samples$DISEASE) # make block split
Rowsplit <- data.frame(Direction = DGE.top$direction)
Rowsplit$Direction <- factor(Rowsplit$Direction, levels = c("Up-regulated", "Down-regulated"))
col_exp <- colorRamp2(c(-3, -2,-0.5, 0, 0.5, 2, 3), c("#02487d", rev(RColorBrewer::brewer.pal(5, "RdYlBu")),"#870709"))


heat_ordered <-Heatmap(ordered_norm, cluster_rows = F, cluster_columns = F, name = 'Z-score',
                       left_annotation = row_ha, show_row_names = F, show_column_names = F,
                       column_split = Colsplit, col = col_exp,
                       row_split = Rowsplit)
                       

ht_list =study_ha%v%  ha %v% heat_ordered
png(filename = glue::glue("{FIGDIR}/heatmap_protein_coding.png"), height = 20, width = 30, units = "cm", res = 500)
  draw(ht_list)
dev.off()

pdf(glue::glue("{FIGDIR}/heatmap_protein_coding.pdf"), height = (20 / 2.54), width = (30 / 2.54))
draw(ht_list)
dev.off()

save(ht_list, file = glue::glue("{FIGDIR}/heatmap_protein_coding.RData"))

svg(filename = glue::glue("{FIGDIR}/heatmap_protein_coding.svg"),height = 8, width = 12)
draw(ht_list)
dev.off()
###########

# CLUSTERIZED FOR TOP GENES ONLY

# rows_keep <- rownames(DGE.top[(DGE.top$log2FoldChange > 3)|(DGE.top$log2FoldChange < (-3)),])

DEG_minimal <- DGE.top %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(20)
  
rows_keep <- DEG_minimal$ID
l2_val <-as.matrix(DEG_minimal$log2FoldChange)
colnames(l2_val)<- "logFC"
col_logFC <- colorRamp2(c(0, 1,2, 2.5, 3, 4),
                        c("#FFFFBF","#FEE090", "#FDAE61", "#F46D43","#D73027","#A50026"))



mat_top <- ordered_norm[rows_keep, ordered_samples$samples]
rownames(mat_top) <- DEG_minimal$gene_name
all(colnames(mat_top) == ordered_samples$samples)

col_exp <- colorRamp2(c(min(mat_top), -1, 0, 1, 2, 5),
                      c("#2C7BB6", "#ABD9E9", "#FFFFBF","#FDAE61", "#D7191C","darkred"))

ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 'GSE101437' = '#a1645c', 'GSE112087'='#599163', 'GSE80183'='#755c91','GSE72509'='#e68a8a')))

row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

h1 <-Heatmap(mat_top, cluster_columns = T, cluster_rows = F, name = 'Z-score', left_annotation = row_ha, 
             col = col_exp, column_km = 2, show_row_names = T, show_column_names = F)

h1_list <- study_ha %v% ha %v% h1

png(filename = glue::glue("{FIGDIR}/heatmapTOP.png"), height = 15, width = 15, units = "cm", res = 600)
draw(h1_list)
dev.off()

pdf(glue::glue("{FIGDIR}/heatmapTOP.pdf"), height = (15 / 2.54), width = (17 / 2.54))
draw(h1_list)
dev.off()

save(h1_list, file = glue::glue("{FIGDIR}/heatmapTOP.RData"))


# row clustering


h2 <-Heatmap(mat_top, cluster_rows = T, name = 'Z-score', left_annotation = row_ha, 
             col = col_exp, column_km = 2, show_column_names = F)

h2_list <- study_ha %v% ha %v% h2
png(filename = glue::glue("{FIGDIR}/heatmapTOP_clusterRows.png"), height = 15, width = 15, units = "cm", res = 600)
draw(h2_list)
dev.off()

save(h2_list, file = glue::glue("{FIGDIR}/heatmapTOP_clusterRows.RData"))

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
#   [1] ggplotify_0.1.3             patchwork_1.3.2            
# [3] RColorBrewer_1.1-3          biomaRt_2.64.0             
# [5] DESeq2_1.48.2               SummarizedExperiment_1.38.1
# [7] Biobase_2.68.0              MatrixGenerics_1.20.0      
# [9] matrixStats_1.5.0           GenomicRanges_1.60.0       
# [11] GenomeInfoDb_1.44.3         IRanges_2.42.0             
# [13] S4Vectors_0.46.0            BiocGenerics_0.54.1        
# [15] generics_0.1.4              shiny_1.13.0               
# [17] rrvgo_1.20.0                ggvenn_0.1.19              
# [19] circlize_0.4.17             gprofiler2_0.2.4           
# [21] ComplexHeatmap_2.25.2       lubridate_1.9.5            
# [23] forcats_1.0.1               stringr_1.6.0              
# [25] dplyr_1.2.0                 purrr_1.2.1                
# [27] readr_2.2.0                 tidyr_1.3.2                
# [29] tibble_3.3.1                tidyverse_2.0.0            
# [31] ggplot2_4.0.2              
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.18.0       jsonlite_2.0.0          shape_1.4.6.1          
# [4] umap_0.2.10.0           magrittr_2.0.4          farver_2.1.2           
# [7] rmarkdown_2.30          GlobalOptions_0.1.3     fs_1.6.6               
# [10] vctrs_0.7.1             memoise_2.0.1           askpass_1.2.1          
# [13] progress_1.2.3          S4Arrays_1.8.1          htmltools_0.5.9        
# [16] curl_7.0.0              gridGraphics_0.5-1      SparseArray_1.8.1      
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