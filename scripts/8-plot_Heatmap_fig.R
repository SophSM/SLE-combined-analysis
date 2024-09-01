# Heatmaps

#####
library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
load("/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"
load("/mnt/Citosina/amedina/backup/lupus/sofi/vsd2.RData")
#####

all_data <- all_data %>%
  mutate(study =  case_when(
    study == "SRP322015" ~ "GSE175839",
    study == "SRP168421" ~ "GSE122459",
    study == "SRP311059" ~ "GSE169080",
    study == "SRP296987" ~ "GSE162828",
    study == "SRP111941" ~ "GSE101437",
    study == "SRP136102" ~ "GSE112087",
    study == "SRP073191" ~ "GSE80183",
    TRUE ~"GSE72509"))

DGE <- as.data.frame(results(dds2))
dim(DGE) # 49465

DGE_names <- merge(df_names, DGE, by = c('log2FoldChange', 'pvalue', 'padj')) # getting gene names
dim(DGE_names) # 18982
DGE_names<-DGE_names[,-c(10,11,12)]


norm_counts <- as.data.frame(assay(vsd2))
dim(norm_counts) # 49465
zscore <- t(scale(t(norm_counts)))
# zscore <- scale(norm_counts)
###

# order dataframes

DGE_names <- DGE_names[order(DGE_names$log2FoldChange, decreasing = TRUE),]
ordered_norm <- zscore[ order(match(rownames(zscore), rownames(DGE_names))), ]
dim(ordered_norm)  # 49465

# get top genes rows
DGE.top <- na.omit(DGE_names[(DGE_names$padj < 0.05) & (abs(DGE_names$log2FoldChange) > 1), ])
dim(DGE.top)# 331

table(DGE.top$log2FoldChange<(-1)) # 41
table(DGE.top$log2FoldChange>(1)) # 290

# match names in DGE dataframe and count matrix
rownames(DGE.top) <- DGE.top$ID
rownames(ordered_norm) <- gsub("\\..*","", rownames(ordered_norm))
df.list <- list(DGE.top, ordered_norm)

# get only rows with common names
common_names = Reduce(intersect, lapply(df.list, row.names))
df.list = lapply(df.list, function(x) { x[row.names(x) %in% common_names,] })

DGE.top <- df.list[[1]]
dim(DGE.top) # 331
ordered_norm <- df.list[[2]]
dim(ordered_norm) # 331

# order norm counts according to logfoldchange in DGE list

ordered_norm <- ordered_norm[rownames(DGE.top),]
all(rownames(ordered_norm)==rownames(DGE.top))
dim(ordered_norm) # 331 318

##############

# FOR ALL GENES


# getting log2FoldChange values
l2_val <- as.matrix(DGE.top$log2FoldChange)

# color map for log fold change
col_logFC <- colorRamp2(c(min(l2_val), 0, 2), c('dodgerblue3','white','firebrick3'))

# remove column names (sample IDs)
# colnames(ordered_norm) <- NULL
mat <- as.matrix(ordered_norm)

row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

# samples ordered (no clustering)
ordered_samples <- all_data[order(all_data$DISEASE),] # reorder all_data

# order count matrix
ordered_norm <- mat[, ordered_samples$samples]
all(colnames(ordered_norm)==ordered_samples$samples)


col_exp <- colorRamp2(c(min(ordered_norm), 0, 2, max(ordered_norm)), c('blue', 'white', 'red','darkred'))
split = data.frame(Samples = ordered_samples$DISEASE) # make block split
ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 'GSE101437' = '#a1645c', 'GSE112087'='#599163', 'GSE80183'='#755c91','GSE72509'='#e68a8a')))
l2_val <- as.matrix(DGE.top$log2FoldChange)
col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c('dodgerblue3','white','firebrick3'))
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

heat_ordered <-Heatmap(ordered_norm, cluster_rows = F, cluster_columns = F, name = 'Z-score',
                       left_annotation = row_ha, show_row_names = F, show_column_names = F,
                       column_split = split, col = col_exp)

ht_list =study_ha%v%  ha %v% heat_ordered
png(filename = paste0(outdir, 'heatmapFULL.png'), height = 20, width = 30, units = "cm", res = 300)
  draw(ht_list)
dev.off()

save(ht_list, file = paste0(outdir,"heatmapFULL.RData"))


###########

# CLUSTERIZED FOR TOP GENES ONLY

# rows_keep <- rownames(DGE.top[(DGE.top$log2FoldChange > 3)|(DGE.top$log2FoldChange < (-3)),])

rows_keep <- head(rownames(DGE.top), 20)
l2_val <-as.matrix(DGE.top[rows_keep,]$log2FoldChange)
colnames(l2_val)<- "logFC"
col_logFC <- colorRamp2(c(0,  min(l2_val),  max(l2_val)), c('white','#ffb3b3','#bd0202'))

mat <- mat[, ordered_samples$samples]
all(colnames(mat) == ordered_samples$samples)
ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 'GSE101437' = '#a1645c', 'GSE112087'='#599163', 'GSE80183'='#755c91','GSE72509'='#e68a8a')))

row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

h1 <-Heatmap(mat[rows_keep,], cluster_columns = T, cluster_rows = F, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Z-score', left_annotation = row_ha, 
             col = col_exp, column_km = 2, show_row_names = T, show_column_names = F)

h1_list <- study_ha %v% ha %v% h1

png(filename = paste0(outdir, 'heatmapTOP.png'), height = 15, width = 15, units = "cm", res = 300)
draw(h1_list)
dev.off()
save(h1_list, file = paste0(outdir,"heatmapTOP.RData"))


# row clustering

h2 <-Heatmap(mat[rows_keep,], cluster_rows = T, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Normalized counts', left_annotation = row_ha, 
             col = col_exp, column_km = 2, show_column_names = F)

h2_list <- study_ha %v% ha %v% h2
png(filename = paste0(outdir, 'heatmapTOP_clusterRows.png'), height = 15, width = 15, units = "cm", res = 300)
draw(h2_list)
dev.off()

save(h2_list, file = paste0(outdir,"heatmapTOP_clusterRows.RData"))

#######
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
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils    
# [8] datasets  methods   base     

# other attached packages:
#   [1] circlize_0.4.14             RColorBrewer_1.1-2         
# [3] ComplexHeatmap_2.9.3        forcats_0.5.1              
# [5] stringr_1.4.0               dplyr_1.0.10               
# [7] purrr_0.3.4                 readr_1.4.0                
# [9] tidyr_1.2.1                 tibble_3.1.3               
# [11] ggplot2_3.3.5               tidyverse_1.3.0            
# [13] DESeq2_1.30.1               SummarizedExperiment_1.20.0
# [15] Biobase_2.50.0              MatrixGenerics_1.2.1       
# [17] matrixStats_1.0.0           GenomicRanges_1.42.0       
# [19] GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [21] S4Vectors_0.28.1            BiocGenerics_0.36.1        

# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7           fs_1.5.0               lubridate_1.7.9.2     
# [4] bit64_4.0.5            doParallel_1.0.17      httr_1.4.2            
# [7] tools_4.0.2            backports_1.2.1        utf8_1.2.2            
# [10] R6_2.5.0               DBI_1.1.1              colorspace_2.0-2      
# [13] GetoptLong_1.0.5       withr_2.4.2            tidyselect_1.2.0      
# [16] bit_4.0.4              compiler_4.0.2         cli_3.6.0             
# [19] rvest_0.3.6            xml2_1.3.2             DelayedArray_0.16.3   
# [22] scales_1.1.1           genefilter_1.72.1      digest_0.6.27         
# [25] XVector_0.30.0         pkgconfig_2.0.3        dbplyr_2.2.1          
# [28] fastmap_1.1.0          GlobalOptions_0.1.2    rlang_1.0.6           
# [31] readxl_1.3.1           rstudioapi_0.13        RSQLite_2.2.7         
# [34] shape_1.4.6            generics_0.1.0         jsonlite_1.7.2        
# [37] BiocParallel_1.24.1    RCurl_1.98-1.3         magrittr_2.0.1        
# [40] GenomeInfoDbData_1.2.4 Matrix_1.3-4           Rcpp_1.0.7            
# [43] munsell_0.5.0          fansi_0.5.0            lifecycle_1.0.3       
# [46] stringi_1.6.2          zlibbioc_1.36.0        blob_1.2.1            
# [49] crayon_1.4.1           lattice_0.20-41        haven_2.3.1           
# [52] splines_4.0.2          annotate_1.68.0        hms_1.1.0             
# [55] locfit_1.5-9.4         pillar_1.6.2           rjson_0.2.20          
# [58] geneplotter_1.68.0     codetools_0.2-18       reprex_1.0.0          
# [61] XML_3.99-0.6           glue_1.4.2             modelr_0.1.8          
# [64] png_0.1-7              vctrs_0.5.1            foreach_1.5.2         
# [67] cellranger_1.1.0       gtable_0.3.0           clue_0.3-59           
# [70] assertthat_0.2.1       cachem_1.0.5           xtable_1.8-4          
# [73] broom_0.7.9            survival_3.5-5         iterators_1.0.13      
# [76] AnnotationDbi_1.52.0   memoise_2.0.0          cluster_2.1.0         
# [79] ellipsis_0.3.2    