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

#####

vsd2 <- vst(dds2, blind=FALSE)
mat <- assay(vsd2)
mm <- model.matrix(~DISEASE, colData(vsd2))
mat <- limma::removeBatchEffect(mat, batch=batch, design=mm)
assay(vsd2) <- mat

# save(vsd2, file=paste0(outdir, "vsd2.RData"))

DGE <- as.data.frame(results(dds2))
dim(DGE)

DGE_names <- merge(df_names, DGE, by = c('log2FoldChange', 'pvalue', 'padj')) # getting gene names
dim(DGE_names)
DGE_names<-DGE_names[,-c(10,11,12)]


norm_counts <- as.data.frame(assay(vsd2))
dim(norm_counts)

# order dataframes

DGE_names <- DGE_names[order(DGE_names$log2FoldChange, decreasing = TRUE),]
ordered_norm <- norm_counts[ order(match(rownames(norm_counts), rownames(DGE_names))), ]

# get top genes rows
DGE.top <- na.omit(DGE_names[(DGE_names$padj < 0.05) & (abs(DGE_names$log2FoldChange) > 1), ])
dim(DGE.top)# 331

table(DGE.top$log2FoldChange<(-1))
table(DGE.top$log2FoldChange>(1))
      # DGE.top <- DGE_names[ (DGE_names$baseMean.x > 50) & (DGE_names$padj < 0.05) & (abs(DGE_names$log2FoldChange) > 0.5), ] # 721

rownames(DGE.top) <- DGE.top$ID
rownames(ordered_norm) <- gsub("\\..*","", rownames(ordered_norm))
df.list <- list(DGE.top, ordered_norm)
common_names = Reduce(intersect, lapply(df.list, row.names))
df.list = lapply(df.list, function(x) { x[row.names(x) %in% common_names,] })

DGE.top <- df.list[[1]]
ordered_norm <- df.list[[2]]

ordered_norm <- ordered_norm[ order(match(rownames(ordered_norm), rownames(DGE.top))), ]

# getting log2 value for each gene we are keeping
l2_val <- as.matrix(DGE.top[rows_keep,]$log2FoldChange)
colnames(l2_val)<- "logFC"
min(as.vector(l2_val))

# color map for log fold change

col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c('#34a624','white','#a709eb'))


colnames(ordered_norm) <- NULL
mat <- as.matrix(ordered_norm)

# color map for expression
col_exp <- colorRamp2(c(min(mat), 5, mean(mat), 10, max(mat)), c('white', 'blue', 'yellow', 'red', 'darkred'))


#col_exp <- colorRamp2(c(min(mat), 5, 10, 15, 20, max(mat)), c('white', 'yellow', 'blue','green', 'orange', 'red'))

ha <- HeatmapAnnotation(Samples = all_data$DISEASE,
                        col = list(Samples = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b')))

row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))


h1 <-Heatmap(mat[rows_keep,], cluster_rows = F, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Norm counts', right_annotation = row_ha, top_annotation = ha, column_km = 2, border = T)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/"

# with 8 genes
png(paste0(outdir,'heatmap.png'), width = 2000, height  = 1000, res = 300)
h1
dev.off()

# with 20  genes

png(paste0(outdir,'heatmap20.png'), width = 2000, height  = 2000, res = 300)
h1
dev.off()

# with 40 genes

png(paste0(outdir,'heatmap40.png'), width = 2000, height  = 2000, res = 300)
h1
dev.off()


# with study annotation

study_ha <-HeatmapAnnotation(Study = all_data$study,
                             col = list(Study = c('SRP062966' = '#f5a142', 'SRP073191' = '#f5ef42', 'SRP111941' = '#2ef0e9', 'SRP136102' = '#f02eb3', 'SRP168421' = '#a1645c', 'SRP296987'='#599163', 'SRP311059'='#755c91','SRP322015'='#e68a8a')))

h1 <-Heatmap(mat[rows_keep,], cluster_rows = F, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Normalized counts', right_annotation = row_ha, bottom_annotation = study_ha, top_annotation = ha, column_km = 2, border = T)

png(paste0(outdir,'heatmap40_study.png'), width = 2000, height  = 2000, res = 300)
h1
dev.off()

# all genes

rownames(mat) <- NULL
l2_val <- as.matrix(DGE.top$log2FoldChange)
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))
h2 <-Heatmap(mat, cluster_rows = F, name = 'Normalized counts', bottom_annotation=study_ha, left_annotation = row_ha, top_annotation = ha, column_km = 2, border = F)

save(h2, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/heatmap_all_object.RData") 

png(paste0(outdir,'heatmapAll_studies.png'), width = 5000, height  = 10000, res = 300)
h2
dev.off()

#####
# heatmap with samples ordered (no clustering)

ordered_samples <- all_data[order(all_data$DISEASE),] # reorder all_data
ordered_samples$number <- rownames(ordered_samples)

# order count matrix
mat <- as.matrix(ordered_norm)
colnames(mat) <- seq(1:318)
ordered_norm <- mat[,order(match(colnames(mat), ordered_samples$number)) ]
colnames(ordered_norm) <- NULL
col_exp <- colorRamp2(c(min(mat), 5, mean(mat), 10, max(mat)), c('white', 'blue', 'yellow', 'red', 'darkred'))

split = data.frame(Samples = ordered_samples$DISEASE) # make block split
ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('SRP062966' = '#f5a142', 'SRP073191' = '#f5ef42', 'SRP111941' = '#2ef0e9', 'SRP136102' = '#f02eb3', 'SRP168421' = '#a1645c', 'SRP296987'='#599163', 'SRP311059'='#755c91','SRP322015'='#e68a8a')))
l2_val <- as.matrix(DGE.top$log2FoldChange)
col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c('#34a624','white','#a709eb'))
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

heat_ordered <-Heatmap(ordered_norm, cluster_rows = F, cluster_columns = F, 
             row_labels =NULL, name = 'Normalized counts',
             left_annotation = row_ha, bottom_annotation = study_ha, top_annotation = ha,
             column_split = split, col = col_exp)

save(heat_ordered, file = paste0(outdir,'heat_ordered.RData'))

png(paste0(outdir,'heatmap_ordered40.png'), width = 2000, height  = 2000, res = 300)
h1
dev.off()

#####

# heatmap with row (gene) clustering

ha <- HeatmapAnnotation(Samples = all_data$DISEASE,
                        col = list(Samples = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b')))
rownames(mat) <- NULL
colnames(mat) <- NULL

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
h3 <-Heatmap(mat, name = 'Norm counts', top_annotation = ha, column_km = 2, border = T,
             col = col_fun)

png(paste0(outdir,'heatmapAll_Rows.png'), width = 5000, height  = 10000, res = 300)
h3
dev.off()

h3 <-Heatmap(mat[rows_keep,], name = 'Norm counts', top_annotation = ha, column_km = 2, border = T)
png(paste0(outdir,'heatmap_Rows.png'), width = 5000, height  = 10000, res = 300)
h3
dev.off()

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