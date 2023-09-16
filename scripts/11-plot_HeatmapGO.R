
# Heatmap for GO terms and genes
######
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
######

load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData")
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"

data <- df_names %>% 
  mutate(Expression = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))


keytypes(org.Hs.eg.db)
original_gene_list = data$log2FoldChange
names(original_gene_list) <- data$ID
ids <- bitr(names(original_gene_list), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')

# GET TOP GENES

top <- 70
top_genes <- bind_rows(data %>% 
                         filter(Expression == 'Up-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top),
                       data %>% 
                         filter(Expression == 'Down-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top)
)

DGE_genes <- data[data$Expression != 'Unchanged',]

### lista para el lab
list_lab <- DGE_genes[,c(8,10,3)]
write.csv(list_lab, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/topGenes_list.csv')



list_lab_up <- list_lab[list_lab$Expression == 'Up-regulated',]
list_lab_up <- list_lab_up[order(-list_lab_up$log2FoldChange), ]
rownames(list_lab_up) <- seq(1:nrow(list_lab_up))

list_lab_down <- list_lab[list_lab$Expression == 'Down-regulated',]
list_lab_down <- list_lab_down[order(list_lab_down$log2FoldChange), ]
rownames(list_lab_down) <- seq(1:nrow(list_lab_down))

write.csv(list_lab_up, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/topGenesUp_list.csv')
write.csv(list_lab_down, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/topGenesDown_list.csv')

######

sigGenes <- top_genes$gene_name
original_gene_list <- top_genes$log2FoldChange
names(original_gene_list) <- top_genes$ID
# names(original_gene_list) <- top_genes$gene_name
ids <- bitr(names(original_gene_list), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
top_genes <- top_genes[!duplicated(top_genes$ID),]
data2 = top_genes[top_genes$ID %in% dedup_ids$ENSEMBL,]
data2$Y = dedup_ids$ENTREZID
names(original_gene_list) <- data2$Y

# GO
go_gene_list<-na.omit(original_gene_list)
go_gene_list = sort(go_gene_list, decreasing = TRUE)
save(go_gene_list, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/GO_gene_list.RData')

gseGO_res <- gseGO(geneList = go_gene_list,
                   ont = 'ALL',
                   OrgDb = org.Hs.eg.db,
                   minGSSize = 3,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   keyType = 'ENTREZID')
save(gseGO_res, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/GO_results.RData")

########
load("/mnt/Citosina/amedina/ssalazar/meta/combined/GO_results.RData")
# GO with all DE genes

sigGenes <- DGE_genes$gene_name
original_gene_list <- DGE_genes$log2FoldChange
names(original_gene_list) <- DGE_genes$ID

# names(original_gene_list) <- top_genes$gene_name
ids <- bitr(names(original_gene_list), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
DGE_genes <- DGE_genes[!duplicated(DGE_genes$ID),]
data2 = DGE_genes[DGE_genes$ID %in% dedup_ids$ENSEMBL,]
data2$Y = dedup_ids$ENTREZID
names(original_gene_list) <- data2$Y

go_gene_list<-na.omit(original_gene_list)
go_gene_list = sort(go_gene_list, decreasing = TRUE)
save(go_gene_list, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/GO_allDEG_list.RData')

gseGO_res <- gseGO(geneList = go_gene_list,
                   ont = 'ALL',
                   OrgDb = org.Hs.eg.db,
                   minGSSize = 3,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   keyType = 'ENTREZID')
save(gseGO_res, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/GO_allDEG_results.RData")

#########
# load

load("/mnt/Citosina/amedina/ssalazar/meta/combined/GO_allDEG_results.RData")

########

symbolGO_results <- setReadable(gseGO_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
symbolGO_results <- symbolGO_results@result

# order symbolGO_results by p.value
symbolGO_results <- symbolGO_results[order(symbolGO_results$p.adjust),]
gseSubset <- subset(symbolGO_results, p.adjust <= 0.05 )


annGSEA <- data.frame(row.names = sigGenes)
for (j in 1:length(sigGenes)) {
  gene <- sigGenes[j]
  pattern <- gene
  for (k in 1:nrow(gseSubset)){ # parse through GO terms 
    if (any(grepl(pattern, gseSubset$core_enrichment[k]))) {
      annGSEA[j,k] <- 1
    } else {
      annGSEA[j,k] <- 0
    }
  }
}
colnames(annGSEA) <- gseSubset[,3]
# remove terms with no overlapping genes
annGSEA <- annGSEA[,apply(annGSEA, 2, mean)!=0]

# remove genes with no overlapping terms
annGSEA <- annGSEA[apply(annGSEA, 1, mean)!=0,]

annGSEA[1:5,1:5]

save(annGSEA, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/annGSEA_allDEG_GO.RData')
########
load('/mnt/Citosina/amedina/ssalazar/meta/combined/annGSEA_allDEG_GO.RData')

# match the order of rownames in top_genes with annGSEA
rownames(top_genes) <- top_genes$gene_name
topTableAligned <- top_genes[which(rownames(top_genes) %in% rownames(annGSEA)),]
topTableAligned <- topTableAligned[match(rownames(annGSEA), rownames(topTableAligned)),]
all(rownames(topTableAligned) == rownames(annGSEA)) # TRUE

### all DEG

rownames(DGE_genes) <- DGE_genes$gene_name
topTableAligned <- DGE_genes[which(rownames(DGE_genes) %in% rownames(annGSEA)),]
topTableAligned <- topTableAligned[match(rownames(annGSEA), rownames(topTableAligned)),]
all(rownames(topTableAligned) == rownames(annGSEA)) # TRUE

# colour bar for fold changes for sigGenes
dfFoldChangeGenes <- data.frame(
  topTableAligned[which(rownames(topTableAligned) %in% rownames(annGSEA)), 'log2FoldChange'])

# merge both
dfGeneAnno <- data.frame(dfFoldChangeGenes)
colnames(dfGeneAnno) <- c('Log2FC')

colours <- colorRamp2(c(min(dfGeneAnno$Log2FC),0, max(dfGeneAnno$Log2FC)), c('royalblue','white','#b50f04'))

haGenes <- rowAnnotation(
  df = dfGeneAnno,
  col = list(Log2FC =colours),
  width = unit(1,'cm'),
  annotation_name_side = 'top')

# bottom annotation with enriched terms

haTerms <- HeatmapAnnotation(
  text = anno_text(
    colnames(annGSEA),
    rot = 45,
    just = 'right',
    gp = gpar(fontsize = 12)),
  annotation_height = unit(8, 'cm'),
  annotation_name_side = 'left')

# bottom annotation with GO source for each term

term_source <- symbolGO_results %>% dplyr::select(ONTOLOGY)
 
sourceTermsha <- HeatmapAnnotation(
  Source = term_source$ONTOLOGY,
  col = list(Source = c('BP' = '#3C6997', 'CC' = '#d9c621', "MF" = '#b30039'))
)

# bottom annotation with p.value of each term
df.pvals <- data.frame(gseSubset$p.adjust)
colnames(df.pvals) <- c('p.value')

pvalHa <- HeatmapAnnotation(
  p.value = df.pvals$p.value,
  col = list(p.value = colorRamp2(c(min(df.pvals$p.value), median(df.pvals$p.value), max(df.pvals$p.value)), c('purple4','purple', '#fcd7f6'))))
# HEATMAP

annGSEA <- as.matrix(annGSEA)
hmapGSEA <- Heatmap(annGSEA,
                    name = 'GO enrichment',
                    # split = dfGeneAnno[,2],
                    
                    col = c('0' = 'white', '1' = 'forestgreen'),
                    
                    rect_gp = gpar(col = 'grey85'),
                    
                    cluster_rows = TRUE,
                    show_row_dend = TRUE,
                    row_title = 'Top Genes',
                    row_title_side = 'left',
                    row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_title_rot = 90,
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_names_side = 'left',
                    row_dend_width = unit(35, 'mm'),
                    
                    cluster_columns = FALSE,
                    show_column_dend = TRUE,
                    column_title = 'Enriched terms',
                    column_title_side = 'top',
                    column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                    column_title_rot = 0,
                    show_column_names = FALSE,
                    
                    show_heatmap_legend = FALSE,
                    
                    # clustering_distance_columns = 'euclidean',
                    # clustering_method_columns = 'ward.D2',
                    clustering_distance_rows = 'euclidean',
                    clustering_method_rows = 'ward.D2',
                    
                    bottom_annotation = c(sourceTermsha, pvalHa, haTerms)
)

png(filename = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/GOenriched_allDEG_HMap_wTerm.png", width = 5000, height = 5000, res = 300)
draw(hmapGSEA + haGenes,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()

png(filename = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/GOenriched_allDEG_HMap.png", width = 5000, height = 5000, res = 300)
draw(hmapGSEA + haGenes,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
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
#   [1] circlize_0.4.14        ComplexHeatmap_2.9.3   org.Hs.eg.db_3.12.0   
# [4] AnnotationDbi_1.52.0   IRanges_2.24.1         S4Vectors_0.28.1      
# [7] Biobase_2.50.0         BiocGenerics_0.36.1    forcats_0.5.1         
# [10] stringr_1.4.0          dplyr_1.0.10           purrr_0.3.4           
# [13] readr_1.4.0            tidyr_1.2.1            tibble_3.1.3          
# [16] tidyverse_1.3.0        ggplot2_3.3.5          clusterProfiler_3.18.1
# [19] DOSE_3.16.0            enrichplot_1.10.2     

# loaded via a namespace (and not attached):
#   [1] fgsea_1.16.0        colorspace_2.0-2    rjson_0.2.20       
# [4] ellipsis_0.3.2      qvalue_2.22.0       GlobalOptions_0.1.2
# [7] fs_1.5.0            clue_0.3-59         rstudioapi_0.13    
# [10] farver_2.1.0        graphlayouts_0.7.1  ggrepel_0.9.1      
# [13] bit64_4.0.5         fansi_0.5.0         scatterpie_0.1.6   
# [16] lubridate_1.7.9.2   xml2_1.3.2          codetools_0.2-18   
# [19] splines_4.0.2       doParallel_1.0.17   cachem_1.0.5       
# [22] GOSemSim_2.16.1     polyclip_1.10-0     jsonlite_1.7.2     
# [25] broom_0.7.9         cluster_2.1.0       GO.db_3.12.1       
# [28] dbplyr_2.2.1        png_0.1-7           ggforce_0.3.3      
# [31] BiocManager_1.30.21 compiler_4.0.2      httr_1.4.2         
# [34] rvcheck_0.1.8       backports_1.2.1     assertthat_0.2.1   
# [37] Matrix_1.3-4        fastmap_1.1.0       cli_3.6.0          
# [40] tweenr_1.0.2        tools_4.0.2         igraph_1.2.8       
# [43] gtable_0.3.0        glue_1.4.2          reshape2_1.4.4     
# [46] DO.db_2.9           fastmatch_1.1-0     Rcpp_1.0.7         
# [49] cellranger_1.1.0    vctrs_0.5.1         iterators_1.0.13   
# [52] ggraph_2.0.5        rvest_0.3.6         lifecycle_1.0.3    
# [55] MASS_7.3-53         scales_1.1.1        tidygraph_1.2.0    
# [58] hms_1.1.0           RColorBrewer_1.1-2  memoise_2.0.0      
# [61] gridExtra_2.3       downloader_0.4      stringi_1.6.2      
# [64] RSQLite_2.2.7       foreach_1.5.2       BiocParallel_1.24.1
# [67] shape_1.4.6         rlang_1.0.6         pkgconfig_2.0.3    
# [70] matrixStats_1.0.0   lattice_0.20-41     cowplot_1.1.1      
# [73] shadowtext_0.0.8    bit_4.0.4           tidyselect_1.2.0   
# [76] plyr_1.8.6          magrittr_2.0.1      R6_2.5.0           
# [79] generics_0.1.0      DBI_1.1.1           pillar_1.6.2       
# [82] haven_2.3.1         withr_2.4.2         modelr_0.1.8       
# [85] crayon_1.4.1        utf8_1.2.2          viridis_0.5.1      
# [88] GetoptLong_1.0.5    readxl_1.3.1        data.table_1.14.0  
# [91] blob_1.2.1          reprex_1.0.0        digest_0.6.27      
# [94] munsell_0.5.0       viridisLite_0.4.0 