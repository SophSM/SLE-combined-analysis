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
load("/mnt/Citosina/amedina/ssalazar/meta/combined/vsd2.RData")
#####

vsd2 <- vst(dds2, blind=FALSE)
mat <- assay(vsd2)
mm <- model.matrix(~DISEASE, colData(vsd2))
mat <- limma::removeBatchEffect(mat, batch=batch, design=mm)
assay(vsd2) <- mat

# save(vsd2, file=paste0(outdir, "vsd2.RData"))

DGE <- as.data.frame(results(dds2))
dim(DGE) # 49465

DGE_names <- merge(df_names, DGE, by = c('log2FoldChange', 'pvalue', 'padj')) # getting gene names
dim(DGE_names) # 18982
DGE_names<-DGE_names[,-c(10,11,12)]


norm_counts <- as.data.frame(assay(vsd2))
dim(norm_counts) # 49465

###

# order dataframes

DGE_names <- DGE_names[order(DGE_names$log2FoldChange, decreasing = TRUE),]
ordered_norm <- norm_counts[ order(match(rownames(norm_counts), rownames(DGE_names))), ]
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
ordered_norm <- ordered_norm[ order(match(rownames(ordered_norm), rownames(DGE.top))), ] 
dim(ordered_norm) # 331

##############

# FOR ALL GENES


# getting log2FoldChange values
l2_val <- as.matrix(DGE.top$log2FoldChange)

# color map for log fold change
col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c('#34a624','white','#a709eb'))

# remove column names (sample IDs)
colnames(ordered_norm) <- NULL
mat <- as.matrix(ordered_norm)

# remove row names (gene names)
rownames(mat) <- NULL
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

# samples ordered (no clustering)
ordered_samples <- all_data[order(all_data$DISEASE),] # reorder all_data
ordered_samples$number <- rownames(ordered_samples)

# order count matrix
colnames(mat) <- seq(1:318)
ordered_norm <- mat[,order(match(colnames(mat), ordered_samples$number)) ]
colnames(ordered_norm) <- NULL
rownames(ordered_norm) <- NULL

col_exp <- colorRamp2(c(min(ordered_norm), 5, mean(ordered_norm), 10, max(ordered_norm)), c('white', 'blue', 'yellow', 'red', 'darkred'))

split = data.frame(Samples = ordered_samples$DISEASE) # make block split
ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('SRP062966' = '#f5a142', 'SRP073191' = '#f5ef42', 'SRP111941' = '#2ef0e9', 'SRP136102' = '#f02eb3', 'SRP168421' = '#a1645c', 'SRP296987'='#599163', 'SRP311059'='#755c91','SRP322015'='#e68a8a')))
l2_val <- as.matrix(DGE.top$log2FoldChange)
col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c('#34a624','white','#a709eb'))
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

heat_ordered <-Heatmap(ordered_norm, cluster_rows = F, cluster_columns = F, name = 'Normalized counts',
                       left_annotation = row_ha, 
                       bottom_annotation = study_ha, top_annotation = ha,
                       column_split = split, col = col_exp)
save(heat_ordered, file = paste0(outdir,'heat_ordered.RData'))


###########

# CLUSTERIZED FOR TOP GENES ONLY
rownames(DGE.top) <- 1:dim(DGE.top)[1]
# rows_keep <- rownames(DGE.top[(DGE.top$log2FoldChange > 3)|(DGE.top$log2FoldChange < (-3)),])

rows_keep <- c(1:20)
l2_val <-as.matrix(DGE.top[rows_keep,] $log2FoldChange)
colnames(l2_val)<- "logFC"

# remove column names (sample IDs)
colnames(ordered_norm) <- NULL
mat <- as.matrix(ordered_norm)

# expression color
col_exp <- colorRamp2(c(min(mat[rows_keep,]), 5, mean(mat[rows_keep,]), 10, max(mat[rows_keep,])), c('white', 'blue', 'yellow', 'red', 'darkred'))

study_ha <-HeatmapAnnotation(Study = all_data$study,
                             col = list(Study = c('SRP062966' = '#f5a142', 'SRP073191' = '#f5ef42', 'SRP111941' = '#2ef0e9', 'SRP136102' = '#f02eb3', 'SRP168421' = '#a1645c', 'SRP296987'='#599163', 'SRP311059'='#755c91','SRP322015'='#e68a8a')))
col_logFC <- colorRamp2(c(min(l2_val) - 1 ,max(l2_val)), c('white','#a709eb'))

ha <- HeatmapAnnotation(Samples = all_data$DISEASE,
                        col = list(Samples = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b')))

row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))

h1 <-Heatmap(mat[rows_keep,], cluster_rows = F, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Normalized counts', left_annotation = row_ha, 
             top_annotation = ha, bottom_annotation = study_ha, col = col_exp,
             column_km = 2)

save(h1, file = paste0(outdir,"cluster_heatmap.RData"))


# row clustering

h2 <-Heatmap(mat[rows_keep,], cluster_rows = T, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Normalized counts', left_annotation = row_ha, 
             top_annotation = ha, bottom_annotation = study_ha, col = col_exp,
             column_km = 2)
save(h2, file = paste0(outdir,"cluster_heatmap_rows.RData"))
