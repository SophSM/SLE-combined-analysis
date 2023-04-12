# Heatmaps

#####

library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
load("/mnt/Citosina/amedina/ssalazar/meta/out/LRT-dds.RData")
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
DGE.top <- DGE_names[ (DGE_names$baseMean.x > 50) & (DGE_names$padj < 0.05) & (abs(DGE_names$log2FoldChange) > 0.5), ] # 721

num_keep <- 20 # keep top n genes
rows_keep <- c(seq(1:num_keep), seq(((nrow(DGE.top)+1)-(num_keep)), nrow(DGE.top)))

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

# color map for expression

col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c('#34a624','white','#a709eb'))


colnames(ordered_norm) <- NULL
mat <- as.matrix(ordered_norm)

ha <- HeatmapAnnotation(Samples = all_data$DISEASE,
                        col = list(Samples = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b')))

row_ha <- rowAnnotation(logFC = l2_val, col = list(logFC =col_logFC))


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

h1 <-Heatmap(mat[rows_keep,], cluster_rows = F, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Norm counts', right_annotation = row_ha, bottom_annotation = study_ha, top_annotation = ha, column_km = 2, border = T)

png(paste0(outdir,'heatmap40_study.png'), width = 2000, height  = 2000, res = 300)
h1
dev.off()

# all genes

rownames(mat) <- NULL
l2_val <- as.matrix(DGE.top$log2FoldChange)
row_ha <- rowAnnotation(logFC = l2_val, col = list(logFC =col_logFC))
h2 <-Heatmap(mat, cluster_rows = F, name = 'Norm counts', bottom_annotation=study_ha, right_annotation = row_ha, top_annotation = ha, column_km = 2, border = T)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/"

png(paste0(outdir,'heatmapAll_studies.png'), width = 5000, height  = 10000, res = 300)
h2
dev.off()

#####
# heatmap with samples ordered (no clustering)

ordered_samples <- all_data[order(all_data$DISEASE),] # reorder all_data
ordered_samples$number <- rownames(ordered_samples)

# order count matrix
colnames(mat) <- seq(1:318)
ordered_norm <- mat[,order(match(colnames(mat), ordered_samples$number)) ]
colnames(ordered_norm) <- NULL

ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b')))
h1 <-Heatmap(mat[rows_keep,], cluster_rows = F, cluster_columns = F, row_labels = DGE.top[rows_keep,]$gene_name, name = 'Norm counts', right_annotation = row_ha, bottom_annotation = study_ha, top_annotation = ha)

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