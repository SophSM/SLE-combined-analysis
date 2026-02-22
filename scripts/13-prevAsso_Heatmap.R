# Heatmap of previously associated genes to SLE
# Sofia Salazar
# 21 ago 2024
# ------------------------


library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(DESeq2)
name_short<-function(table){
  table$ID <- gsub("\\..*","", table$ID)
  return(table)
}
# Data

DIR = "/Users/sofiasalazar/Desktop/SLE-combined-analysis"
FIGDIR = glue::glue("{DIR}/figures_updated")
DGElist = data.table::fread(glue::glue("{DIR}/out/DGElist_all.csv"))
metadata <- read.table(glue::glue("{DIR}/data/all_data.csv"), header = T, row.names = 1,
                       sep = ",")
load(glue::glue("{DIR}/data/vsd2.RData"))

DGElist <- DGElist %>%
  filter(transcript_biotype == "protein_coding") %>%
  filter(gene_name != "")


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

ordered_samples <- metadata[order(metadata$DISEASE),] # reorder all_data

associated_genes <- read.csv(glue::glue("{DIR}/data/FIG2A-genes-pvals.csv"), header = T,
                             row.names = "X")

associated_genes <- associated_genes %>%
  arrange(p.value) %>% head(20)
DGElist_interest <- DGElist %>%
  filter(gene_name %in% associated_genes$Gene)

norm_counts <- as.data.frame(assay(vsd2))
norm_counts <- tibble::rownames_to_column(norm_counts, "ID")
norm_counts_id <- name_short(norm_counts)

norm_counts_name <- inner_join(norm_counts_id, DGElist_interest, by = 'ID')

all(norm_counts_name$gene_name == DGElist_interest$gene_name)

norm_counts_name<- norm_counts_name %>%
  column_to_rownames("gene_name") %>%
  dplyr::select(-c(baseMean, log2FoldChange,lfcSE,pvalue,padj,ID,
                   transcript_biotype))



norm_counts_name <- as.matrix(norm_counts_name)
zscore <- t(scale(t(norm_counts_name)))


# Heatmap


mat <- zscore[,ordered_samples$samples]

ordered_genes <- DGElist_interest %>%
  column_to_rownames("gene_name") %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"))
ordered_genes <- ordered_genes[rownames(mat),]
all(rownames(ordered_genes) == rownames(mat))

ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 'GSE101437' = '#a1645c', 'GSE112087'='#599163', 'GSE80183'='#755c91','GSE72509'='#e68a8a')))

col_lfc <- colorRamp2(c(min(ordered_genes$log2FoldChange), 0, 0.5,  max(ordered_genes$log2FoldChange)),
                      c("blue", "#FFFFBF", "#ffb5b0","red"))
lfc_anno = rowAnnotation("log2FC" = ordered_genes$log2FoldChange, col = list("log2FC" = col_lfc))
# l2_val <- ordered_genes$log2FoldChange
split = data.frame(Samples = ordered_samples$DISEASE) # make block split

col_exp <- colorRamp2(c(-4, -3,-1, 0,  1, 3, 4), 
                      c("#02487d","#4575B4", "#91BFDB", "#FFFFBF", "#FC8D59", "#D73027", 'darkred'))

[1] "#4575B4" "#91BFDB" "#E0F3F8" "#FFFFBF" "#FEE090" "#FC8D59" "#D73027"
                      c('blue', "lightblue1", 'white', 'red','darkred'))

"#4575B4" "#74ADD1" "#ABD9E9" "#E0F3F8" "#FFFFBF" "#FEE090" "#FDAE61"
[8] "#F46D43" "#D73027"-

Rsplit <- data.frame(Direction = ordered_genes$direction)

Rsplit$Direction <- factor(Rsplit$Direction, levels = c("Up-regulated", "Down-regulated"))


heat_ordered <-Heatmap(mat, cluster_rows = T, cluster_columns = T, name = 'Z-score',
                       col = col_exp,
                       left_annotation = lfc_anno, show_row_names = T,
                       column_gap = unit(2, "mm"),
                       row_gap = unit(2, "mm"),
                       show_column_names = F,
                       column_split = split, row_split = Rsplit, show_column_dend = F)

heat_list <- study_ha %v% ha %v% heat_ordered
png(glue::glue("{FIGDIR}/prevAsso_heatmap.png"), height = 15, width = 20, units = "cm", res = 500)
  draw(heat_list)
dev.off()

heat_clust <-Heatmap(mat, cluster_rows = F, cluster_columns = T, name = 'Z-score',
                       left_annotation = lfc_anno, show_row_names = T, show_column_names = F,
                       show_column_dend = T, column_km = 2, 
                     col = col_exp,
                      column_gap = unit(2, "mm"),
                      row_gap = unit(2, "mm"),
                      row_split = Rsplit)

heat_list2 <- study_ha %v% ha %v% heat_clust

png(glue::glue("{FIGDIR}/prevAsso_heatmapClust.png"), height = 15, width = 17, units = "cm", res = 300)
draw(heat_list2)
dev.off()

