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

workdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/"
fig_dir = glue::glue("{workdir}/figures/")
load("/mnt/Citosina/amedina/backup/lupus/sofi/vsd2.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")

all_data <-  read.csv(glue::glue("{workdir}/all_data.csv"), header = T)

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

ordered_samples <- all_data[order(all_data$DISEASE),] # reorder all_data

associated_genes <- read.csv(glue::glue("{workdir}/FIG2A-genes-pvals.csv"), header = T,
                             row.names = "X")

associated_genes <- associated_genes %>%
  arrange(p.value) %>% head(20)

norm_counts <- as.data.frame(assay(vsd2))
norm_counts <- tibble::rownames_to_column(norm_counts, "ID")
norm_counts_id <- name_short(norm_counts)
norm_counts_name <- inner_join(norm_counts_id, df_names, by = 'ID')

all(df_names$gene_name == norm_counts_name$gene_name)

norm_counts_name<- norm_counts_name[,-c(1,320:327)]
rownames(norm_counts_name) <- df_names$gene_name


counts <- as.matrix(norm_counts_name)
zscore <- t(scale(t(counts)))

dge_df <- df_names %>% filter(gene_name %in% associated_genes$Gene) %>% 
  arrange(padj)

# Heatmap
col_lfc <- colorRamp2(c(min(dge_df$log2FoldChange), 0, 0.5,  max(dge_df$log2FoldChange)),
                      c("blue", "white", "#ffb5b0","red"))
lfc_anno = rowAnnotation("log2FC" = dge_df$log2FoldChange, col = list("log2FC" = col_lfc))

mat <- zscore[dge_df$gene_name,ordered_samples$samples]
ha <- HeatmapAnnotation(Samples = ordered_samples$DISEASE,
                        col = list(Samples = c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')))
study_ha <-HeatmapAnnotation(Study = ordered_samples$study,
                             col = list(Study = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 'GSE101437' = '#a1645c', 'GSE112087'='#599163', 'GSE80183'='#755c91','GSE72509'='#e68a8a')))

split = data.frame(Samples = ordered_samples$DISEASE) # make block split
col_exp <- colorRamp2(c(min(mat),-5, 0, 3, max(mat)), c('blue', "lightblue1", 'white', 'red','darkred'))


heat_ordered <-Heatmap(mat, cluster_rows = T, cluster_columns = T, name = 'Z-score',
                       left_annotation = lfc_anno, show_row_names = T, show_column_names = F,
                       column_split = split,show_column_dend = F, top_annotation = c(study_ha,ha))


png(glue::glue("{fig_dir}/prevAsso_heatmap.png"), height = 15, width = 15, units = "cm", res = 300)
  draw(heat_ordered)
dev.off()

heat_clust <-Heatmap(mat, cluster_rows = T, cluster_columns = T, name = 'Z-score',
                       left_annotation = lfc_anno, show_row_names = T, show_column_names = F,
                       show_column_dend = T, column_km = 2, top_annotation = c(study_ha,ha))


png(glue::glue("{fig_dir}/prevAsso_heatmapClust.png"), height = 15, width = 15, units = "cm", res = 300)
draw(heat_clust)
dev.off()

