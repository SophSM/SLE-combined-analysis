# Heatmap of reproductive-relevant genes
# Sofia Salazar
# --------------------

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
DIR = "/Users/sofiasalazar/Desktop/SLE-combined-analysis"
FIGDIR = glue::glue("{DIR}/figures_updated")
DGElist = data.table::fread(glue::glue("{DIR}/out/DGElist_all.csv"))
metadata <- read.table(glue::glue("{DIR}/data/all_data.csv"), header = T, row.names = 1,
                       sep = ",")

load(glue::glue("{DIR}/data/vsd2.RData"))

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

DGE_protein_coding <- DGE_protein_coding %>%
  arrange(desc(log2FoldChange))
ordered_norm <- zscore[DGE_protein_coding$ID, ]

# ---------------------------subset-------------------
significant_reproductive <- c("F5", "F2","F12","F8","VKORC1","PROZ","SERPINF2","F10")
genes_reproductive <- c("F10", "F11", "F12", "F13A1", "F13B", "F2", "F5", "F7",
                        "F8", "F9", "FGA", "FGB", "FGG", "GP1BA", "HRG", "KLKB1",
                        "KNG1", "LMAN1", "MCFD2", "MTHFR", "PLAT", "PROC", "PROS1",
                        "PROZ", "SERPINC1", "SERPIND1", "SERPINE1", "SERPINF2",
                        "THBD", "VKORC1", "VWF")

DGE.reproductive<- DGElist %>%
  filter(transcript_biotype == "protein_coding") %>%
  filter(gene_name %in% genes_reproductive) %>%
  # filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated")) %>%
  group_by(direction) %>%
  arrange(desc(log2FoldChange)) %>%
  ungroup()

dim(DGE.reproductive) # 31 9

DGE.reproductive %>%
  filter(padj < 0.05) %>% nrow() # signficant

# get only rows with common names
common_names = intersect(DGE.reproductive$ID, rownames(zscore) )


# order norm counts according to logfoldchange in DGE list

ordered_norm <- ordered_norm[common_names,]
all(rownames(ordered_norm)==DGE.reproductive$ID)
dim(ordered_norm) # 31 318

# ---- heatmap----
l2_val <- as.matrix(DGE.reproductive$log2FoldChange)

# color map for log fold change
col_logFC <- colorRamp2(c(-0.4, -0.1,0, 0.1, 0.4, 0.7), c("#91BFDB","#E0F3F8", "#FFFFBF","#FEE090","#FC8D59", "#D73027"))


# "#D73027" "#FC8D59" "#FEE090" "#FFFFBF" "#E0F3F8" "#91BFDB" "#4575B4"
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
l2_val <- as.matrix(DGE.reproductive$log2FoldChange)
row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC =col_logFC))
Colsplit <- data.frame(Samples = ordered_samples$DISEASE) # make block split
Rowsplit <- data.frame(Direction = DGE.reproductive$direction)
Rowsplit$Direction <- factor(Rowsplit$Direction, levels = c("Up-regulated", "Down-regulated"))

col_exp <- colorRamp2(c(-4, -2,-0.5, 0, 0.5, 2, 4), c("#02487d", rev(RColorBrewer::brewer.pal(5, "RdYlBu")),"#870709"))


heat_ordered <-Heatmap(ordered_norm, cluster_rows = F, cluster_columns = F, name = 'Z-score',
                       left_annotation = row_ha, show_row_names = T, show_column_names = F,
                       column_split = Colsplit,
                       col = col_exp,
                       row_split = Rowsplit, row_labels = DGE.reproductive$gene_name,
                      row_title = " ")


ht_list =study_ha%v%  ha %v% heat_ordered
png(filename = glue::glue("{FIGDIR}/heatmap_reproductive_genes.png"), 
    height = 15, width = 15, units = "cm", res = 500)
draw(ht_list)
dev.off()
