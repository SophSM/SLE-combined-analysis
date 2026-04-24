library(DESeq2)
library(tidyverse)


# genes interest
DGE_df <-data.table::fread("/Users/sofiasalazar/Desktop/SLE-combined-analysis/out/DGElist_protein_coding.csv")
genes_interest <- c("F10", "F11", "F12", "F13A1", "F13B", "F2", "F5", "F7",
                    "F8", "F9", "FGA", "FGB", "FGG", "GP1BA", "HRG", "KLKB1",
                    "KNG1", "LMAN1", "MCFD2", "MTHFR", "PLAT", "PROC", "PROS1",
                    "PROZ", "SERPINC1", "SERPIND1", "SERPINE1", "SERPINF2",
                    "THBD", "VKORC1", "VWF")

DEG_interest <- DGE_df %>%
  filter(gene_name %in% genes_interest) %>%
  dplyr::select(-transcript_biotype) %>%
  unique() %>%
  mutate(direction = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < 0 & padj < 0.05 ~ "Downregulated",
    T ~ "No change"),
         significant = case_when(padj < 0.05 ~ "Yes", 
                                 padj > 0.05 ~ "No",
                                 T ~ "No"),
    label = ifelse(direction != "No change", gene_name, ""))

# write.table(DEG_interest, file = "/Users/sofiasalazar/Desktop/SLE-combined-analysis/out/DEGs_ninive_list.txt",
            # sep = "\t", quote = F, row.names = F)

DEG_interest %>%
  filter(padj < 0.05) %>% 
  dplyr::select(c(gene_name, log2FoldChange)) %>%
  arrange(desc(abs(log2FoldChange)))
# --- Volcano plot ---

top_genes <- DEG_interest %>%
  filter(padj < 0.05)
png(filename  = "/Users/sofiasalazar/Desktop/SLE-combined-analysis/figures_updated/volcano_reproductiveGenes.png",
    height = 10, width = 10, res = 500, units = "cm")
ggplot(DEG_interest, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
  geom_point() +
  theme_minimal() +
  ggrepel::geom_label_repel(data = top_genes,
                            mapping = aes(x = log2FoldChange, 
                                          y = -log(padj,10), label = gene_name),
                            size = 4) +
  scale_color_manual(values = c("Upregulated" = "firebrick", 
                                "Downregulated"="dodgerblue",
                                "No change" = "gray15")) +
  geom_hline(yintercept = -log10(0.05), color = "gray30", lty = "dashed") +
  labs(x = "log2 FC", y = "-log10(Adj. p-value)", color = "Direction") +
  ylim(0, 12) +
  guides(color = "none")
dev.off()

library(ComplexHeatmap)
library(circlize)

rownames(DEG_interest ) = NULL
mat_deg <- DEG_interest %>%
  dplyr::select(c(gene_name, log2FoldChange)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

pval_annot <- rowAnnotation("Is significant" = DEG_interest$significant,
                            col = list("Is significant" = c("Yes" = "#007300", "No"="lightgray")),
                            gp = gpar(col = "black", lwd = 2), show_annotation_name = F)

col_fun = colorRamp2(breaks = c(-0.7, -0.5, -0.1, -0.05, 0, 0.05, 0.1, 0.5, 1),
                     colors = c(RColorBrewer::brewer.pal(9, "RdYlBu")))
Heatmap(mat_deg, name = "log2FC", col = col_fun,
        rect_gp = gpar(col = "black", lwd = 2), left_annotation = pval_annot,
        show_column_names = F, show_row_dend = F)
