# Violin Plots
#####
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(perm)
DIR = "/Users/sofiasalazar/Desktop/SLE-combined-analysis"
FIGDIR = glue::glue("{DIR}/figures_updated")
DGElist = data.table::fread(glue::glue("{DIR}/out/DGElist_all.csv"))
load(glue::glue("{DIR}/data/vsd2.RData"))
metadata <- read.table(glue::glue("{DIR}/data/all_data.csv"), header = T, row.names = 1,
                       sep = ",")
#----------------------------------------------------------#

norm_counts_name <- read.csv('/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/normcounts_name.csv',
                             header = T, row.names = "X")

DGElist_clean <-  DGElist %>%
  filter(gene_name != "") %>%
  filter(!is.na(gene_name)) %>%
  filter(gene_name != "NA") %>%
  filter(transcript_biotype == "protein_coding")
norm_counts <- as.data.frame(assay(vsd2))

norm_counts_name <- norm_counts %>%
  rownames_to_column("geneID") %>% 
  mutate(geneID = str_remove(geneID, "\\..*$")) %>%
  left_join(DGElist_clean %>%
              dplyr::select(c(ID, gene_name)) %>%
              unique(), by = c("geneID"="ID")) %>%
  dplyr::select(-geneID)
  
zscore <- t(scale(t(norm_counts)))
zscore <- as.data.frame(zscore) %>%
  rownames_to_column("geneID") %>% 
  mutate(geneID = str_remove(geneID, "\\..*$")) %>%
  left_join(DGElist_clean %>%
              dplyr::select(c(ID, gene_name)) %>%
              unique(), by = c("geneID"="ID")) %>%
  dplyr::select(-geneID)


##########
# OUR SELECTED INTERFERON SIGNATURE GENES

# compute pvalues
all_pvalues <- c()
list.genes <- c("CCL2", "IFIT1", "RSAD2", "IFI44L", 
                "IFI44", "USP18", "ISG15", "IFI27", "SIGLEC1")

for (gene in list.genes){
  counts.gene <- norm_counts_name %>% 
    filter(gene_name == gene) %>%
    dplyr::select(-gene_name)
  counts.gene<- as.data.frame(t(counts.gene))
  
  # test
  expression = counts.gene[metadata$samples,1]
  sample = metadata$DISEASE
  df.gene <- data.frame(expression, sample)
  
  control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
  sle.e<- df.gene[df.gene$sample=='SLE',]$expression
  test <- permTS(control.e, sle.e)
  pval <- (test$p.value)
  print(gene)
  all_pvalues <- c(all_pvalues,pval)
  
}


# Adjust p values

adj_pvalues <- p.adjust(all_pvalues, method = "fdr", n = length(all_pvalues))
length(adj_pvalues)

pvals_df <- data.frame(pval = adj_pvalues, gene_name = list.genes)

violin_df <- norm_counts_name %>%
  filter(gene_name %in% pvals_df$gene) %>%
  # column_to_rownames("gene_name") %>%
  pivot_longer(names_to = "sample", cols = -gene_name) %>%
  left_join(metadata, by = c("sample" = "samples")) %>%
  left_join(pvals_df, by = c("gene_name" = "gene_name"))

df_mean <- violin_df %>%
  group_by(DISEASE, gene_name) %>%
  summarize(average = mean(value)) %>%
  ungroup()
violin_df <- violin_df %>%
  left_join(df_mean, by = c("gene_name" = "gene_name", "DISEASE"="DISEASE")) %>%
  arrange(pval)

violin_df$gene_name <- factor(violin_df$gene_name, levels = unique(violin_df$gene_name))

p_inferferon <- ggplot(violin_df, aes(x = DISEASE, y = value, fill = DISEASE)) +
  geom_violin(trim=FALSE) +
  facet_wrap(~as.factor(gene_name), ncol = 3, scales = "free_y") +
  geom_boxplot(width=0.2, color = 'black', fill=NA) +
  geom_point(df_mean, mapping = aes(x = DISEASE, y = average),
             color = 'white', shape = 18, size = 3) +
  scale_fill_manual(values=c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')) +
  geom_point(aes(x = DISEASE, y = average), color = 'white', shape = 18, size = 3) +
  geom_line(aes(x = DISEASE, y = average, group = 1), linetype = "dashed") +
  theme_minimal() +
  # geom_text(aes(x = 1.2, y = Inf, 
  #               label = paste0("p-value = ", signif(pval, digits = 2)),
  #               vjust = 1)) +
  geom_text(
    data = pvals_df,
    aes(x = 1.2,              # center between groups
        y = Inf,
        label = paste0("p = ", signif(pval, 2))),
    inherit.aes = FALSE,
    vjust = 1.2
  ) + 
  labs(x = "", y = "Normalized counts", fill = "DISEASE") +
  theme(plot.title = element_text(hjust = 0.5, size = 26),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 40, vjust = 0.5),  
        axis.text.y = element_text(size = 15))

# png(filename = glue::glue("{DIR}/figures_updated/violinplot_interferonSign.png"),
    # height = 20, width = 25, units = "cm", res = 500)
print(p_inferferon)
# dev.off()

ggsave(filename = glue::glue("{DIR}/figures_updated/violinplot_interferonSign.svg"),
       height = 20, width = 25, units = "cm", dpi = 500, plot = p_inferferon)
######
#######

# Table of all previously associated genes
# table genes
list.genes_prev <- c("ATG5", "BANK1","BLK","C1QA","C1QB","C1QC","C2", "C4A","C4B",
                    "CRP","ETS1","FAM167A", "FCGR2A", 
                    "GTF2I", "GTF2IRD1","HIP1","HLA-DQA1","HLA-DQB1",
                    "HLA-DRB1", "IKZF1", "IL12RB2","IRAK1","IRF5", 
                    "ITGAM","JAZF1","KCP","LBH", "LYN",
                    "PHRF1", "PTPN22","RASGRP3","SLC15A4", 
                    "SPATA48","STAT4","TNFAIP3","TNFSF4","TNIP1", "TNPO3","TREX1",
                    "TYK2","UBE2L3","WDFY4")
all_pvalues_prev <- data.frame()


for (gene in list.genes_prev){
  if(gene %in% norm_counts_name$gene_name){
    counts.gene <- norm_counts_name %>% 
    filter(gene_name == gene) %>%
    dplyr::select(-gene_name)
    counts.gene<- as.data.frame(t(counts.gene))
  
    # test
    expression = counts.gene[metadata$samples,1]
    sample = metadata$DISEASE
    df.gene <- data.frame(expression, sample)
    
    control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
    sle.e<- df.gene[df.gene$sample=='SLE',]$expression
    test <- permTS(control.e, sle.e)
    pval <- (test$p.value)
    print(gene)
    o <- data.frame(gene_name = gene, pval = pval )
    all_pvalues_prev <- rbind(all_pvalues_prev,o)
  }
}

#####
# Adjust p values


all_pvalues_prev$padj <- p.adjust(all_pvalues_prev$pval, method = "fdr", n = nrow(all_pvalues_prev))
length(adj_pvalues)


write.csv(all_pvalues_prev, file = glue::glue("{DIR}/data/FIG2A-genes-pvals.csv"),
          row.names = F, quote = F)

# ------------------------
# SELECTED PREVIOUSLY ASSOCIATED GENES

# Complement: C2, C1QC, C1QB, C1QA
# Removal of immune complexes: FCGR2A
# Inflammation: CRP, HLA-DQA1, TNFAIP3
# Autophagy: ATG5

violin_prev_df <- norm_counts_name %>%
  filter(gene_name %in% c("C2", "C1QC", "C1QB", "C1QA",
                          "FCGR2A", "TREX1", "HLA-DQA1","TNFAIP3",
                          "ATG5")) %>%
  pivot_longer(names_to = "sample", cols = -gene_name) %>%
  left_join(metadata, by = c("sample" = "samples")) %>%
  left_join(all_pvalues_prev, by = c("gene_name" = "gene_name"))

sub_pvalues_prev <- all_pvalues_prev %>%
  filter(gene_name %in% violin_prev_df$gene_name)

df_mean <- violin_prev_df %>%
  group_by(DISEASE, gene_name) %>%
  summarize(average = mean(value)) %>%
  ungroup()
violin_prev_df <- violin_prev_df %>%
  left_join(df_mean, by = c("gene_name" = "gene_name", "DISEASE"="DISEASE")) %>%
  arrange(pval)

violin_prev_df$gene_name <- factor(violin_prev_df$gene_name, levels = unique(violin_prev_df$gene_name))

p_prev <- ggplot(violin_prev_df, aes(x = DISEASE, y = value, fill = DISEASE)) +
  geom_violin(trim=FALSE) +
  facet_wrap(~as.factor(gene_name), ncol = 3, scales = "free_y") +
  geom_boxplot(width=0.2, color = 'black', fill=NA) +
  geom_point(df_mean, mapping = aes(x = DISEASE, y = average),
             color = 'white', shape = 18, size = 3) +
  scale_fill_manual(values=c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')) +
  geom_point(aes(x = DISEASE, y = average), color = 'white', shape = 18, size = 3) +
  geom_line(aes(x = DISEASE, y = average, group = 1), linetype = "dashed") +
  theme_minimal() +
  # geom_text(aes(x = 1.2, y = Inf, 
  #               label = paste0("p-value = ", signif(pval, digits = 2)),
  #               vjust = 1)) +
  geom_text(
    data = sub_pvalues_prev,
    aes(x = 1.2,              # center between groups
        y = Inf,
        label = paste0("p = ", signif(padj, 2))),
    inherit.aes = FALSE,
    vjust = 1.2
  ) + 
  labs(x = "", y = "Normalized counts", fill = "DISEASE") +
  theme(plot.title = element_text(hjust = 0.5, size = 26),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 40, vjust = 0.5),  
        axis.text.y = element_text(size = 15))

# png(filename = glue::glue("{DIR}/figures_updated/violinplot_prevAssociated.png"),
# height = 20, width = 25, units = "cm", res = 500)
print(p_prev)
# dev.off()

ggsave(filename = glue::glue("{DIR}/figures_updated/violinplot_prevAssociated.svg"),
       height = 20, width = 25, units = "cm", dpi = 500, plot = p_prev)

#############
## Reporductive genes
genes_reproductive <- c("F10", "F11", "F12", "F13A1", "F13B", "F2", "F5", "F7",
                        "F8", "F9", "FGA", "FGB", "FGG", "GP1BA", "HRG", "KLKB1",
                        "KNG1", "LMAN1", "MCFD2", "MTHFR", "PLAT", "PROC", "PROS1",
                        "PROZ", "SERPINC1", "SERPIND1", "SERPINE1", "SERPINF2",
                        "THBD", "VKORC1", "VWF")

significant_reproductive <- c("F5", "F2","F12","F8","VKORC1","PROZ","SERPINF2","F10")

pvalues_reproductive <- data.frame()


for (gene in significant_reproductive){
  if(gene %in% norm_counts_name$gene_name){
    counts.gene <- norm_counts_name %>% 
      filter(gene_name == gene) %>%
      dplyr::select(-gene_name)
    counts.gene<- as.data.frame(t(counts.gene))
    
    # test
    expression = counts.gene[metadata$samples,1]
    sample = metadata$DISEASE
    df.gene <- data.frame(expression, sample)
    
    control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
    sle.e<- df.gene[df.gene$sample=='SLE',]$expression
    test <- permTS(control.e, sle.e)
    pval <- (test$p.value)
    print(gene)
    o <- data.frame(gene_name = gene, pval = pval )
    pvalues_reproductive <- rbind(pvalues_reproductive,o)
  }
}
pvalues_reproductive$padj <- p.adjust(pvalues_reproductive$pval, method = "fdr", n = nrow(pvalues_reproductive))

violin_reproductive_df <- norm_counts_name %>%
  filter(gene_name %in% pvalues_reproductive$gene_name) %>%
  pivot_longer(names_to = "sample", cols = -gene_name) %>%
  left_join(metadata, by = c("sample" = "samples")) %>%
  left_join(pvalues_reproductive, by = c("gene_name" = "gene_name"))


df_mean <- violin_reproductive_df %>%
  group_by(DISEASE, gene_name) %>%
  summarize(average = mean(value)) %>%
  ungroup()
violin_reproductive_df <- violin_reproductive_df %>%
  left_join(df_mean, by = c("gene_name" = "gene_name", "DISEASE"="DISEASE")) %>%
  arrange(pval)

violin_reproductive_df$gene_name <- factor(violin_reproductive_df$gene_name, levels = unique(violin_reproductive_df$gene_name))

p_reproductive <- ggplot(violin_reproductive_df, aes(x = DISEASE, y = value, fill = DISEASE)) +
  geom_violin(trim=FALSE) +
  facet_wrap(~as.factor(gene_name), ncol = 4, scales = "free_y") +
  geom_boxplot(width=0.2, color = 'black', fill=NA) +
  geom_point(df_mean, mapping = aes(x = DISEASE, y = average),
             color = 'white', shape = 18, size = 3) +
  scale_fill_manual(values=c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')) +
  geom_point(aes(x = DISEASE, y = average), color = 'white', shape = 18, size = 3) +
  geom_line(aes(x = DISEASE, y = average, group = 1), linetype = "dashed") +
  theme_minimal() +
  # geom_text(aes(x = 1.2, y = Inf, 
  #               label = paste0("p-value = ", signif(pval, digits = 2)),
  #               vjust = 1)) +
  geom_text(
    data = pvalues_reproductive,
    aes(x = 1.2,              # center between groups
        y = Inf,
        label = paste0("p = ", signif(padj, 2))),
    inherit.aes = FALSE,
    vjust = 1.2
  ) + 
  labs(x = "", y = "Normalized counts", fill = "DISEASE") +
  theme(plot.title = element_text(hjust = 0.5, size = 26),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 40, vjust = 0.5),  
        axis.text.y = element_text(size = 15))

png(filename = glue::glue("{DIR}/figures_updated/violinplot_reproductive.png"),
height = 18, width = 25, units = "cm", res = 500)
print(p_reproductive)
dev.off()
# ----
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
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] perm_1.0-0.4                DESeq2_1.48.2               SummarizedExperiment_1.38.1
# [4] Biobase_2.68.0              MatrixGenerics_1.20.0       matrixStats_1.5.0          
# [7] GenomicRanges_1.60.0        GenomeInfoDb_1.44.3         IRanges_2.42.0             
# [10] S4Vectors_0.46.0            BiocGenerics_0.54.1         generics_0.1.4             
# [13] ComplexHeatmap_2.25.2       lubridate_1.9.5             forcats_1.0.1              
# [16] stringr_1.6.0               dplyr_1.2.0                 purrr_1.2.1                
# [19] readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1               
# [22] ggplot2_4.0.2               tidyverse_2.0.0             readxl_1.4.5               
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        farver_2.1.2            S7_0.2.1               
# [4] digest_0.6.39           timechange_0.4.0        lifecycle_1.0.5        
# [7] cluster_2.1.8.2         magrittr_2.0.4          compiler_4.5.0         
# [10] rlang_1.1.7             tools_4.5.0             utf8_1.2.6             
# [13] data.table_1.18.2.1     knitr_1.51              S4Arrays_1.8.1         
# [16] labeling_0.4.3          DelayedArray_0.34.1     RColorBrewer_1.1-3     
# [19] pkgload_1.5.0           abind_1.4-8             BiocParallel_1.42.2    
# [22] withr_3.0.2             colorspace_2.1-2        scales_1.4.0           
# [25] iterators_1.0.14        dichromat_2.0-0.1       cli_3.6.5              
# [28] crayon_1.5.3            otel_0.2.0              rstudioapi_0.18.0      
# [31] httr_1.4.8              tzdb_0.5.0              rjson_0.2.23           
# [34] parallel_4.5.0          cellranger_1.1.0        XVector_0.48.0         
# [37] vctrs_0.7.1             Matrix_1.7-4            jsonlite_2.0.0         
# [40] hms_1.1.4               GetoptLong_1.1.0        clue_0.3-67            
# [43] locfit_1.5-9.12         foreach_1.5.2           glue_1.8.0             
# [46] codetools_0.2-20        stringi_1.8.7           shape_1.4.6.1          
# [49] gtable_0.3.6            UCSC.utils_1.4.0        pillar_1.11.1          
# [52] GenomeInfoDbData_1.2.14 circlize_0.4.17         R6_2.6.1               
# [55] doParallel_1.0.17       evaluate_1.0.5          lattice_0.22-9         
# [58] png_0.1-8               Rcpp_1.1.1              SparseArray_1.8.1      
# [61] xfun_0.56               pkgconfig_2.0.3         GlobalOptions_0.1.3    