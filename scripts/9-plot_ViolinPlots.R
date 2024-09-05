# Violin Plots
#####
library(ggplot2)
library(DESeq2)
library(tidyverse)
load('/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData')
load("/mnt/Citosina/amedina/backup/lupus/sofi/vsd2.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")
outdir = '/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/'
#####

#####
# In cluster

head(df_names)

# name norm matrix with gene name
norm_counts <- as.data.frame(assay(vsd2))
norm_counts <- tibble::rownames_to_column(norm_counts, "ID")
name_short<-function(table){
  table$ID <- gsub("\\..*","", table$ID)
  return(table)
}

norm_counts_id <- name_short(norm_counts)

norm_counts_name <- inner_join(norm_counts_id, df_names, by = 'ID')

norm_counts_name<- norm_counts_name[,-c(1,320:327)]
rownames(norm_counts_name) <- df_names$gene_name

write.csv(norm_counts_name, "/mnt/Citosina/amedina/ssalazar/meta/combined/normcounts_name.csv")

#################################################
#################################################
#################################################

# ------------------In local----------------------

library(tidyverse)
library(perm)
norm_counts_name <- read.csv('/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/normcounts_name.csv',
                             header = T, row.names = "X")
all_data <-  read.csv('/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/all_data.csv')


##########
# OUR SELECTED INTERFERON SIGNATURE GENES

# compute pvalues
all_pvalues <- list()
list.genes <- c("CCL2", "IFIT1", "RSAD2", "IFI44L", 
                "IFI44", "USP18", "ISG15", "IFI27", "SIGLEC1")
for (i in 1:length(list.genes)){
  gene <- list.genes[i]
  counts.gene <- norm_counts_name[rownames(norm_counts_name)==gene,]
  counts.gene<- as.data.frame(t(counts.gene))
  
  # test
  expression = counts.gene[,1]
  sample = all_data$DISEASE
  df.gene <- data.frame(expression, sample)
  
  control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
  sle.e<- df.gene[df.gene$sample=='SLE',]$expression
  test <- permTS(control.e, sle.e)
  pval <- (test$p.value)
  print(gene)
  all_pvalues[i] <- pval
  
}

l <- length(all_pvalues)
l

# Adjust p values


adj_pvalues <- p.adjust(unlist(all_pvalues), method = "fdr", n = length(all_pvalues))
length(adj_pvalues)

pvals_df <- data.frame(pval = adj_pvalues, gene = list.genes)
pvals_df <- pvals_df %>% arrange(pval)
for (i in 1:length(pvals_df$pval)){
  gene <- pvals_df$gene[i]
  counts.gene <- norm_counts_name[rownames(norm_counts_name)==gene,]
  counts.gene<- as.data.frame(t(counts.gene))
  expression = counts.gene[,1]
  sample = all_data$DISEASE
  check_data <- all(all_data$samples == rownames(counts.gene))
  if (check_data){
    df.gene <- data.frame(expression, sample)
    
    control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
    sle.e<- df.gene[df.gene$sample=='SLE',]$expression
    
    pval <- pvals_df$pval[i]
    if(pval < 0.05){
      print(gene)
      # plot
      df_mean <- df.gene %>%
        group_by(sample) %>%
        summarize(average = mean(expression)) %>%
        ungroup()
      
      p <- ggplot(df.gene, aes(x=sample, y=expression, fill=sample)) + 
        geom_violin(trim=FALSE) +
        labs(title=gene,x = NULL, y="Normalized counts", fill = "Samples") +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        ylim(0,17)
      
      p1 <- p + geom_boxplot(width=0.2, color = 'black', fill=NA) +
        # geom_jitter(shape=16, position=position_jitter(0.2)) +
        scale_fill_manual(values=c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')) +
        geom_point(df_mean, mapping = aes(x = sample, y = average), color = 'white', shape = 18, size = 3) +
        geom_line(df_mean, mapping = aes(x = sample, y = average, group = 1), linetype = "dashed")  +
        theme_classic()
      
      # p2 <- p1 + theme(legend.position="none") +
      #   annotate("text",
      #            x = 1:length(table(df.gene$sample)),
      #            y = 16,
      #            label = paste0('n = ', table(df.gene$sample)),
      #            col = "black",
      #            vjust = - 1)
      
      p3 <- p1 + annotate("text",
                          x = 1,
                          y = 14.5,
                          label = paste0('p-value = ', signif(pval, digits = 2)),
                          col = "black",
                          vjust = - 1,
                          size = 6)
      p3 <- p3 + theme(plot.title = element_text(hjust = 0.5, size = 18),
                       plot.background = element_rect(fill = "white"),
                       text = element_text(size = 16),
                       axis.title = element_text(size = 16),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 16),
                       axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5),  
                       axis.text.y = element_text(size = 16))
      # ggsave(paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot",gene,"-ptest.png"), dpi = 300, plot = p3)  
      save(p3, file = paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violin/violinplot-",gene,'.RData'))
    }
  } else{
    print("Found discrepancy")
  }
}
######
#######

# Table of all previously associated genes
# table genes
list.genes_tab <- c("ATG5", "BANK1","BLK","C1QA","C1QB","C1QC","C2", "C4A","C4B",
                    "CRP","ETS1","FAM167A", "FCGR2A", 
                    "GTF2I", "GTF2IRD1","HIP1","HLA-DQA1","HLA-DQB1",
                    "HLA-DRB1", "IKZF1", "IL12RB2","IRAK1","IRF5", 
                    "ITGAM","JAZF1","KCP","LBH", "LYN",
                    "PHRF1", "PTPN22","RASGRP3","SLC15A4", 
                    "SPATA48","STAT4","TNFAIP3","TNFSF4","TNIP1", "TNPO3","TREX1",
                    "TYK2","UBE2L3","WDFY4")


for (i in 1:length(list.genes_tab)){
  gene <- list.genes_tab[i]
  counts.gene <- norm_counts_name[rownames(norm_counts_name)==gene,]
  counts.gene<- as.data.frame(t(counts.gene))
  
  # test
  expression = counts.gene[,1]
  
  if(!is.null(expression)){
  sample = all_data$DISEASE
  df.gene <- data.frame(expression, sample)
  
  control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
  sle.e<- df.gene[df.gene$sample=='SLE',]$expression
  test <- permTS(control.e, sle.e)
  p <- (test$p.value)
  all_pvalues[l+i] <- p
  print(gene)
  }
  
}
length(all_pvalues)
#####
# Adjust p values

adj_pvalues <- p.adjust(unlist(all_pvalues), method = "fdr", n = length(all_pvalues))
length(adj_pvalues)
######


g <- list()
pvals <- list()
for (i in 1:length(list.genes_tab)){
  gene <- list.genes_tab[i]
  p <- adj_pvalues[length(list.genes)+i]
  if(p < 0.05){
    print(gene)
    g[i] <- gene
    pvals[i] <- p
  }
}

g <- unlist(g)
pvals <- unlist(pvals)

table <- data.frame(as.character(g), as.numeric(pvals))
colnames(table) <- c("Gene", "p.value")
write.csv(table, file = "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/FIG2A-genes-pvals.csv")

# ------------------------
# SELECTED PREVIOUSLY ASSOCIATED GENES

# Complement: C2, C1QC, C1QB, C1QA
# Removal of immune complexes: FCGR2A
# Inflammation: CRP, HLA-DQA1, TNFAIP3
# Autophagy: ATG5

table <- read.csv(file = "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/FIG2A-genes-pvals.csv",
                  header = T, row.names = "X")

prev_asso_df <- table %>% filter(Gene %in% c("C2", "C1QC", "C1QB", "C1QA",
                                             "FCGR2A", "CRP", "HLA-DQA1","TNFAIP3",
                                             "ATG5")) %>% arrange(p.value)

list.genes_new <- prev_asso_df$Gene
arr_pvalues <- prev_asso_df$p.value

for (i in 1:length(list.genes_new)){
  gene <- list.genes_new[i]
  counts.gene <- norm_counts_name[rownames(norm_counts_name)==gene,]
  counts.gene<- as.data.frame(t(counts.gene))
  
  expression = counts.gene[,1]
  sample = all_data$DISEASE
  df.gene <- data.frame(expression, sample)
  
  control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
  sle.e<- df.gene[df.gene$sample=='SLE',]$expression
  pvalue <- arr_pvalues[i]
  if(pvalue < 0.05){
    print(gene)
    # plot
    df_mean <- df.gene %>%
      group_by(sample) %>%
      summarize(average = mean(expression)) %>%
      ungroup()
    p <- ggplot(df.gene, aes(x=sample, y=expression, fill=sample)) + 
      geom_violin(trim=FALSE) +
      labs(title=gene,x = NULL, y="Normalized counts", fill = "Samples") +
      theme(plot.title = element_text(hjust = 0.5, size = 20)) +
      ylim(0,17)
    
    p1 <- p + geom_boxplot(width=0.3, color = 'black', fill=NA) +
      # geom_jitter(shape=16, position=position_jitter(0.2)) +
      scale_fill_manual(values=c('CONTROL' = '#96d4ccff', 'SLE' = '#b493b4ff')) +
      geom_point(df_mean, mapping = aes(x = sample, y = average), color = 'white', shape = 18, size = 3) +
      geom_line(df_mean, mapping = aes(x = sample, y = average, group = 1), linetype = "dashed")  +
      theme_classic()
    
    # p2 <- p1 + theme(legend.position="none") +
    #   annotate("text",
    #            x = 1:length(table(df.gene$sample)),
    #            y = 16,
    #            label = paste0('n = ', table(df.gene$sample)),
    #            col = "black",
    #            vjust = - 1)
    
    p3 <- p1 + annotate("text",
                        x = 1,
                        y = 14.5,
                        label = paste0('p-value = ', signif(pvalue, digits = 2)),
                        col = "black",
                        vjust = - 1,
                        size = 6)
    p3 <- p3 + theme(plot.title = element_text(hjust = 0.5, size = 18),
                     plot.background = element_rect(fill = "white"),
                     text = element_text(size = 16),
                     axis.title = element_text(size = 16),
                     legend.title = element_text(size = 16),
                     legend.text = element_text(size = 16),
                     axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5),  
                     axis.text.y = element_text(size = 16))
    # ggsave(paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot",gene,"-ptest.png"), dpi = 300, plot = p3)  
    save(p3, file = paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violin/violinplot-",gene,'.RData'))
  }
}

p3
######
sessionInfo()

# R version 4.2.2 (2022-10-31)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.2.1

# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#   [1] perm_1.0-0.2                lubridate_1.9.2            
# [3] forcats_1.0.0               stringr_1.5.0              
# [5] dplyr_1.1.2                 purrr_1.0.1                
# [7] readr_2.1.4                 tidyr_1.3.0                
# [9] tibble_3.2.1                tidyverse_2.0.0            
# [11] DESeq2_1.36.0               SummarizedExperiment_1.26.1
# [13] Biobase_2.56.0              MatrixGenerics_1.8.1       
# [15] matrixStats_1.0.0           GenomicRanges_1.48.0       
# [17] GenomeInfoDb_1.32.4         IRanges_2.30.1             
# [19] S4Vectors_0.34.0            BiocGenerics_0.44.0        
# [21] ggplot2_3.3.5              

# loaded via a namespace (and not attached):
#   [1] httr_1.4.6             bit64_4.0.5            splines_4.2.2         
# [4] blob_1.2.4             GenomeInfoDbData_1.2.8 yaml_2.3.7            
# [7] pillar_1.9.0           RSQLite_2.3.1          lattice_0.21-8        
# [10] glue_1.6.2             limma_3.52.4           digest_0.6.32         
# [13] RColorBrewer_1.1-3     XVector_0.36.0         colorspace_2.1-0      
# [16] htmltools_0.5.5        Matrix_1.5-4.1         XML_3.99-0.14         
# [19] pkgconfig_2.0.3        genefilter_1.78.0      zlibbioc_1.42.0       
# [22] xtable_1.8-4           scales_1.2.1           tzdb_0.4.0            
# [25] BiocParallel_1.30.4    timechange_0.2.0       annotate_1.74.0       
# [28] KEGGREST_1.36.3        generics_0.1.3         cachem_1.0.8          
# [31] withr_2.5.0            cli_3.6.1              survival_3.5-5        
# [34] magrittr_2.0.3         crayon_1.5.2           memoise_2.0.1         
# [37] evaluate_0.21          fansi_1.0.4            tools_4.2.2           
# [40] data.table_1.14.8      hms_1.1.3              lifecycle_1.0.3       
# [43] munsell_0.5.0          locfit_1.5-9.8         DelayedArray_0.22.0   
# [46] AnnotationDbi_1.58.0   Biostrings_2.64.1      compiler_4.2.2        
# [49] rlang_1.1.1            grid_4.2.2             RCurl_1.98-1.12       
# [52] rstudioapi_0.14        bitops_1.0-7           rmarkdown_2.22        
# [55] gtable_0.3.3           codetools_0.2-19       DBI_1.1.3             
# [58] R6_2.5.1               knitr_1.43             fastmap_1.1.1         
# [61] bit_4.0.5              utf8_1.2.3             stringi_1.7.12        
# [64] parallel_4.2.2         Rcpp_1.0.11            vctrs_0.6.3           
# [67] geneplotter_1.74.0     png_0.1-8              tidyselect_1.2.0      
# [70] xfun_0.39    