# Violin Plots
#####
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(perm)
load('/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData')
load("/mnt/Citosina/amedina/ssalazar/meta/combined/vsd2.RData")
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

###### 

## In local

norm_counts_name <- read.csv('/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/normcounts_name.csv')
all_data <-  read.csv('/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/all_data.csv')

rownames(norm_counts_name) <- norm_counts_name$X
norm_counts_name<-norm_counts_name[,-c(1)]
all_pvalues <- list()
list.genes <- c("IFI27", "OTOF", "IFI44L","SIGLEC1","USP18", "IFI44", "IFIT1", "SPATS2L")
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
#######


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
#####


# VIOLIN PLOTS PANEL A

for (i in 1:length(list.genes)){
  gene <- list.genes[i]
  counts.gene <- norm_counts_name[rownames(norm_counts_name)==gene,]
  counts.gene<- as.data.frame(t(counts.gene))
  expression = counts.gene[,1]
  sample = all_data$DISEASE
  df.gene <- data.frame(expression, sample)
  
  control.e <- df.gene[df.gene$sample=='CONTROL',]$expression
  sle.e<- df.gene[df.gene$sample=='SLE',]$expression
  
  pval <- adj_pvalues[i]
  if(pval < 0.05){
    print(gene)
    # plot
    df_mean <- df.gene %>%
      group_by(sample) %>%
      summarize(average = mean(expression)) %>%
      ungroup()
  
    p <- ggplot(df.gene, aes(x=sample, y=expression, fill=sample)) + 
      geom_violin(trim=FALSE) +
      labs(title=gene,x = NULL, y="Normalized counts") +
      ylim(0,21)
  
    p1 <- p + geom_boxplot(width=0.15, color = 'black', fill=NA) +
      # geom_jitter(shape=16, position=position_jitter(0.2)) +
      scale_fill_manual(values=c("#e65a5a", "#4d80c4")) +
      geom_point(df_mean, mapping = aes(x = sample, y = average), color = 'white', shape = 18, size = 3) +
      geom_line(df_mean, mapping = aes(x = sample, y = average, group = 1), linetype = "dashed")  +
      theme_classic()
  
    p2 <- p1 + theme(legend.position="none") +
      annotate("text",
             x = 1:length(table(df.gene$sample)),
             y = 19,
             label = paste0('n = ', table(df.gene$sample)),
             col = "black",
             vjust = - 1)
  
    p3 <- p2 + annotate("text",
                      x = 1,
                      y = 17,
                      label = paste0('p-value = ', signif(pval, digits = 4)),
                      col = "black",
                      vjust = - 1)
  
  # ggsave(paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot",gene,"-ptest.png"), dpi = 300, plot = p3)  
    save(p3, file = paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-",gene,'.RData'))
  }

}
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

list.genes_new <- head(table %>% arrange(p.value),8)$Gene
arr_pvalues <- head(table %>% arrange(p.value),8)$p.value
####
# Volcano plots PANEL B
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
      labs(title=gene,x = NULL, y="Normalized counts") +
      ylim(0,17)
    
    p1 <- p + geom_boxplot(width=0.3, color = 'black', fill=NA) +
      # geom_jitter(shape=16, position=position_jitter(0.2)) +
      scale_fill_manual(values=c("#e65a5a", "#4d80c4")) +
      geom_point(df_mean, mapping = aes(x = sample, y = average), color = 'white', shape = 18, size = 3) +
      geom_line(df_mean, mapping = aes(x = sample, y = average, group = 1), linetype = "dashed")  +
      theme_classic()
    
    p2 <- p1 + theme(legend.position="none") +
      annotate("text",
               x = 1:length(table(df.gene$sample)),
               y = 16,
               label = paste0('n = ', table(df.gene$sample)),
               col = "black",
               vjust = - 1)
    
    p3 <- p2 + annotate("text",
                        x = 1,
                        y = 14.5,
                        label = paste0('p-value = ', signif(pvalue, digits = 4)),
                        col = "black",
                        vjust = - 1)
    
    # ggsave(paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot",gene,"-ptest.png"), dpi = 300, plot = p3)  
    save(p3, file = paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-",gene,'.RData'))
  }
}

