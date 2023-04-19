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

# Boxplots
norm_counts_name <- read.csv('/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/normcounts_name.csv')
all_data <-  read.csv('/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/all_data.csv')

rownames(norm_counts_name) <- norm_counts_name$X
norm_counts_name<-norm_counts_name[,-c(1)]

list.genes <- c("IFI27","BANK1","BLK","C1QA","C1QB","C1QC","C2","FCGR2A","FCGR3A","IRF5","ITGAM","LYN","STAT4","SLC1A7","TREX1","DNASE1L3","TLR7","TLR9","IFI44L","IL12RB2",'HLA-A','HLA-B','HLA-C')

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
  
  # wilcox.test(x = control.e, y = sle.e, alternative = "two.sided", mu = 0,
  # paired = FALSE, conf.int = 0.95)
  if(test$p.value < 0.05){
    print(gene)
    # plot
    
    p <- ggplot(df.gene, aes(x=sample, y=expression, fill=sample)) + 
      geom_violin(trim=FALSE) +
      labs(title=gene,x = NULL, y="Normalized counts") +
      ylim(0,18)
    
    p1 <- p + geom_boxplot(width=0.05, color = 'black', fill='gray') +
      # geom_jitter(shape=16, position=position_jitter(0.2)) +
      scale_fill_manual(values=c("#e65a5a", "#4d80c4")) +
      stat_summary(fun.y=mean, geom = 'point', shape = 23, size = 2) +
      theme_classic()
    
    p2 <- p1 + theme(legend.position="none") +
      annotate("text",
               x = 1:length(table(df.gene$sample)),
               y = 17,
               label = paste0('n = ', table(df.gene$sample)),
               col = "black",
               vjust = - 1)
    
    p3 <- p2 + annotate("text",
                        x = 1,
                        y = 15,
                        label = paste0('Perm test, p = ', signif(test$p.value, digits = 4)),
                        col = "black",
                        vjust = - 1)
    
    ggsave(paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot",gene,"-ptest.png"), dpi = 300, plot = p3)  
    
  }
}