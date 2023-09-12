# Volcano plot

#####

library(ggplot2)
library(tidyverse)

load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"

#####


dge <- dge[order(dge$log2FoldChange),]
colnames(dge) <- c("GeneName", "log2FC", "pvalAdj", "Expression")
write.csv(dge, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/DGE_list.csv")

data <- df_names %>% 
  mutate(Expression = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))

volcanoplot <- ggplot(data, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-adj")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

data %>% 
  count(Expression) %>% 
  knitr::kable()


top <- 20 # no. of highlighted genes in plot
top_genes <- bind_rows(data %>% 
                         filter(Expression == 'Up-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top),
                       data %>% 
                         filter(Expression == 'Down-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top)
)

dim(data[data$Expression == 'Up-regulated',])
dim(data[data$Expression == 'Down-regulated',])

top_genes %>% 
  knitr::kable()

volcanoplot_names <-  volcanoplot +
  ggrepel::geom_label_repel(data = top_genes,
                            mapping = aes(log2FoldChange, -log(padj,10), label = gene_name),
                            size = 2) + theme_classic() + theme(legend.position = 'bottom')

save(volcanoplot_names, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/volcanoplot_object.RData")

ggsave(paste0(outdir,"volcanoPlotWithTopGenes.png"),
       plot = volcanoplot_names, dpi = 300)
dev.off()

######
sessionInfo()

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS:   /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRblas.so
# LAPACK: /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRlapack.so

# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.10    purrr_0.3.4    
# [5] readr_1.4.0     tidyr_1.2.1     tibble_3.1.3    tidyverse_1.3.0
# [9] ggplot2_3.3.5  

# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.7        cellranger_1.1.0  pillar_1.6.2      compiler_4.0.2   
# [5] dbplyr_2.2.1      tools_4.0.2       jsonlite_1.7.2    lubridate_1.7.9.2
# [9] lifecycle_1.0.3   gtable_0.3.0      pkgconfig_2.0.3   rlang_1.0.6      
# [13] reprex_1.0.0      rstudioapi_0.13   DBI_1.1.1         cli_3.6.0        
# [17] haven_2.3.1       xml2_1.3.2        withr_2.4.2       httr_1.4.2       
# [21] fs_1.5.0          generics_0.1.0    vctrs_0.5.1       hms_1.1.0        
# [25] grid_4.0.2        tidyselect_1.2.0  glue_1.4.2        R6_2.5.0         
# [29] fansi_0.5.0       readxl_1.3.1      modelr_0.1.8      magrittr_2.0.1   
# [33] scales_1.1.1      backports_1.2.1   ellipsis_0.3.2    rvest_0.3.6      
# [37] assertthat_0.2.1  colorspace_2.0-2  utf8_1.2.2        stringi_1.6.2    
# [41] munsell_0.5.0     broom_0.7.9       crayon_1.4.1 