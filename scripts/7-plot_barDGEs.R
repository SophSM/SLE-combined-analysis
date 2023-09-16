# Barplots for all DEGs
#####
library(tidyverse)
library(ggplot2)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"
DGE.list <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/DGElist_withNames_noncoding.csv")
DGE.list <- DGE.list[,-1]
#####


upreg <- na.omit(subset(DGE.list, padj < 0.05 & (log2FoldChange > 1)))
downreg <- na.omit(subset(DGE.list, padj < 0.05 & (log2FoldChange < -1 )))

levels(as.factor(upreg$transcript_biotype))

count.types1 <- as.data.frame(table(upreg$transcript_biotype))
keep1 <- c("lncRNA","miRNA","protein_coding","protein_coding_CDS_not_defined","snoRNA","snRNA")
df1.small <-count.types1[count.types1$Var1 %in% keep1,]


levels(as.factor(downreg$transcript_biotype))
count.types2 <- as.data.frame(table(downreg$transcript_biotype))
keep2 <- c("lncRNA","protein_coding","protein_coding_CDS_not_defined")
df2.small <-count.types2[count.types2$Var1 %in% keep2,]


## barplots

p <- ggplot(data = df1.small, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', fill = '#c21f35') + 
  geom_text(aes(label = Freq, x = Var1, y = Freq, hjust = -0.5, size = 2)) +
  coord_flip() +
  xlab(NULL) + 
  theme_classic() +
  theme(legend.position="none",
        axis.ticks.y=element_blank())

outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"
save(p, file=paste0(outdir,"barplot_up_nc.RData"))
# ggsave(paste0(outdir, "biotypes_barplot.png"), plot = p, dpi = 300)

p2 <- ggplot(data = df2.small, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', fill = '#195ae6') + 
  geom_text(aes(label = Freq, x = Var1, y = Freq, hjust = -0.5, size = 2)) +
  coord_flip() +
  xlab(NULL) +
  theme_classic() +
  theme(legend.position="none",
        axis.ticks.y=element_blank())
# ggsave(paste0(outdir, "biotypes_barplot_2.png"), plot = p2, dpi = 300)

save(p2, file=paste0(outdir,"barplot_down_nc.RData"))

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
# [5] readr_1.4.0     tidyr_1.2.1     tibble_3.1.3    ggplot2_3.3.5  
# [9] tidyverse_1.3.0

# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.7        cellranger_1.1.0  pillar_1.6.2      compiler_4.0.2   
# [5] dbplyr_2.2.1      tools_4.0.2       jsonlite_1.7.2    lubridate_1.7.9.2
# [9] lifecycle_1.0.3   gtable_0.3.0      pkgconfig_2.0.3   rlang_1.0.6      
# [13] reprex_1.0.0      rstudioapi_0.13   DBI_1.1.1         cli_3.6.0        
# [17] haven_2.3.1       withr_2.4.2       xml2_1.3.2        httr_1.4.2       
# [21] fs_1.5.0          generics_0.1.0    vctrs_0.5.1       hms_1.1.0        
# [25] grid_4.0.2        tidyselect_1.2.0  glue_1.4.2        R6_2.5.0         
# [29] fansi_0.5.0       readxl_1.3.1      modelr_0.1.8      magrittr_2.0.1   
# [33] backports_1.2.1   scales_1.1.1      ellipsis_0.3.2    rvest_0.3.6      
# [37] assertthat_0.2.1  colorspace_2.0-2  utf8_1.2.2        stringi_1.6.2    
# [41] munsell_0.5.0     broom_0.7.9       crayon_1.4.1