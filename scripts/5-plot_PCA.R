# PCA
#####
library(DESeq2)
library(ggplot2)
library(tidyverse)
load('/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData')
load("/mnt/Citosina/amedina/backup/lupus/sofi/vsd2.RData")
outdir = '/mnt/Citosina/amedina/ssalazar/meta/combined/figures/'

######

mat <- as.matrix(assay(vsd2)) # col = samples, rows  = genes
pc <- prcomp(t(mat)) # col = genes , rows = samples
dtp <- data.frame('DISEASE' = all_data$DISEASE, pc$x[,1:2])

pcaData <- plotPCA(vsd2, 'DISEASE', returnData = TRUE)

pcaData <- tibble::rownames_to_column(pcaData, "study")

percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=group, shape = group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c('CONTROL' = '#a9e536', 'SLE' = '#f5704b'))+
  stat_ellipse() + theme_classic() +
  labs(color = 'Samples', shape = "Samples") +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"))

save(pca_plot, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/pca_object.RData")

ggsave(paste0(outdir,"PCA-ellipse.png"), width = 3000, height = 3000,
       units = 'px', dpi = 300, bg = "white", plot = pca_plot)

# Color by study
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

pcaData <- plotPCA(vsd2, 'DISEASE', returnData = TRUE)
pcaData <- tibble::rownames_to_column(pcaData, "samples")


pc_df <- pcaData %>%
  inner_join(all_data, by = c("samples" = "samples"))
percentVar <- round(100 * attr(pcaData, "percentVar"))


pca_plot2 <- ggplot(pc_df, aes(PC1, PC2, color=study, shape = DISEASE.x)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c('GSE175839' = '#f5a142', 'GSE122459' = '#f5ef42', 
                                'GSE169080' = '#2ef0e9', 'GSE162828' = '#f02eb3', 
                                'GSE101437' = '#a1645c', 'GSE112087'='#599163', 
                                'GSE80183'='#755c91','GSE72509'='#e68a8a')) +
  stat_ellipse(aes(group = DISEASE.x)) + theme_classic() +
  labs(color = 'Study', shape = "Samples") +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"))

ggsave(paste0(outdir,"PCA-studies.png"), width = 3000, height = 3000,
       units = 'px', dpi = 300, bg = "white", plot = pca_plot2)
#######
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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#   [1] ggplot2_3.3.5               DESeq2_1.30.1              
# [3] SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [5] MatrixGenerics_1.2.1        matrixStats_1.0.0          
# [7] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
# [9] IRanges_2.24.1              S4Vectors_0.28.1           
# [11] BiocGenerics_0.36.1        

# loaded via a namespace (and not attached):
#   [1] genefilter_1.72.1      locfit_1.5-9.4         tidyselect_1.2.0      
# [4] splines_4.0.2          lattice_0.20-41        generics_0.1.0        
# [7] colorspace_2.0-2       vctrs_0.5.1            utf8_1.2.2            
# [10] blob_1.2.1             XML_3.99-0.6           survival_3.5-5        
# [13] rlang_1.0.6            pillar_1.6.2           withr_2.4.2           
# [16] glue_1.4.2             DBI_1.1.1              BiocParallel_1.24.1   
# [19] bit64_4.0.5            RColorBrewer_1.1-2     GenomeInfoDbData_1.2.4
# [22] lifecycle_1.0.3        zlibbioc_1.36.0        munsell_0.5.0         
# [25] gtable_0.3.0           memoise_2.0.0          geneplotter_1.68.0    
# [28] fastmap_1.1.0          AnnotationDbi_1.52.0   fansi_0.5.0           
# [31] Rcpp_1.0.7             xtable_1.8-4           scales_1.1.1          
# [34] cachem_1.0.5           DelayedArray_0.16.3    annotate_1.68.0       
# [37] XVector_0.30.0         bit_4.0.4              dplyr_1.0.10          
# [40] grid_4.0.2             cli_3.6.0              tools_4.0.2           
# [43] bitops_1.0-7           magrittr_2.0.1         RCurl_1.98-1.3        
# [46] RSQLite_2.2.7          tibble_3.1.3           pkgconfig_2.0.3       
# [49] crayon_1.4.1           ellipsis_0.3.2         Matrix_1.3-4          
# [52] assertthat_0.2.1       httr_1.4.2             R6_2.5.0              
# [55] compiler_4.0.2     