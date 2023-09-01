# Undirected acylic graph

#####
library('DOSE')
library('DESeq2')
library('enrichplot')
library('ggnewscale')
library(org.Hs.eg.db)
library(ggplot2)
# load("/Users/sofiasalazar/clusterliigh/meta/combined/namedDGElist.RData")
load("/Users/sofiasalazar/clusterliigh/meta/combined/GO_gene_list.RData")
load("/Users/sofiasalazar/clusterliigh/meta/combined/GO_allDEG_results.RData")
outdir = "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/"
#####

# In local

symbolGO_results <- setReadable(gseGO_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

cnt_enrichment <- cnetplot(symbolGO_results, showCategory = 5, node_label = 'all', foldChange = go_gene_list, colorEdge = TRUE )
cnt_enrichment$labels$size <- "No. genes associated with term"
cnt_enrichment$labels$edge_colour <- "Term association"
ggsave(paste0(outdir,"cnt_enrichmentUP.png"), width = 3000, height = 5000, units = 'px', dpi = 300, bg = "white", plot = cnt_enrichment)
dev.off()

save(cnt_enrichment, file = paste0(outdir, "cnt_enrichment.RData"))

p3 <- cnetplot(symbolGO_results, foldChange=go_gene_list, circular = TRUE, colorEdge = TRUE) 
ggsave(paste0(outdir,"cnt_circle_enrichmentUP.png"), width = 3000, height = 5000, units = 'px', dpi = 300, bg = "white", plot = p3)
dev.off()

#####
sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.2.1

# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#   [1] ggplot2_3.3.5               org.Hs.eg.db_3.15.0        
# [3] AnnotationDbi_1.58.0        ggnewscale_0.4.9           
# [5] enrichplot_1.16.2           DESeq2_1.36.0              
# [7] SummarizedExperiment_1.26.1 Biobase_2.56.0             
# [9] MatrixGenerics_1.8.1        matrixStats_1.0.0          
# [11] GenomicRanges_1.48.0        GenomeInfoDb_1.32.4        
# [13] IRanges_2.30.1              S4Vectors_0.34.0           
# [15] BiocGenerics_0.44.0         DOSE_3.22.1                

# loaded via a namespace (and not attached):
#   [1] fgsea_1.22.0           colorspace_2.1-0       ggtree_3.4.4          
# [4] qvalue_2.28.0          XVector_0.36.0         aplot_0.1.10          
# [7] rstudioapi_0.14        farver_2.1.1           graphlayouts_1.0.0    
# [10] ggrepel_0.9.3          bit64_4.0.5            scatterpie_0.2.1      
# [13] fansi_1.0.4            codetools_0.2-19       splines_4.2.2         
# [16] cachem_1.0.8           GOSemSim_2.22.0        geneplotter_1.74.0    
# [19] knitr_1.43             polyclip_1.10-4        jsonlite_1.8.5        
# [22] annotate_1.74.0        GO.db_3.15.0           png_0.1-8             
# [25] ggforce_0.4.1          compiler_4.2.2         httr_1.4.6            
# [28] Matrix_1.5-4.1         fastmap_1.1.1          lazyeval_0.2.2        
# [31] limma_3.52.4           cli_3.6.1              tweenr_2.0.2          
# [34] htmltools_0.5.5        tools_4.2.2            igraph_1.4.3          
# [37] gtable_0.3.3           glue_1.6.2             GenomeInfoDbData_1.2.8
# [40] reshape2_1.4.4         DO.db_2.9              dplyr_1.1.2           
# [43] fastmatch_1.1-3        Rcpp_1.0.11            vctrs_0.6.3           
# [46] Biostrings_2.64.1      ape_5.7-1              nlme_3.1-162          
# [49] ggraph_2.1.0           xfun_0.39              stringr_1.5.0         
# [52] lifecycle_1.0.3        XML_3.99-0.14          zlibbioc_1.42.0       
# [55] MASS_7.3-60            scales_1.2.1           tidygraph_1.2.3       
# [58] parallel_4.2.2         RColorBrewer_1.1-3     yaml_2.3.7            
# [61] memoise_2.0.1          gridExtra_2.3          ggfun_0.0.9           
# [64] yulab.utils_0.0.6      stringi_1.7.12         RSQLite_2.3.1         
# [67] genefilter_1.78.0      tidytree_0.4.2         BiocParallel_1.30.4   
# [70] rlang_1.1.1            pkgconfig_2.0.3        bitops_1.0-7          
# [73] evaluate_0.21          lattice_0.21-8         purrr_1.0.1           
# [76] treeio_1.20.2          patchwork_1.1.2        shadowtext_0.1.2      
# [79] bit_4.0.5              tidyselect_1.2.0       plyr_1.8.8            
# [82] magrittr_2.0.3         R6_2.5.1               generics_0.1.3        
# [85] DelayedArray_0.22.0    DBI_1.1.3              pillar_1.9.0          
# [88] withr_2.5.0            survival_3.5-5         KEGGREST_1.36.3       
# [91] RCurl_1.98-1.12        tibble_3.2.1           crayon_1.5.2          
# [94] utf8_1.2.3             rmarkdown_2.22         viridis_0.6.3         
# [97] locfit_1.5-9.8         grid_4.2.2             data.table_1.14.8     
# [100] blob_1.2.4             digest_0.6.32          xtable_1.8-4          
# [103] tidyr_1.3.0            gridGraphics_0.5-1     munsell_0.5.0         
# [106] viridisLite_0.4.2      ggplotify_0.1.0    