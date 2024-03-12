# figures with cowplot
library(cowplot)
library(ComplexHeatmap)
library(ggplot2)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"


#######

# FIG 1

# heatmap ordered: heat_ordered
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/heat_ordered.RData") 
heatmap <- grid.grabExpr(draw(ht_list))
# volcano plot: volcanoplot_names
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/volcanoplot_object.RData")
# pca: pca_plot
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/pca_object.RData")

top_row <-plot_grid(pca_plot, volcanoplot_names, labels= c('A', 'B'), align = 'h', axis = 'b', rel_widths = c(1.3, 1.5)) 
bottom_row <- plot_grid(heatmap, labels = 'C')

p1 <- plot_grid(top_row, bottom_row, ncol = 1, align = 'v', axis = 'l',rel_heights = c(1,1))


ggsave(dpi  = 300 , paste0(outdir, "figure1_ordered.png"), plot = p1, width = 4000, height = 5000, units = 'px')


######

# FIG 2
indir <- '/mnt/Citosina/amedina/ssalazar/meta/combined/figures/'
list.genes <- c("TREX1" , "LBH" ,   "C2" ,    "C1QC"  , "C1QB"  , "IRF5" ,  "PHRF1" , "FCGR2A")
# in local
load(paste0(indir,"violinplot-TREX1.RData"))
trex1 <- p3

load(paste0(indir, "violinplot-LBH.RData"))
lbh <- p3

load(paste0(indir, "violinplot-C2.RData"))
c2 <- p3

load(paste0(indir, "violinplot-C1QC.RData"))
c1qc <- p3

load(paste0(indir, "violinplot-C1QB.RData"))
c1qb <- p3

load(paste0(indir, "violinplot-IRF5.RData"))
irf5 <- p3

load(paste0(indir, "violinplot-PHRF1.RData"))
phrf1 <- p3

load(paste0(indir, "violinplot-FCGR2A.RData"))
fcgr2a <- p3


fig2_A <- plot_grid(trex1, lbh, c2, c1qc, c1qb, irf5, phrf1, fcgr2a, labels = NULL, ncol = 2)


##

list.genes <- c("IFI27", "OTOF", "IFI44L","SIGLEC1","USP18", "IFI44", "IFIT1", "SPATS2L")

load(paste0(indir, "violinplot-IFI27.RData"))

ifi27 <- p3

load(paste0(indir, "violinplot-OTOF.RData"))
otof <- p3


load(paste0(indir, "violinplot-IFI44L.RData"))
ifi44l <- p3

load(paste0(indir, "violinplot-SIGLEC1.RData"))
siglec1 <- p3

load(paste0(indir, "violinplot-USP18.RData"))
usp18 <- p3

load(paste0(indir, "violinplot-IFI44.RData"))
ifi44 <- p3


load(paste0(indir, "violinplot-IFIT1.RData"))
ifit1 <- p3

load(paste0(indir, "violinplot-SPATS2L.RData"))
spats2l <- p3

fig2_B <- plot_grid(ifi27, otof, ifi44l, siglec1, usp18, ifi44,ifit1,spats2l, labels = NULL, ncol = 2)



fig2 <- plot_grid(fig2_B, NULL, fig2_A, labels = c("A","","B"), ncol = 3, nrow = 1,
                  rel_widths = c(1, 0.3, 1))

ggsave(dpi  = 300 , paste0(indir, "figure2.png"), plot = fig2, width = 5000, height = 4000, units = 'px')

######

# supplementary DEGs barplots

load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/barplot_up_nc.RData")
up <- p
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/barplot_down_nc.RData")
down <- p2

sup <- plot_grid(up, down, labels = "AUTO", ncol = 1)
ggsave(dpi  = 300 , file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/sup_DEGs.png", plot = sup, width = 5000, units = 'px')

#######

# supplementary heatmap top DEGs
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/cluster_heatmap.RData") # h1
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/cluster_heatmap_rows.RData") # h2
heatmap1 <- grid.grabExpr(draw(h1_list))
heatmap2 <- grid.grabExpr(draw(h2_list))

sup2 <- plot_grid(heatmap1, heatmap2, labels ="AUTO", ncol = 2)

ggsave(dpi  = 300 , file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/sup_heatmap.png", plot = sup2, width = 5000, units = 'px')

########

# FIG 3

# in local

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/barplotGO_DOWN.RData") # g.down
load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/barplotGO_UP.RData") # g.up
load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/cnt_enrichment.RData") # cnt_enrichment

cnt <- cnt_enrichment + theme(plot.background = element_rect(color = 'white'))
top_row <-plot_grid(g.down, g.up, labels= c('A', 'B'), align = 'h')
bottom_row <- plot_grid(cnt_enrichment, ncol = 1, labels = 'C')

ggsave(dpi  = 300 , paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/", "figure3A.png"), plot = top_row, width = 6000, height = 3000, units = 'px')


fig3 <- plot_grid(top_row, cnt, ncol = 1, align = 'v', axis = 'l',rel_heights = c(1.3,1))

ggsave(dpi  = 300 , paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/", "figure3.png"), plot = fig3, width = 6000, height = 6000, units = 'px')

#####
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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#   [1] ggplot2_3.3.5        ComplexHeatmap_2.9.3 cowplot_1.1.1       

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-2  pillar_1.6.2        compiler_4.0.2     
# [4] iterators_1.0.13    digest_0.6.27       clue_0.3-59        
# [7] lifecycle_1.0.3     tibble_3.1.3        gtable_0.3.0       
# [10] pkgconfig_2.0.3     png_0.1-7           rlang_1.0.6        
# [13] foreach_1.5.2       DBI_1.1.1           cli_3.6.0          
# [16] parallel_4.0.2      withr_2.4.2         cluster_2.1.0      
# [19] dplyr_1.0.10        S4Vectors_0.28.1    generics_0.1.0     
# [22] vctrs_0.5.1         GlobalOptions_0.1.2 IRanges_2.24.1     
# [25] stats4_4.0.2        tidyselect_1.2.0    glue_1.4.2         
# [28] R6_2.5.0            GetoptLong_1.0.5    fansi_0.5.0        
# [31] magrittr_2.0.1      scales_1.1.1        codetools_0.2-18   
# [34] ellipsis_0.3.2      matrixStats_1.0.0   BiocGenerics_0.36.1
# [37] assertthat_0.2.1    shape_1.4.6         circlize_0.4.14    
# [40] colorspace_2.0-2    utf8_1.2.2          doParallel_1.0.17  
# [43] munsell_0.5.0       crayon_1.4.1        rjson_0.2.20     


