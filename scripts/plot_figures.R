# figures with cowplot
library(cowplot)
library(ComplexHeatmap)
library(ggplot2)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"


#######

# FIG 1

# heatmap: h2
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/heatmap_all_object.RData") 
heatmap <- grid.grabExpr(draw(h2))
# heatmap ordered: heat_ordered
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/heat_ordered.RData") 
heatmap <- grid.grabExpr(draw(heat_ordered))
# volcano plot: volcanoplot_names
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/volcanoplot_object.RData")
# pca: pca_plot
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/pca_object.RData")

top_row <-plot_grid(pca_plot, volcanoplot_names, labels= c('A', 'B'), align = 'h', axis = 'b', rel_widths = c(1.3, 1.5)) 
bottom_row <- plot_grid(heatmap, labels = 'C')

p1 <- plot_grid(top_row, bottom_row, ncol = 1, align = 'v', axis = 'l',rel_heights = c(1,1))


ggsave(dpi  = 300 , paste0(outdir, "figure1.png"), plot = p1, width = 4000, height = 5000, units = 'px')
ggsave(dpi  = 300 , paste0(outdir, "figure1_ordered.png"), plot = p1, width = 4000, height = 5000, units = 'px')


######

# FIG 2

list.genes <- c("TREX1" , "LBH" ,   "C2" ,    "C1QC"  , "C1QB"  , "IRF5" ,  "PHRF1" , "FCGR2A")
# in local
load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-TREX1.RData")
trex1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-LBH.RData")
lbh <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-C2.RData")
c2 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-C1QC.RData")
c1qc <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-C1QB.RData")
c1qb <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IRF5.RData")
irf5 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-PHRF1.RData")
phrf1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-FCGR2A.RData")
fcgr2a <- p3


fig2_A <- plot_grid(trex1, lbh, c2, c1qc, c1qb, irf5, phrf1, fcgr2a, labels = NULL, ncol = 2)


##

list.genes <- c("IFI27", "OTOF", "IFI44L","SIGLEC1","USP18", "IFI44", "IFIT1", "SPATS2L")

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFI27.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-IFI27.RData")

ifi27 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-OTOF.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-OTOF.RData")

otof <- p3


load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFI44L.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-IFI44L.RData")

ifi44l <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-SIGLEC1.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-SIGLEC1.RData")

siglec1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-USP18.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-USP18.RData")
usp18 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFI44.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-IFI44.RData")

ifi44 <- p3


load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFIT1.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-IFIT1.RData")

ifit1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-SPATS2L.RData")
# load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/violinplot-SPATS2L.RData")

spats2l <- p3

fig2_B <- plot_grid(ifi27, otof, ifi44l, siglec1, usp18, ifi44,ifit1,spats2l, labels = NULL, ncol = 2)



fig2 <- plot_grid(fig2_B, NULL, fig2_A, labels = c("A","","B"), ncol = 3, nrow = 1,
                  rel_widths = c(1, 0.3, 1))

ggsave(dpi  = 300 , paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/", "figure2.png"), plot = fig2, width = 5000, height = 4000, units = 'px')

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
heatmap1 <- grid.grabExpr(draw(h1))
heatmap2 <- grid.grabExpr(draw(h2))

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



