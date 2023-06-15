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

list.genes <- c("BANK1", "BLK", "C2","IRF5", "ITGAM", "STAT4", "TNFAIP3", "TNFSF4")

# in local
load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-BANK1.RData")
bank1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-BLK.RData")
blk <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-C2.RData")
c2 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IRF5.RData")
irf5 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-ITGAM.RData")
itgam <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-STAT4.RData")
stat4 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-TNFAIP3.RData")
tnfaip3 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-TNFSF4.RData")
tnfsf4 <- p3


fig2_A <- plot_grid(bank1, blk, c2, irf5, itgam, stat4, tnfaip3, tnfsf4, labels = NULL, ncol = 2)

####

# boxplots

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-BANK1.RData")
bank1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-BLK.RData")
blk <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-C2.RData")
c2 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-IRF5.RData")
irf5 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-ITGAM.RData")
itgam <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-STAT4.RData")
stat4 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-TNFAIP3.RData")
tnfaip3 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-TNFSF4.RData")
tnfsf4 <- p3


fig2_A <- plot_grid(bank1, blk, c2, irf5, itgam, stat4, tnfaip3, tnfsf4, labels = NULL, ncol = 4)
ggsave(dpi  = 300 , paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/", "boxplots1.png"), plot = fig2_A, width = 4000, height = 2000, units = 'px')

##

list.genes_new <- c("IFI27", "OTOF", "USP18", "IFI44L", "IFI44", "RSAD2", "ISG15", "IFIT1")

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFI27.RData")
ifi27 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-OTOF.RData")
otof <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-USP18.RData")
usp18 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFI44L.RData")
ifi44l <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-RSAD2.RData")
rsad2 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-ISG15.RData")
isg15 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFIT1.RData")
ifit1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/violinplot-IFI44.RData")
ifi44 <- p3

fig2_B <- plot_grid(ifi27, otof, usp18, ifi44l, ifi44, rsad2,isg15,ifit1, labels = NULL, ncol = 2)



fig2 <- plot_grid(fig2_A, NULL, fig2_B, labels = c("A","","B"), ncol = 3, nrow = 1,
                  rel_widths = c(1, 0.3, 1))

ggsave(dpi  = 300 , paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/", "figure2.png"), plot = fig2, width = 5000, height = 4000, units = 'px')

# 
# boxplot

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-IFI27.RData")
ifi27 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-OTOF.RData")
otof <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-USP18.RData")
usp18 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-IFI44L.RData")
ifi44l <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-RSAD2.RData")
rsad2 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-ISG15.RData")
isg15 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-IFIT1.RData")
ifit1 <- p3

load("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/boxplot-IFI44.RData")
ifi44 <- p3

fig2_B <- plot_grid(ifi27, otof, usp18, ifi44l, ifi44, rsad2,isg15,ifit1, labels = NULL, ncol = 4)
ggsave(dpi  = 300 , paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/", "boxplots2.png"), plot = fig2_B, width = 4000, height = 2000, units = 'px')

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
bottom_row <- plot_grid(cnt, labels = 'C')

fig3 <- plot_grid(top_row, bottom_row, ncol = 1, align = 'v', axis = 'l',rel_heights = c(1.3,1))

ggsave(dpi  = 300 , paste0("/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/", "figure3.png"), plot = fig3, width = 6000, height = 6000, units = 'px')



