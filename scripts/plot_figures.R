# figures with cowplot
library(cowplot)
library(ComplexHeatmap)
library(ggplot2)
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"


# heatmap: h2
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/heatmap_all_object.RData") 
heatmap <- grid.grabExpr(draw(h2))
# volcano plot: volcanoplot_names
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/volcanoplot_object.RData")
# pca: pca_plot
load("/mnt/Citosina/amedina/ssalazar/meta/combined/figures/pca_object.RData")

top_row <-plot_grid(pca_plot, volcanoplot_names, labels= c('A', 'B'), align = 'h', axis = 'b', rel_widths = c(1.3, 1.5)) 
bottom_row <- plot_grid(heatmap, labels = 'C')

p1 <- plot_grid(top_row, bottom_row, ncol = 1, align = 'v', axis = 'l',rel_heights = c(1,1))

ggsave(dpi  = 300 , paste0(outdir, "figure1.png"), plot = p1, width = 4000, height = 5000, units = 'px')

