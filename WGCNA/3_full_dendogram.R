library(WGCNA)
indir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/"
outdir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/figures/"

load(paste0(indir, "full_blockwiseModules.RData"))
####
#1. Plot module dendogram

sizeGrWindow(12, 9)

plotColors <- labels2colors(net_unmer$colors)

unique(plotColors)
png(file = paste0(outdir, "clusterDendogram.png"))
  plotDendroAndColors(net_unmer$dendrograms[[1]], plotColors[net_unmer$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
