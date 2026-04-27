library(WGCNA)

indir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/"
outdir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/"

load(paste0(indir, "full_datExpr.RData"))

#####
# 1. Blockwise module detection

enableWGCNAThreads(30)

start.time <- Sys.time()
net_unmer <- blockwiseModules(datExpr, maxBlockSize = 8000, corType = "pearson",
                 power = 20, TOMType = "signed", saveTOMs = FALSE,
                 minModuleSize = 30, nThreads = 30, verbose = 5)

end.time <- Sys.time()
time.taken <- end.time - start.time
print("Done, time taken: ")
print(time.taken)
print("Saving...")
save(net_unmer, file = paste0(outdir, "full_blockwiseModules.RData"))
print("DONE!")