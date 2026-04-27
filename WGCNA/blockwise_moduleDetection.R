library(WGCNA)
library(DESeq2)

enableWGCNAThreads(15)
all_data <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/all_data.csv")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/less_samplesDatExpr.RData")
outdir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/"

all_data <- all_data[,c(2,3)]
# create a dataframe analogous to expression data that will hold the clinical traits
samples_ID <- rownames(datExpr)
traitRows <- match(samples_ID, all_data$samples)
datTraits = all_data[traitRows,]
rownames(datTraits) = all_data[traitRows, 1]

start.time <- Sys.time()
net_unmer <- blockwiseModules(datExpr, power = 22,
                              TOMType = "signed", minModuleSize = 30,
                              networkType = "signed",
                              maxBlockSize = 8000,
                              mergeCutHeight = FALSE, deepSplit = 2,
                              corType = "pearson",
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = FALSE,
                              verbose = 4)
end.time <- Sys.time()
time.taken <- end.time - start.time
print("Done, time:")
print(time.taken)
print("Saving...")
save(net_unmer, file = paste0(outdir, "blockwiseModules.RData"))
