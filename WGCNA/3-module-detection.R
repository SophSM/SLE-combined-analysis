## Reduced WGCNA script
library(DESeq2)
library(WGCNA)
load("/mnt/Citosina/amedina/ssalazar/meta/combined/named-WGCNA-counts-df.RData")
datExpr <- as.data.frame(t(as.matrix(reduced.vsd)))
all_data <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/all_data.csv")
## 
# REDUCED SAMPLES
load("/mnt/Citosina/amedina/ssalazar/meta/combined/less_samplesDatExpr.RData")
all_data <- all_data[,c(2,3)]
# create a dataframe analogous to expression data that will hold the clinical traits
samples_ID <- rownames(datExpr)
traitRows <- match(samples_ID, all_data$samples)
datTraits = all_data[traitRows,]
rownames(datTraits) = all_data[traitRows, 1]


####
# ALL SAMPLES
sampleTree = hclust(dist(datExpr), method = "average")
samples_ID <- rownames(datExpr)
traitRows <- match(samples_ID, all_data$samples)
datTraits = all_data[traitRows,-1]
rownames(datTraits) = all_data[traitRows, 2]
#####
enableWGCNAThreads(30)

print("Calculating adjacencies")
start.time <- Sys.time()
softPower = 22
adjacency = adjacency(datExpr, power = softPower)
save(adjacency, file = "/mnt/Citosina/amedina/lupus/named-reduced-sle-adjacencyWGCNA.RData")
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
print("Done")


print("Calculating TOM matrix...")
start.time <- Sys.time()
TOM = TOMsimilarity(adjacency)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
print("Done")
save(TOM, file = "/mnt/Citosina/amedina/lupus/named-reduced-sle-TOM-WGCNA.RData")
print("Saved")

######
######
load("/mnt/Citosina/amedina/lupus/named-reduced-sle-TOM-WGCNA.RData")

dissTOM = 1-TOM
print("Hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes...")
geneTree = hclust(as.dist(dissTOM), method = "average")
print("Done")

png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/named-reduced-gene-dendogram.png")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

print("Dendogram saved")

#### 
print("Initiating module detection...")

minModuleSize = 30 # min number of genes per module
print("Min module size:")
print(minModuleSize)

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 1, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize, cutHeight = 0.97)
print("Done module identification")
# Convert numeric lables into colors

dynamicColors = labels2colors(dynamicMods)
length(levels(as.factor(dynamicColors)))
# 27 with cutHeight 0.9, mod size 4
# 103 with cutHeight 0.95, mod size 4
# 15 with cutHeight 0.85, mod size 4

# 27 with cutHeight 0.99, mod size 30
# 9 with cutHeight 0.97, mod size 30
table(dynamicColors)

print("Plotting dendogram with module colors...")
# Plot the dendrogram and colors underneath

png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/named-sle-colored-gene-dendogram.png")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

print("Done")
#### Merging of modules whose expression profiles are similar

print("Calculate eigengenes...")
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

print("Done")
# Plot the result

png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/named-reduced-module-eigengenes.png")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()
