# WGCNA - Choose soft threshold parameter (Step 1) Author: Sofia Salazar Modified by: Evelia Coss Date: 7/Nov/2023 Input: VSD and metadata qlogin cd 
# /mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/ module load r/4.0.2

########

library(WGCNA)
library(DESeq2)
library(tidyverse)

########
# Make a VSD count matrix with only SLE samples and desired transcript annotations
indir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/"
outdir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/figures/"
#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================
# Read the gene counts table
load(paste0(indir, "vsd2.RData")) # normalized counts vsd2 
all_data <- read.csv(paste0(indir, "all_data.csv")) # meta-data
#all_data <- read.csv("all_data.csv")
datExpr <- as.data.frame(t(as.matrix(assay(vsd2)))) 
dim(datExpr) # 318 49465 
gsg = goodSamplesGenes(datExpr, verbose = 3) 
gsg$allOK # should be TRUE
# in not true, remove offending genes and samples
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
###### 
## trait data
# create a dataframe analogous to expression data that will hold the clinical traits
samples_ID <- rownames(datExpr)
traitRows <- match(samples_ID, all_data$samples)
datTraits <- all_data[traitRows,]
rownames(datTraits) <- all_data[traitRows, 1]
# Plot
sampleTree <- hclust(dist(datExpr), method = "average") 
png(file="/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/figures/named-sle-wgcna-samples_7Nov2023.png") 
par(cex = 0.6); 
par(mar = c(0,4,2,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) 
abline(h = 140, col="red") 
dev.off()
######
## removing sample outliers
clust <- cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust) 
keepSamples = (clust==1)
datExpr = datExpr[keepSamples,] 
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================
# 2. Network construction and module detection
enableWGCNAThreads(30) # allow multi-threading
# Network construction using a soft threshold Choose a set of soft-thresholding powers (beta)
powers = c(c(1:10), seq(from = 10, to=30, by=2))
# unsigned -> nodes with positive & negative correlation are treated equally signed -> nodes with negative correlation are considered *unconnected*, 
# treated as zero Call the network topology analysis function correccion bicor
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed", corFnc = bicor)
# Plot topology results
png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/figures/samples_topology_7Nov2023.png")
# png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/named-sle-network-topology.png")
par(mfrow = c(1,2)); 
cex1 = 0.9;
# ---------- Select Beta ------------------ Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")); text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") # the closer to 1, means the network will exhibit a scale free topology
# Mean connectivity as a function of the soft-thresholding power connectivity = how correlated a gene is with all other network genes
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red") 
dev.off() 
rm(sft) 

beta=12 
Degreedata <- softConnectivity(datExpr,power=beta, minNSamples = 5 )-1 
par(mfrow=c(1,1)) pdf(file = 
"/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/figures/Degreedata_Beta.pdf",
    width = 12 , height = 10) scaleFreePlot(Degreedata, AF1=paste("DOL3, power=",beta), truncated=T) 
dev.off()
##
save(datExpr, all_data, datTraits, Degreedata, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/samplesDatExpr_7Nov2023.RData")
