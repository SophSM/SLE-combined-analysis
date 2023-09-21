# WGCNA reduced
######
library(DESeq2)
library(WGCNA)
######

# Get soft power threshold for module detection, clean samples.

#load("/mnt/Citosina/amedina/ssalazar/meta/combined/named-WGCNA-counts-df.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/named-sle.vsd.RData")
reduced.vsd <- sle.vsd
dim(reduced.vsd)
datExpr <- as.data.frame(t(as.matrix(reduced.vsd)))
dim(datExpr) # 318 33887
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK # should be TRUE

# in not true, remove offending genes and samples
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr), method = "average")

png(file="/mnt/Citosina/amedina/ssalazar/meta/combined/figures/named-sle-wgcna-samples.png")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 140, col="red")
dev.off()
######
## removing sample outliers

clust = cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

###### 
# Load trait data
all_data <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/all_data.csv")
dim(all_data)

all_data <- all_data[all_data$DISEASE == "SLE", ]
all_data <- all_data[,c(2,3)]
# create a dataframe analogous to expression data that will hold the clinical traits
samples_ID <- rownames(datExpr)
traitRows <- match(samples_ID, all_data$samples)
datTraits = all_data[traitRows,]
rownames(datTraits) = all_data[traitRows, 1]

########
#####
# 2. Network construction and module detection
enableWGCNAThreads(15) # allow multi-threading


# Network construction using a soft threshold

# Choose a set of soft-thresholding powers (beta)
powers = c(c(1:10), seq(from = 10, to=30, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot topology results
png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/less_samples_topology.png")

# png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/named-sle-network-topology.png")

par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") # the closer to 1, means the network will exhibit a scale free topology


# Mean connectivity as a function of the soft-thresholding power
# connectivity = how correlated a gene is with all other network genes
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

##
save(datExpr, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/less_samplesDatExpr.RData")
