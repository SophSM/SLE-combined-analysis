# WGCNA

########
library(DESeq2)
library(WGCNA)
########
# 1. Data input

# loading expression data (vsd matrix)

load("/mnt/Citosina/amedina/ssalazar/meta/combined/vsd2.RData")
dim(assay(vsd2)) # 49465 transcripts = rows, 318 samples = columns

# transpose: samples = rows, transcripts = columns
datExpr <- as.data.frame(t(as.matrix(assay(vsd2))))
dim(datExpr) # 318 49465

# check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr, verbose = 3);
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

# cluster samples to check for outliers

sampleTree = hclust(dist(datExpr), method = "average");

png(file="/mnt/Citosina/amedina/ssalazar/meta/combined/figures/wgcna-samples.png")
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

####
# 1.a load trait data
all_data <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/all_data.csv")
dim(all_data)

# create a dataframe analogous to expression data that will hold the clinical traits
samples_ID <- rownames(datExpr)
traitRows <- match(samples_ID, all_data$samples)
datTraits = all_data[traitRows,-1]
rownames(datTraits) = all_data[traitRows, 2]

########
#####
# 2. Network construction and module detection
enableWGCNAThreads(15) # allow multi-threading


# Network construction using a soft threshold

# Choose a set of soft-thresholding powers (beta)
powers = c(c(1:10), seq(from = 12, to=22, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot topology results
png(file = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/network-topology.png")

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

#####
# 2.a Co-expression similarity and adjacency
softPower = 20 # calculate adjacencies, using soft thresholding power 6
adjacency = adjacency(datExpr, power = softPower) # around 15 min with 15 threads

save(adjacency, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/adjacencyWGCNA.RData")

#####
# 2.b Topological Oveflap Matrix (TOM) - 1st step of MODULE DETECTION
# To minimize effects of noise and spurious associations, we transform the 
# adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

save(TOM, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/TOM-WGCNA.RData")

