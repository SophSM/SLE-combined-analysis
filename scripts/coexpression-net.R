# WGCNA

########
library(DESeq2)
library(WGCNA)
########
# 1. Data input

# loading expression data (vsd matrix)

load("/mnt/Citosina/amedina/ssalazar/meta/combined/vsd2.RData")
dim(assay(vsd2)) # 49465 transcripts = rows, 318 samples = columns
head(names(assay(vsd2),10))

# transpose: samples = rows, transcripts = columns
datExpr <- as.data.frame(t(as.matrix(assay(vsd2))))
dim(datExpr) # 318 49465
head(datExpr, 1)
head(rownames(datExpr), 10)

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
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

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
# 2. Network construction and module detection
enableWGCNAThreads() # allow multi-threading

## 
# Network construction using a soft threshold

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

