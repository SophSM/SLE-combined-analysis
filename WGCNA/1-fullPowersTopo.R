library(WGCNA)
library(DESeq2)

indir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/"
outdir <- "/mnt/Citosina/amedina/ssalazar/meta/combined/coexpression/figures/"

load(paste0(indir, "vsd2.RData"))
all_data <- read.csv(paste0(indir, "all_data.csv"))
#####
# 1. Data cleanup

datExpr <- as.data.frame(t(as.matrix(assay(vsd2))))

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

# Assess sample outliers
png(file=paste0(outdir, "tree_all_samples.png"))
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
      cex.axis = 1.5, cex.main = 2)
  abline(h = 180, col="red")
dev.off()

# Remove sample outliers

clust = cutreeStatic(sampleTree, cutHeight = 180, minSize = 10)
table(clust)
# 0   1 
# 9 309 
keepSamples = (clust==1) # removing 9 samples
datExpr = datExpr[keepSamples,]
nGenes = ncol(datExpr) # 49465 genes
nSamples = nrow(datExpr) # 309 samples

# Load trait data
dim(all_data)
all_data <- all_data[,c(2,3)] # 318   2
dim(all_data)

# Create a dataframe analogous to expression data that will hold the clinical traits
samples_ID <- rownames(datExpr)
traitRows <- match(samples_ID, all_data$samples)
datTraits = all_data[traitRows,]
rownames(datTraits) = all_data[traitRows, 1]

dim(datTraits)

######
# 2. Soft power determination

enableWGCNAThreads(12)
powers = c(c(1:10), seq(from = 10, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Plot toplogy
png(file = paste0(outdir,"full_topology.png"))

  par(mfrow = c(1,2));
  cex1 = 0.9;

  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
  abline(h=0.90,col="red") # r squared cutoff
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
dev.off()

#######
# 3. Output

save(datExpr, file = paste0(indir, "coexpression/full_datExpr.RData"))
write.csv(datTraits, file = paste0(indir, "coexpression/datTraits.csv"))
