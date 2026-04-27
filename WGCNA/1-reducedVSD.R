# WGCNA reduced

########
library(DESeq2)
library(tidyverse)
########

# Make a VSD count matrix with only SLE samples and desired transcript annotations

load("/mnt/Citosina/amedina/ssalazar/meta/combined/vsd2.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData")
DGE <- res2

transcripts <- rownames(vsd2) 
new <- list()
for(i in 1:length(transcripts)){
  new[i] <- gsub("\\..*","", transcripts[i])
}
transcripts <- (unlist(new))

df <- as.data.frame(DGE)
df <- tibble::rownames_to_column(df, "ID")

name_short<-function(table){
  table$ID <- gsub("\\..*","", table$ID)
  return(table)
}

df <- name_short(df)

library("biomaRt")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","transcript_biotype"),
               filters = c("ensembl_gene_id"), 
               values = df$ID,
               mart = mart)
names(genes)[1] <- "ID"
names(genes)[2] <- "gene_name"

keep <- c("lncRNA","miRNA","protein_coding","protein_coding_CDS_not_defined","snoRNA","snRNA")
new.genes <- subset(genes, transcript_biotype %in% keep)

dim(new.genes[unique(new.genes$ID),])
info <- new.genes[!duplicated(new.genes$ID),]
dim(info)
info <- (info[info$gene_name != "",])

save(info, file= "/mnt/Citosina/amedina/ssalazar/meta/combined/named-transcripts-info.RData")
keep.rownames <- info$ID
length(keep.rownames)

counts <- as.data.frame(assays(vsd2))
counts <- counts[,-c(1,2)]

rownames(counts) <- transcripts


reduced.vsd <- counts[rownames(counts) %in% keep.rownames,]
dim(reduced.vsd)

r <- rownames(reduced.vsd)

all(r %in% keep.rownames) # TRUE

save(reduced.vsd, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/named-WGCNA-counts-df.RData")

####
# vsd with only SLE samples

load("/mnt/Citosina/amedina/ssalazar/meta/combined/WGCNA-counts-df.RData")
all_data <- read.csv("/mnt/Citosina/amedina/ssalazar/meta/combined/all_data.csv")
sle <- all_data[all_data$DISEASE == "SLE",]$samples

sle.vsd <- reduced.vsd[,colnames(reduced.vsd) %in% sle]
dim(sle.vsd)

save(sle.vsd, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/named-sle.vsd.RData")

