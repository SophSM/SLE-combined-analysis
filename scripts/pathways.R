

#####
library(limma)
library("AnnotationDbi")
library(org.Hs.eg.db)
library(tidyverse)
#####

# get genes and pathways 
tab <- getGeneKEGGLinks(species='hsa')
tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,
                     column="SYMBOL", keytype="ENTREZID")

ids <- str_remove_all(tab$PathwayID, "path:")
tab$PathwayID <- ids
head(tab)

# get pathway names
pth_names <- getKEGGPathwayNames(species="hsa")

# GENE - PATHWAY relation data frame
gene_path <- inner_join(tab,pth_names,by="PathwayID")
head(gene_path)

pathways <- levels(as.factor(gene_path$Description))
write.csv(as.data.frame(pathways), "/mnt/Citosina/amedina/ssalazar/meta/combined/pathways.csv")

# load heatmap information
load('/mnt/Citosina/amedina/ssalazar/meta/combined/annGSEA_allDEG_GO.RData')
genes_interest <- rownames(annGSEA)

gene_path <- gene_path %>% dplyr::select(Symbol, Description)
  
  
gene_path_interest <- gene_path[gene_path$Symbol %in% genes_interest,]

gene_path_interest %>% dplyr::count(Symbol)
interest_paths <- levels(as.factor(gene_path_interest$Description))

sle_subset <- gene_path_interest %>% subset(Description == 'Systemic lupus erythematosus - Homo sapiens (human)')
