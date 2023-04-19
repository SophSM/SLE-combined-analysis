
# Heatmap for GO terms and genes
######
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
######

load("/mnt/Citosina/amedina/ssalazar/meta/combined/namedDGElist.RData")
load("/mnt/Citosina/amedina/ssalazar/meta/combined/LRT-dds.RData")
outdir = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/"

data <- df_names %>% 
  mutate(Expression = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))


keytypes(org.Hs.eg.db)
original_gene_list = data$log2FoldChange
names(original_gene_list) <- data$ID
ids <- bitr(names(original_gene_list), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')

# GET TOP GENES

top <- 70
top_genes <- bind_rows(data %>% 
                         filter(Expression == 'Up-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top),
                       data %>% 
                         filter(Expression == 'Down-regulated') %>% 
                         arrange(padj, desc(abs(log2FoldChange))) %>% 
                         head(top)
)

DGE_genes <- data[data$Expression != 'Unchanged',]

### lista para el lab
list_lab <- DGE_genes[,c(8,10,3)]
write.csv(list_lab, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/topGenes_list.csv')



list_lab_up <- list_lab[list_lab$Expression == 'Up-regulated',]
list_lab_up <- list_lab_up[order(-list_lab_up$log2FoldChange), ]
rownames(list_lab_up) <- seq(1:nrow(list_lab_up))

list_lab_down <- list_lab[list_lab$Expression == 'Down-regulated',]
list_lab_down <- list_lab_down[order(list_lab_down$log2FoldChange), ]
rownames(list_lab_down) <- seq(1:nrow(list_lab_down))

write.csv(list_lab_up, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/topGenesUp_list.csv')
write.csv(list_lab_down, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/topGenesDown_list.csv')

######

sigGenes <- top_genes$gene_name
original_gene_list <- top_genes$log2FoldChange
names(original_gene_list) <- top_genes$ID
# names(original_gene_list) <- top_genes$gene_name
ids <- bitr(names(original_gene_list), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
top_genes <- top_genes[!duplicated(top_genes$ID),]
data2 = top_genes[top_genes$ID %in% dedup_ids$ENSEMBL,]
data2$Y = dedup_ids$ENTREZID
names(original_gene_list) <- data2$Y

# GO
go_gene_list<-na.omit(original_gene_list)
go_gene_list = sort(go_gene_list, decreasing = TRUE)
save(go_gene_list, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/GO_gene_list.RData')

gseGO_res <- gseGO(geneList = go_gene_list,
                   ont = 'ALL',
                   OrgDb = org.Hs.eg.db,
                   minGSSize = 3,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   keyType = 'ENTREZID')
save(gseGO_res, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/GO_results.RData")

########

# GO with all DE genes

sigGenes <- DGE_genes$gene_name
original_gene_list <- DGE_genes$log2FoldChange
names(original_gene_list) <- DGE_genes$ID

# names(original_gene_list) <- top_genes$gene_name
ids <- bitr(names(original_gene_list), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
DGE_genes <- DGE_genes[!duplicated(DGE_genes$ID),]
data2 = DGE_genes[DGE_genes$ID %in% dedup_ids$ENSEMBL,]
data2$Y = dedup_ids$ENTREZID
names(original_gene_list) <- data2$Y

go_gene_list<-na.omit(original_gene_list)
go_gene_list = sort(go_gene_list, decreasing = TRUE)
save(go_gene_list, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/GO_allDEG_list.RData')

gseGO_res <- gseGO(geneList = go_gene_list,
                   ont = 'ALL',
                   OrgDb = org.Hs.eg.db,
                   minGSSize = 3,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   keyType = 'ENTREZID')
save(gseGO_res, file = "/mnt/Citosina/amedina/ssalazar/meta/combined/GO_allDEG_results.RData")

#########
# load

load("/mnt/Citosina/amedina/ssalazar/meta/combined/GO_allDEG_results.RData")

########

symbolGO_results <- setReadable(gseGO_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
symbolGO_results <- symbolGO_results@result

gseSubset <- subset(symbolGO_results, p.adjust <= 0.05 )

annGSEA <- data.frame(row.names = sigGenes)
for (j in 1:length(sigGenes)) {
  gene <- sigGenes[j]
  pattern <- gene
  for (k in 1:nrow(gseSubset)){ # parse through GO terms 
    if (any(grepl(pattern, gseSubset$core_enrichment[k]))) {
      annGSEA[j,k] <- 1
    } else {
      annGSEA[j,k] <- 0
    }
  }
}
colnames(annGSEA) <- gseSubset[,3]
# remove terms with no overlapping genes
annGSEA <- annGSEA[,apply(annGSEA, 2, mean)!=0]

# remove genes with no overlapping terms
annGSEA <- annGSEA[apply(annGSEA, 1, mean)!=0,]

annGSEA[1:5,1:5]

save(annGSEA, file = '/mnt/Citosina/amedina/ssalazar/meta/combined/annGSEA_allDEG_GO.RData')
########

# match the order of rownames in top_genes with annGSEA
rownames(top_genes) <- top_genes$gene_name
topTableAligned <- top_genes[which(rownames(top_genes) %in% rownames(annGSEA)),]
topTableAligned <- topTableAligned[match(rownames(annGSEA), rownames(topTableAligned)),]
all(rownames(topTableAligned) == rownames(annGSEA)) # TRUE

### all DEG

rownames(DGE_genes) <- DGE_genes$gene_name
topTableAligned <- DGE_genes[which(rownames(DGE_genes) %in% rownames(annGSEA)),]
topTableAligned <- topTableAligned[match(rownames(annGSEA), rownames(topTableAligned)),]
all(rownames(topTableAligned) == rownames(annGSEA)) # TRUE

# colour bar for fold changes for sigGenes
dfFoldChangeGenes <- data.frame(
  topTableAligned[which(rownames(topTableAligned) %in% rownames(annGSEA)), 'log2FoldChange'])

# merge both
dfGeneAnno <- data.frame(dfFoldChangeGenes)
colnames(dfGeneAnno) <- c('Log2FC')

colours <- colorRamp2(c(min(dfGeneAnno$Log2FC),0, max(dfGeneAnno$Log2FC)), c('royalblue','white','#b50f04'))

haGenes <- rowAnnotation(
  df = dfGeneAnno,
  col = list(Log2FC =colours),
  width = unit(1,'cm'),
  annotation_name_side = 'top')

# bottom annotation with enriched terms

haTerms <- HeatmapAnnotation(
  text = anno_text(
    colnames(annGSEA),
    rot = 45,
    just = 'right',
    gp = gpar(fontsize = 12)),
  annotation_height = unit(8, 'cm'),
  annotation_name_side = 'left')

# bottom annotation with GO source for each term

term_source <- symbolGO_results %>% dplyr::select(ONTOLOGY)
 
sourceTermsha <- HeatmapAnnotation(
  Source = term_source$ONTOLOGY,
  col = list(Source = c('BP' = '#3C6997', 'CC' = '#d9c621', "MF" = '#b30039'))
)

# HEATMAP

annGSEA <- as.matrix(annGSEA)
hmapGSEA <- Heatmap(annGSEA,
                    name = 'GO enrichment',
                    # split = dfGeneAnno[,2],
                    
                    col = c('0' = 'white', '1' = 'forestgreen'),
                    
                    rect_gp = gpar(col = 'grey85'),
                    
                    cluster_rows = TRUE,
                    show_row_dend = TRUE,
                    row_title = 'Top Genes',
                    row_title_side = 'left',
                    row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_title_rot = 90,
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_names_side = 'left',
                    row_dend_width = unit(35, 'mm'),
                    
                    cluster_columns = TRUE,
                    show_column_dend = TRUE,
                    column_title = 'Enriched terms',
                    column_title_side = 'top',
                    column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                    column_title_rot = 0,
                    show_column_names = FALSE,
                    
                    show_heatmap_legend = FALSE,
                    
                    clustering_distance_columns = 'euclidean',
                    clustering_method_columns = 'ward.D2',
                    clustering_distance_rows = 'euclidean',
                    clustering_method_rows = 'ward.D2',
                    
                    bottom_annotation = c(sourceTermsha, haTerms)
)

png(filename = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/GOenriched_allDEG_HMap_wTerm.png", width = 5000, height = 5000, res = 300)
draw(hmapGSEA + haGenes,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()

png(filename = "/mnt/Citosina/amedina/ssalazar/meta/combined/figures/GOenriched_allDEG_HMap.png", width = 5000, height = 5000, res = 300)
draw(hmapGSEA + haGenes,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()
