# Undirected acylic graph

#####
library('DOSE')
library('DESeq2')
library('enrichplot')
library('ggnewscale')
library(org.Hs.eg.db)
library(ggplot2)
# load("/Users/sofiasalazar/clusterliigh/meta/combined/namedDGElist.RData")
load("/Users/sofiasalazar/clusterliigh/meta/combined/GO_gene_list.RData")
load("/Users/sofiasalazar/clusterliigh/meta/combined/GO_allDEG_results.RData")
outdir = "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/figures/"
#####

# In local

symbolGO_results <- setReadable(gseGO_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

cnt_enrichment <- cnetplot(symbolGO_results, showCategory = 5, node_label = 'all', foldChange = go_gene_list, colorEdge = TRUE )
cnt_enrichment$labels$size <- "No. genes associated with term"
cnt_enrichment$labels$edge_colour <- "Term association"
ggsave(paste0(outdir,"cnt_enrichmentUP.png"), width = 3000, height = 5000, units = 'px', dpi = 300, bg = "white", plot = cnt_enrichment)
dev.off()

save(cnt_enrichment, file = paste0(outdir, "cnt_enrichment.RData"))

p3 <- cnetplot(symbolGO_results, foldChange=go_gene_list, circular = TRUE, colorEdge = TRUE) 
ggsave(paste0(outdir,"cnt_circle_enrichmentUP.png"), width = 3000, height = 5000, units = 'px', dpi = 300, bg = "white", plot = p3)
dev.off()