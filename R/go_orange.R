library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(pathview)
library(clusterProfiler)

#cols <- c("gainsboro", rgb(250, 164, 26, max = 255), rgb(81, 184, 72, max = 255), "grey25", "mediumorchid4", rgb(245, 126, 32, max=255), rgb(10, 140, 68, max = 255))

probe_bg <- scan('results/probe_bg.txt', what = '')
load("data/annotations_with_gencode_1500_TSS_1000.RData", verbose=TRUE); anno <- anno_with_k_tss_p_named;
anno <- subset(anno[,c("IlmnID", "gene_name")], IlmnID %in% probe_bg)

bg <- setdiff(anno$gene_name, NA)

m3_genes <- scan('results/M3_genes.txt', what = '')
gene_entrez_link <- read.table("data/gene2entrez.txt", T, '\t', as.is = T)
e_bg <- as.character(gene_entrez_link[,2][which(gene_entrez_link[,1] %in% bg)])

e_enrich <- gene_entrez_link[,2][which(gene_entrez_link[,1] %in% m3_genes)]
go <- enrichGO(gene = as.character(e_enrich),  readable = T, pvalueCutoff = .2, ont = "BP", universe = e_bg)
if(!is.null(go)){
  write.table(summary(go), file = 'results/m3_hypo_upreg_genes_go_table.txt' , sep = '\t', row.names = F)
}
