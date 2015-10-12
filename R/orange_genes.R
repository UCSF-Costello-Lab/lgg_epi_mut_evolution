load("results/logit_limma_3.RData")
load("results/log_limma_3.RData")

p_cutoff <- .05
go_p_cutoff <- .001

load("data/annotations_with_gencode_1500_TSS_1000.RData", verbose=TRUE); anno <- anno_with_k_tss_p_named;

log_df <- subset(log_df, rownames(log_df) %in% anno$genes[anno$IlmnID %in% rownames(logit_df)])

fp_count_control <- function(p_methyl, p_gene_constant ){
  temp <- length(p_methyl)/length(p_gene_constant)*p_methyl
  temp[temp>1] <- 1
  return(  temp )
}

p_methyl_adjusted <- fp_count_control(logit_df$temp, log_df$temp)

cond_rna <- which((log_df$temp < p_cutoff) & (log_df[,'g4_fc_log_2'] >= log(2,2)) & (log_df[,'g4_fc_log_2'] - log_df[,'g2_fc_log_2']) >1 )
cond_methyl <- which((p_methyl_adjusted < p_cutoff) & (logit_df[,'g4_beta_d'] <= -.2) & (logit_df[,'g4_beta_d'] - logit_df[,'g2_beta_d']) < -.15)


dat <- log_df[cond_rna,  ]
methyl_dat <- logit_df[cond_methyl ,  ]  

cat(as.character(rownames(methyl_dat)), sep = '\n', file = 'results/sig_probes.txt')

anno_in_sig_gene <- unique(anno[which(anno$genes %in% row.names(dat)),])
sig_genes_hypeMehtyl_upExpr <- unique(anno_in_sig_gene[which(anno_in_sig_gene$IlmnID %in% row.names(methyl_dat)),"gene_name"])


cat(as.character(sig_genes_hypeMehtyl_upExpr), sep = '\n', file = "results/M3_genes.txt")

cat(as.character(row.names(logit_df)), sep = '\n', file = "results/probe_bg.txt")
cat(as.character(row.names(log_df)), sep = '\n', file = 'results/rna_bg.txt')
