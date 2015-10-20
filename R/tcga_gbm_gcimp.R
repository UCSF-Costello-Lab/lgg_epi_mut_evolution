
gcimp_genes <- read.csv('data/gcimp_genes.csv', T, as.is =T , skip = 3)

gcimp_genes$Gene.Name[which(gcimp_genes$Gene.Name == 'KIAA0746')] <- 'SEL1L3'
gcimp_genes$Gene.Name[which(gcimp_genes$Gene.Name == 'MOSC2')] <- 'MARC2'
gcimp_genes$Gene.Name[which(gcimp_genes$Gene.Name == 'TMEM22')] <- 'SLC35G2'
gcimp_genes$Gene.Name[which(gcimp_genes$Gene.Name == 'DKFZP586H2123')] <- 'PAMR1'
gcimp_genes$Gene.Name[which(gcimp_genes$Gene.Name == 'FLJ21963')] <- 'ACSS3'
gcimp_genes$Gene.Name[which(gcimp_genes$Gene.Name == 'C7orf46')] <- 'FAM221A'


hypo_genes <- scan('results/M3_genes.txt', what = '')
mean(hypo_genes %in% gcimp_genes$Gene.Name)

load('data/annotations_with_gencode_1500_TSS_1000_2.RData')
length(intersect(gcimp_genes$Gene.Name, anno$gene_name))


sig_probes <- scan('results/sig_probes.txt', what = '')
