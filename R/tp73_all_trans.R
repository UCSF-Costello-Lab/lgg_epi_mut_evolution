setwd('data')

myfile <- "hg19_GencodeCompV19_geneid.gtf"
#library(stringr)
gtf <- read.delim(myfile, header=FALSE, as.is = T)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
gtf_genes_transcripts <- matrix(unlist(strsplit(gsub("transcript_id|gene_id| {1,}", "", gtf$attributes), ";")), ncol = 2, byrow = T)
colnames(gtf_genes_transcripts) <- c("genes", "transcripts")
gtf <- cbind.data.frame(gtf, gtf_genes_transcripts)


anns <- read.table("hg19_GencodeCompV19_ann.txt", F, as.is = T)

trans <- anns[which(anns[,2]=="TP73"),]

gtf_tp73 <- gtf[which(gtf$transcripts %in% trans[,1]),]

trans_starts <- aggregate(start ~ transcripts, data = gtf_tp73, FUN = min)

trans_anno <- merge(trans_starts, trans, by = 1)

cat(as.character(trans_anno[which(trans_anno[,2] == 3607236),1]), sep='\n', file = 'tp73_trun.txt')


cat(as.character(trans_anno[which(trans_anno[,2] %in% c(3569084, 3569129, 3598930) ),1]), sep='\n', file = 'tp73_full.txt')
