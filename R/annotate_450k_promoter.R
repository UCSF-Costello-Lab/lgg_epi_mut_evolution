
load('data/data.betas.anno.Rdata')
anno_full <- data.betas.anno

library(GenomicFeatures)
anno_good <- anno_full[!is.na(anno_full$MAPINFO),]
anno_good_GR <- with(anno_good, GRanges(seqnames = Rle(paste0("chr", CHR)), ranges = IRanges(MAPINFO, MAPINFO)))




gtf <- read.delim('data/hg19_GencodeCompV19_geneid.gtf', header=FALSE)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
#chronly <- c(paste0("chr",1:22), "chrX", "chrY")
chronly <- paste0("chr",1:22)
gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows

gtf_genes_transcripts <- matrix(unlist(strsplit(gsub("transcript_id|gene_id| {1,}", "", gtf$attributes), ";")), ncol = 2, byrow = T)
colnames(gtf_genes_transcripts) <- c("genes", "transcripts")
gtf <- cbind.data.frame(gtf, gtf_genes_transcripts)



k <- 1500
p <- 1000
gtf_ex <- subset(gtf, feature == "exon")
min_ex_forward <- aggregate(start ~ seqname + strand + genes + transcripts, data = gtf_ex, FUN= min, subset= (strand == '+'))
min_ex_reverse <- aggregate(end ~ seqname + strand + genes + transcripts, data = gtf_ex, FUN= max, subset= (strand == '-'))

first_ex_for <- merge(min_ex_forward, gtf_ex, all.x = T)
first_ex_rev <- merge(min_ex_reverse, gtf_ex, all.x = T)

first_ex_for$end <-  first_ex_for$start+p
first_ex_rev$start <- first_ex_rev$end-p
first_ex_for$start <- first_ex_for$start-k
first_ex_rev$end <- first_ex_rev$end+k

all_first_ex <- rbind.data.frame(first_ex_for, first_ex_rev)
all_first_ex_GR <- with(all_first_ex, GRanges(seqnames = Rle(seqname), ranges = IRanges(start, end-1), genes = genes, transcripts = transcripts)  )

anno_in_first_ex <- as.data.frame(findOverlaps(anno_good_GR, all_first_ex_GR))

anno_first_ex_link <- cbind.data.frame(IlmnID = anno_good[anno_in_first_ex[,1],"IlmnID"], genes= all_first_ex[anno_in_first_ex[,2], "genes"])

anno_first_ex_link_uniq <- unique(anno_first_ex_link)
anno_with_k_tss_p <- merge(anno_full, anno_first_ex_link_uniq, all.x = T)

gtf_anno <- read.table('data/hg19_GencodeCompV19_ann.txt', sep = '\t', as.is = T)
gene_names <- unique(gtf_anno[,c('V3', 'V4')])

anno_with_k_tss_p_named <- merge(anno_with_k_tss_p, gene_names, all.x = T, by.x = 'genes', by.y = 'V3')

names(anno_with_k_tss_p_named)[ncol(anno_with_k_tss_p_named)] <- 'gene_name'

save(anno_with_k_tss_p_named, file = paste0("data/annotations_with_gencode_",k ,"_TSS_",p,".RData"))




