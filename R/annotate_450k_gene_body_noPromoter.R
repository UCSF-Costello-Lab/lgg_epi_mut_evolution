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


gtf_ex <- subset(gtf, feature == "exon")
gtf_ex_GR <- with(gtf_ex, GRanges(seqnames = Rle(seqname), ranges = IRanges(start, end-1), genes = genes, transcripts = transcripts)  )


gtf_ex_split_GR <- split(gtf_ex_GR, gtf_ex_GR$transcripts)
gtf_ex_split_GR_ordered <- gtf_ex_split_GR[all_first_ex_GR$transcripts]

#ex_noProm_split_GR <- mapply(function(exon_GR, prom_GR){
#  setdiff(exon_GR, prom_GR)
#}, gtf_ex_split_GR_ordered, all_first_ex_GR, SIMPLIFY=F)

library(parallel)
ex_noProm_split_GR <- mclapply(seq_along(all_first_ex_GR), function(ind) setdiff(gtf_ex_split_GR_ordered[[ind]], all_first_ex_GR[ind]), mc.cores = 24)

ex_noProm_split_GRL <- GRangesList(ex_noProm_split_GR)


anno_in_gene_noProm <- as.data.frame(findOverlaps(anno_good_GR, ex_noProm_split_GRL))

anno_in_gene_link <- cbind.data.frame(IlmnID = anno_good[anno_in_gene_noProm[,1],"IlmnID"], genes= names(gtf_ex_split_GR_ordered)[anno_in_gene_noProm[,2]])

anno_in_gene_link_uniq <- unique(anno_in_gene_link)
anno_with_inGene <- merge(anno_full, anno_in_gene_link_uniq, all.x = T)

gtf_anno <- read.table('data/hg19_GencodeCompV19_ann.txt', sep = '\t', as.is = T)
gene_names <- unique(gtf_anno[,c('V1', 'V4')])

anno_with_inGene_named <- merge(anno_with_inGene, gene_names, all.x = T, by.x = 'genes', by.y = 'V1')

names(anno_with_inGene_named)[ncol(anno_with_inGene_named)] <- 'gene_name'

save(anno_with_inGene_named, file = 'data/annotations_with_genecode_gene_body_noPromoter.RData')
