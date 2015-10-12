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

full_genes <- aggregate(cbind(start , -end) ~ seqname + strand + genes + transcripts, data = gtf, FUN= min)
names(full_genes)[6] <- "end"
full_genes[,6] <- -full_genes[,6]

#### ADD PROMOTER ####
full_genes$start[full_genes$strand=='+'] <- with(full_genes, sapply(start[strand=='+'] - 1500, function(x) max(c(x, 0)))) 
full_genes$end[full_genes$strand=='-'] <- with(full_genes, end[strand=='-'] + 1500) 

### ORDER ####
full_genes <- full_genes[with(full_genes,order(seqname, start)),]
full_genes_GR <- with(full_genes, GRanges(seqnames = Rle(seqname), ranges = IRanges(start, end-1), genes = genes, transcripts = transcripts)  )


anno_in_genes <- as.data.frame(findOverlaps(anno_good_GR, full_genes_GR))

anno_genes_link <- cbind.data.frame(IlmnID = anno_good[anno_in_genes[,1],"IlmnID"], genes= full_genes[anno_in_genes[,2], "genes"])

anno_genes_link_uniq <- unique(anno_genes_link)
anno_with_genes <- merge(anno_full, anno_genes_link_uniq, all.x = T)

gtf_anno <- read.table('data/hg19_GencodeCompV19_ann.txt', sep = '\t', as.is = T)
gene_names <- unique(gtf_anno[,c('V3', 'V4')])

anno_with_genes_named <- merge(anno_with_genes, gene_names, all.x = T, by.x = 'genes', by.y = 'V3')

names(anno_with_genes_named)[ncol(anno_with_genes_named)] <- 'gene_name'

save(anno_with_genes_named, file = paste0("data/annotations_with_gencode_full_gene.RData"))  ### TSS-1.5kb & introns



