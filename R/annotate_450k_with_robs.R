
load('data/data.betas.anno.Rdata')
anno_full <- data.betas.anno

library(GenomicFeatures)
anno_good <- anno_full[!is.na(anno_full$MAPINFO),]
anno_good_GR <- with(anno_good, GRanges(seqnames = Rle(paste0("chr", CHR)), ranges = IRanges(MAPINFO, MAPINFO)))




gtf <- read.delim('data/hg19_GencodeCompV19_geneid.gtf', header=FALSE)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
chronly <- c(paste0("chr",1:22), "chrX", "chrY")
gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows

gtf_genes_transcripts <- matrix(unlist(strsplit(gsub("transcript_id|gene_id| {1,}", "", gtf$attributes), ";")), ncol = 2, byrow = T)
colnames(gtf_genes_transcripts) <- c("genes", "transcripts")
gtf <- cbind.data.frame(gtf, gtf_genes_transcripts)

full_genes <- aggregate(cbind(start , -end) ~ seqname + strand + genes + transcripts, data = gtf, FUN= min)
names(full_genes)[6] <- "end"
full_genes[,6] <- -full_genes[,6]

rob_anno_dir <- "data/rob_anno/ChromHMM_States/"

anno_with_robs_anno <- anno_full
for(i in list.files(rob_anno_dir)){
  rob_anno_table <- read.table(paste0(rob_anno_dir, i), sep='\t', as.is = T)
  cat(i, '\n')
  rob_anno_sites_GR <- with(rob_anno_table, GRanges(seqnames = Rle(V1), ranges = IRanges(V2, V3-1)) )
  anno_in_rob_sites <- as.data.frame(findOverlaps(anno_good_GR, rob_anno_sites_GR))
  
  anno_in_robs_link <- cbind.data.frame(IlmnID = anno_good[anno_in_rob_sites[,1],"IlmnID"], sites= rob_anno_table[anno_in_rob_sites[,2], "V4"])

  anno_in_robs_link_uniq <- unique(anno_in_robs_link)
  anno_with_robs_anno <- merge(anno_with_robs_anno, anno_in_robs_link_uniq, all.x = T, by = "IlmnID")
  names(anno_with_robs_anno)[ncol(anno_with_robs_anno)] <- gsub(".bed$", "", i)
}

anno_with_robs_chip <- anno_with_robs_anno

save(anno_with_robs_chip, file = "data/annotations_with_robs_chip.RData")

anno_good_split <- split(anno_good, anno_good$CHR)
names(anno_good_split) <- paste0("chr",names(anno_good_split) )

full_genes_split <- split(full_genes, full_genes$seqname)

#all_chroms_list <- list()
anno_genes_link <- data.frame(stringsAsFactors = F)
for(i in names(anno_good_split)){
  genes_i <- full_genes_split[[i]]
  if(nrow(genes_i) == 0){ next }
  anno_i <- anno_good_split[[i]]
  closest_list <- list()
  distance_vec <- vector("integer", length = nrow(anno_i))
  for(k in seq_len(nrow(anno_i))){
    dist1 <- abs(anno_i$MAPINFO[k] - genes_i[,5])
    m_d1 <- min(dist1)
    dist2 <- abs(anno_i$MAPINFO[k] - genes_i[,6])
    m_d2 <- min(dist2)
#  min_vec <- vector(mode = "integer")
    if( m_d1 < m_d2 ){
  	  min_vec <- which(dist1 == m_d1)
  	  distance_vec[k] <- m_d1
    }else if( m_d1 > m_d2 ){
      min_vec <- which(dist2 == m_d2)
      distance_vec[k] <- m_d2
    }else{ 
      min_vec <-  unique(c(which(dist1 == m_d1), which(dist2 == m_d2)))
      distance_vec[k] <- m_d1
    }  	
    closest_list[[k]] <-  min_vec
  }
#  all_chroms_list[[i]] <- closest_list
  anno_ind <- rep(seq_along(closest_list), sapply(closest_list, length))
  dist_vec <- rep(distance_vec, sapply(closest_list, length))

   temp <- cbind.data.frame(IlmnID = anno_i[anno_ind,"IlmnID"], genes= genes_i[unlist(closest_list), "genes"], dist = dist_vec)
   
   anno_genes_link <- rbind.data.frame(anno_genes_link, temp)
}

gtf_anno <- read.table('data/hg19_GencodeCompV19_ann.txt', sep = '\t', as.is = T)
gene_names <- unique(gtf_anno[,c('V3', 'V4')])

gene_symbols <- gene_names[,2]
names(gene_symbols) <- gene_names[,1]

anno_genes_link_uniq <- unique(anno_genes_link)
anno_genes_link_uniq[,1] <- as.character(anno_genes_link_uniq[,1])
anno_genes_link_uniq[,2] <- as.character(anno_genes_link_uniq[,2])
anno_genes_link_uniq[,"gene_name"] <- gene_symbols[anno_genes_link_uniq[,2]]

anno_with_nearest_genes <- merge(anno_full, anno_genes_link_uniq, all.x = T, by = "IlmnID")
save(anno_with_nearest_genes, file = "data/annotations_with_nearest_genes.RData")
