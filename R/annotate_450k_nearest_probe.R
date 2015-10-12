
load('data/data.betas.anno.Rdata')
anno_full <- data.betas.anno

bg_probes <- scan('results/probe_bg.txt', what = '')
anno_full <- subset(anno_full, IlmnID %in% bg_probes)

library(GenomicFeatures)
anno_good <- anno_full[!is.na(anno_full$MAPINFO),]
anno_good_GR <- with(anno_good, GRanges(seqnames = Rle(paste0("chr", CHR)), ranges = IRanges(MAPINFO, MAPINFO)))


anno_good_split <- split(anno_good, anno_good$CHR)
names(anno_good_split) <- paste0("chr",names(anno_good_split) )

anno_near_probe_link <- data.frame(stringsAsFactors = F)
#anno_chr <- anno_good_split[[1]]
for(anno_chr in anno_good_split){
  cat(anno_chr$CHR[1], '\n')
  ordr <- sapply(anno_chr$MAPINFO, function(x){
    probe_dist <- x - anno_chr$MAPINFO
    order(abs(probe_dist))[2]
  } )
  temp <- cbind.data.frame(anno_chr$IlmnID, nearest_probe = anno_chr$IlmnID[ordr], dist = (anno_chr$MAPINFO - anno_chr$MAPINFO[ordr]))
  
  anno_near_probe_link <- rbind.data.frame(anno_near_probe_link, temp)
}

save(anno_near_probe_link, file = "data/annotations_nearest_probe_link.RData")




