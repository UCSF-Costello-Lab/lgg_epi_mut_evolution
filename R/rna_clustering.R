cuffnorm_loc <- 'data/cuffnorm_gencode_all'


trans <- function(x) log(x, 2)
pseudocount <- .2

samps <-  read.table("data/rna_samples2.txt", header=T, sep  = "\t", as.is = T)
rec_grade <- aggregate(Grade ~ Patient, data = samps , FUN = max)
names(rec_grade)[2] <- "Grade_rec"
samps_full <- merge(samps, rec_grade)

samp_trun <- with(samps, paste0(Patient, substr(Tumor,1,1), Grade))
names(samp_trun) <- samps$RNA.seq

samp_trun <- samp_trun[sort(names(samp_trun))]

dat <- read.table(file.path(cuffnorm_loc, "genes.fpkm_table"), header = T, sep = '\t', as.is=T)
full_dat <- trans( dat[, sort(grep(paste0(samps$RNA.seq, collapse='|'), names(dat), value = T))] + pseudocount)
names(full_dat) <- samp_trun



sds <- apply(full_dat, 1, sd, na.rm = T)

ord_dat <-  full_dat[order(sds, decreasing = T),]

per_cutoff <- .01

sbs_dat <- ord_dat[1:(per_cutoff*length(sds)),]


pdf("results/plots/rna_dendro_01.pdf")
distance <- dist( t(sbs_dat), method="euclidean")
#distance <- as.dist( (1 - cor(sbs_dat, use = "pairwise.complete.obs"))/2 )
hc <- hclust(distance, method="ward")
plot(hc, main="", hang = -1)
dev.off()


per_cutoff <- .50

sbs_dat <- ord_dat[1:(per_cutoff*length(sds)),]


pdf("results/plots/rna_dendro_50.pdf")
distance <- dist( t(sbs_dat), method="euclidean")
#distance <- as.dist( (1 - cor(sbs_dat, use = "pairwise.complete.obs"))/2 )
hc <- hclust(distance, method="ward")
plot(hc, main="", hang = -1)
dev.off()



