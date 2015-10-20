cuffnorm_loc <- 'data/cuffnorm_gencode_all'


trans <- function(x) log(x, 2)
pseudocount <- .2
samps <-  read.table("data/rna_samples2.txt", header=T, sep  = "\t", as.is = T)
rec_grade <- aggregate(Grade ~ Patient, data = samps , FUN = max)
names(rec_grade)[2] <- "Grade_rec"
samps_full <- merge(samps, rec_grade)

tp73_full <- scan(file = 'data/tp73_full.txt', what = '')
tp73_trunc <- scan(file = 'data/tp73_trun.txt', what = '')

dat <- read.table(file.path(cuffnorm_loc, "isoforms.fpkm_table"), header = T, sep = '\t', as.is=T)
tp73_full_ind <- which(dat[,1] %in% tp73_full)
tp73_trunc_ind <- which(dat[,1] %in% tp73_trunc)
dat <- dat[c(tp73_full_ind, tp73_trunc_ind),]

p_4 <- with(samps_full, RNA.seq[Tumor == "Primary" & Grade_rec == 4])
r_4 <- with(samps_full, RNA.seq[Tumor != "Primary" & Grade_rec == 4])

p_23 <- with(samps_full, RNA.seq[Tumor == "Primary" & Grade_rec != 4])
r_23 <- with(samps_full, RNA.seq[Tumor != "Primary" & Grade_rec != 4])


g4_p <- trans(dat[, grep(paste0(p_4, collapse = '|'), names(dat), value = T)]+pseudocount)
g4_r <- trans(dat[, grep(paste0(r_4, collapse = '|'), names(dat), value = T)]+pseudocount)

g23_p <- trans(dat[, grep(paste0(p_23, collapse = '|'), names(dat), value = T)]+pseudocount)
g23_r <- trans(dat[, grep(paste0(r_23, collapse = '|'), names(dat), value = T)]+pseudocount)

save(g4_p,g4_r, g23_p,g23_r, file = "data/tp73_isoforms.RData")

full_dat <- rbind(t(g23_p), t(g23_r), t(g4_p), t(g4_r))
full_small <- full_dat[,colnames(full_dat) %in% tp73_trunc_ind]
full_big <- full_dat[,colnames(full_dat) %in% tp73_full_ind]

small_mean <- rowMeans(full_small)
big_mean <- rowMeans(full_big)
gbm_groups <- rep( c("g23_p", "g23_r", "g4_p", "g4_r"), c(ncol(g23_p), ncol(g23_r), ncol(g4_p), ncol(g4_r)))

pdf('results/tp73_trunc.pdf')
boxplot(small_mean ~ gbm_groups, pch = 17, col = rep(c('green4','darkorange1'), each = 2), main = "short transcripts", ylim = c(log(.2,2), 0))
dev.off()

pdf('results/tp73_full.pdf')
boxplot(big_mean ~ gbm_groups, pch = 17, col = rep(c('green4','darkorange1'), each = 2), main = "long transcripts", ylim = c(log(.2,2), 0))
dev.off()
