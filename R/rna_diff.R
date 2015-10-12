cuffnorm_loc <- 'data/cuffnorm_gencode_all'


trans <- function(x) log(x, 2)
pseudocount <- .2
samps <-  read.table("data/rna_samples2.txt", header=T, sep  = "\t", as.is = T)
rec_grade <- aggregate(Grade ~ Patient, data = samps , FUN = max)
names(rec_grade)[2] <- "Grade_rec"
samps_full <- merge(samps, rec_grade)


dat <- read.table(file.path(cuffnorm_loc, "genes.fpkm_table"), header = T, sep = '\t', as.is=T)
p_4 <- with(samps_full, RNA.seq[Tumor == "Primary" & Grade_rec == 4])
r_4 <- with(samps_full, RNA.seq[Tumor != "Primary" & Grade_rec == 4])

p_3 <- with(samps_full, RNA.seq[Tumor == "Primary" & Grade_rec == 3])
r_3 <- with(samps_full, RNA.seq[Tumor != "Primary" & Grade_rec == 3])

p_2 <- with(samps_full, RNA.seq[Tumor == "Primary" & Grade_rec == 2])
r_2 <- with(samps_full, RNA.seq[Tumor != "Primary" & Grade_rec == 2])

p_23 <- with(samps_full, RNA.seq[Tumor == "Primary" & Grade_rec != 4])
r_23 <- with(samps_full, RNA.seq[Tumor != "Primary" & Grade_rec != 4])


g4_p <- trans(dat[, grep(paste0(p_4, collapse = '|'), names(dat), value = T)]+pseudocount)
g4_r <- trans(dat[, grep(paste0(r_4, collapse = '|'), names(dat), value = T)]+pseudocount)

g3_p <- trans(dat[, grep(paste0(p_3, collapse = '|'), names(dat), value = T)]+pseudocount)
g3_r <- trans(dat[, grep(paste0(r_3, collapse = '|'), names(dat), value = T)]+pseudocount)

g2_p <- trans(dat[, grep(paste0(p_2, collapse = '|'), names(dat), value = T)]+pseudocount)
g2_r <- trans(dat[, grep(paste0(r_2, collapse = '|'), names(dat), value = T)]+pseudocount)

g23_p <- trans(dat[, grep(paste0(p_23, collapse = '|'), names(dat), value = T)]+pseudocount)
g23_r <- trans(dat[, grep(paste0(r_23, collapse = '|'), names(dat), value = T)]+pseudocount)



g4_d <- g4_r-g4_p
g3_d <- g3_r-g3_p
g2_d <- g2_r-g2_p
g23_d <- g23_r-g23_p



library(limma)
design <- model.matrix(~factor(c(rep(1,ncol(g4_d)), rep(2,ncol(g23_d))))) 
colnames(design) <- c("intercept", "slope") 
eset <- cbind.data.frame(g4_d, g23_d)
fit.g4.g23 <- lmFit(eset, design) 
fit2.g4.g23 <- ebayes(fit.g4.g23)
fit2.g4.g23.pvalues <- fit2.g4.g23$p[,2]
temp <- fit2.g4.g23.pvalues
temp[is.na(temp)] <- 1


fit.all <- lmFit(eset) 
fit2.all <- ebayes(fit.all)
temp_all <- fit2.all$p
temp_all[is.na(temp_all)] <- 1


fit.g4 <- lmFit(g4_d) 
fit2.g4 <- ebayes(fit.g4)
temp4 <- fit2.g4$p
temp4[is.na(temp4)] <- 1


fit.g2 <- lmFit(g2_d) 
fit2.g2 <- ebayes(fit.g2)
temp2 <- fit2.g2$p
temp2[is.na(temp2)] <- 1

fit.g3 <- lmFit(g3_d) 
fit2.g3 <- ebayes(fit.g3)
temp3 <- fit2.g3$p
temp3[is.na(temp3)] <- 1

fit.g23 <- lmFit(g23_d) 
fit2.g23 <- ebayes(fit.g23)
temp23 <- fit2.g23$p
temp23[is.na(temp23)] <- 1



log_df_1 <- cbind.data.frame(all_fc_log_2 = rowMeans(cbind.data.frame(g4_d, g23_d), na.rm=T), fit2.all$t, temp_all)
row.names(log_df_1) <- dat[,1]
save(log_df_1, file = "results/log_limma_1.RData")


log_df_2d <- cbind.data.frame(g4_fc_log_2 = rowMeans(g4_d, na.rm=T), fit2.g4$t, temp4)
row.names(log_df_2d) <- dat[,1]
save(log_df_2d, file = "results/log_limma_2d.RData")


log_df_2a <- cbind.data.frame(g2_fc_log_2 = rowMeans(g2_d, na.rm=T), fit2.g2$t, temp2)
row.names(log_df_2a) <- dat[,1]
save(log_df_2a, file = "results/log_limma_2a.RData")

log_df_2b <- cbind.data.frame(g3_fc_log_2 = rowMeans(g3_d, na.rm=T), fit2.g3$t, temp3)
row.names(log_df_2b) <- dat[,1]
save(log_df_2b, file = "results/log_limma_2b.RData")

log_df_2c <- cbind.data.frame(g23_fc_log_2 = rowMeans(g23_d, na.rm=T), fit2.g23$t, temp23)
row.names(log_df_2c) <- dat[,1]
save(log_df_2c, file = "results/log_limma_2c.RData")


log_df <- cbind.data.frame(g4_fc_log_2 = rowMeans(g4_d, na.rm=T), g2_fc_log_2 = rowMeans(g2_d, na.rm=T),fit2.g4.g23$t, temp)
row.names(log_df) <- dat[,1]
save(log_df, file = "results/log_limma_3.RData")


