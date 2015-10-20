
load('data/sample_inds.RData')
load('data/data.clean.Rdata')


trans <- function(x) log(x/(1-x))

g2_p <- data.clean[,pri.v1.pairs.r1g2]
g3_p <- data.clean[,pri.v1.pairs.r1g3]
g4_p <- data.clean[,pri.v1.pairs.r1g4]
g23_p <- data.clean[,pri.v1.pairs.r1g23]



g2_r <- data.clean[,rec1.v1.pairs.r1g2]
g3_r <- data.clean[,rec1.v1.pairs.r1g3]
g4_r <- data.clean[,rec1.v1.pairs.r1g4]
g23_r <- data.clean[,rec1.v1.pairs.r1g23]

g4_d <- g4_r - g4_p
g3_d <- g3_r - g3_p
g2_d <- g2_r - g2_p
g23_d <- g23_r - g23_p


t_g2_p <- trans(g2_p)
t_g3_p <- trans(g3_p)
t_g4_p <- trans(g4_p)
t_g23_p <- trans(g23_p)

t_g2_r <- trans(g2_r)
t_g3_r <- trans(g3_r)
t_g4_r <- trans(g4_r)
t_g23_r <- trans(g23_r)


t_g4_d <- t_g4_r - t_g4_p
t_g3_d <- t_g3_r - t_g3_p
t_g2_d <- t_g2_r - t_g2_p
t_g23_d <- t_g23_r - t_g23_p


library(limma)
design <- model.matrix(~factor(c(rep(1,ncol(t_g4_d)), rep(2,ncol(t_g23_d))))) 
colnames(design) <- c("intercept", "slope") 
eset <- cbind.data.frame(t_g4_d, t_g23_d)
fit.g4.g23 <- lmFit(eset, design) 
fit2.g4.g23 <- ebayes(fit.g4.g23)
fit2.g4.g23.pvalues <- fit2.g4.g23$p[,2]
temp <- fit2.g4.g23.pvalues
temp[is.na(temp)] <- 1


fit.all <- lmFit(eset) 
fit2.all <- ebayes(fit.all)
temp_all <- fit2.all$p
temp_all[is.na(temp_all)] <- 1


fit.g4 <- lmFit(t_g4_d) 
fit2.g4 <- ebayes(fit.g4)
temp4 <- fit2.g4$p
temp4[is.na(temp4)] <- 1


fit.g2 <- lmFit(t_g2_d) 
fit2.g2 <- ebayes(fit.g2)
temp2 <- fit2.g2$p
temp2[is.na(temp2)] <- 1

fit.g3 <- lmFit(t_g3_d) 
fit2.g3 <- ebayes(fit.g3)
temp3 <- fit2.g3$p
temp3[is.na(temp3)] <- 1

fit.g23 <- lmFit(t_g23_d) 
fit2.g23 <- ebayes(fit.g23)
temp23 <- fit2.g23$p
temp23[is.na(temp23)] <- 1


logit_df_1 <- cbind.data.frame(all_beta_d = rowMeans(cbind.data.frame(g4_d, g23_d), na.rm=T), fit2.all$t, temp_all)
save(logit_df_1, file = "results/logit_limma_1.RData")


logit_df_2d <- cbind.data.frame(g4_beta_d = rowMeans(g4_d, na.rm=T), fit2.g4$t, temp4)
save(logit_df_2d, file = "results/logit_limma_2d.RData")


logit_df_2a <- cbind.data.frame(g2_beta_d = rowMeans(g2_d, na.rm=T), fit2.g2$t, temp2)
save(logit_df_2a, file = "results/logit_limma_2a.RData")

logit_df_2b <- cbind.data.frame(g3_beta_d = rowMeans(g3_d, na.rm=T), fit2.g3$t, temp3)
save(logit_df_2b, file = "results/logit_limma_2b.RData")

logit_df_2c <- cbind.data.frame(g23_beta_d = rowMeans(g23_d, na.rm=T), fit2.g23$t, temp23)
save(logit_df_2c, file = "results/logit_limma_2c.RData")


logit_df <- cbind.data.frame(g4_beta_d = rowMeans(g4_d, na.rm=T), g2_beta_d = rowMeans(g2_d, na.rm=T),fit2.g4.g23$t, temp)
save(logit_df, file = "results/logit_limma_3.RData")



