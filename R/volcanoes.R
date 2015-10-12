
load('results/log_limma_1.RData')
load('results/logit_limma_1.RData')


p_cutoff <- .05

load("data/annotations_with_gencode_1500_TSS_1000.RData", verbose=TRUE); anno <- anno_with_k_tss_p_named;

log_df_1 <- subset(log_df_1, rownames(log_df_1) %in% anno$genes[anno$IlmnID %in% rownames(logit_df_1)])

methyl_sf <- (nrow(logit_df_1)/nrow(log_df_1))


cols <- c("gainsboro", rgb(250, 164, 26, max = 255), rgb(81, 184, 72, max = 255), "grey25", "mediumorchid4", rgb(245, 126, 32, max=255), rgb(10, 140, 68, max = 255))


###############   Model 1   ###############
names(log_df_1)[3] <- "temp"
names(logit_df_1)[3] <- "temp"

cond1 <- which(logit_df_1$temp*methyl_sf < p_cutoff & (logit_df_1[,1] <= -.2) )
cond2 <- which(logit_df_1$temp*methyl_sf < p_cutoff & (logit_df_1[,1] >= .2) )
length(cond1); length(cond2)

pdf(paste0("results/plots/m1_methyl_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(logit_df_1[,1], -log10(logit_df_1$temp), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~beta)), type = "n", main = "450k Methylation Array: Model 1", xlim = c(-1,1), yaxt = 'n')
axis(2, las = 1)

with(logit_df_1[-c(cond1,cond2),], points(all_beta_d, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(logit_df_1[cond1,], points(all_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[2] ))  
with(logit_df_1[cond2,], points(all_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-.2,.2), col=cols[6:7], lwd = 2)
abline(lty=5, h = -log(p_cutoff/methyl_sf, 10), col=c(cols[4]), lwd = 2)

dev.off()



cond1 <- which(log_df_1$temp < p_cutoff & log_df_1$all_fc_log_2 > log(2,2))
cond2 <- which(log_df_1$temp < p_cutoff & log_df_1$all_fc_log_2 < log(1/2,2))

pdf(paste0("results/plots/m1_rna_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

xl <- max(abs(log_df_1[,1]))

plot(log_df_1[,1], -log10(log_df_1[,3]), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~log[2]~FPKM)), type = "n", main = "RNA-seq: Model 1", yaxt = 'n', xlim=c(-xl,xl))
axis(2, las = 1)

with(log_df_1[-c(cond1, cond2),], points(all_fc_log_2, -log10(temp), pch = 20, cex = .7, col = cols[1])) 
with(log_df_1[cond1,], points(all_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[2] )) 
with(log_df_1[cond2,], points(all_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[3] )) 
abline(lty=2, v = c(log(1/2,2),log(2,2)), col=cols[7:6], lwd = 2)
abline(lty=5, h = -log10(p_cutoff), col=c(cols[4]), lwd = 2)
dev.off()



###############   Model 2   ###############
load('results/log_limma_2a.RData')
load('results/logit_limma_2a.RData')
load('results/log_limma_2b.RData')
load('results/logit_limma_2b.RData')
load('results/log_limma_2c.RData')
load('results/logit_limma_2c.RData')
load('results/log_limma_2d.RData')
load('results/logit_limma_2d.RData')


log_df_2a <- subset(log_df_2a, rownames(log_df_2a) %in% anno$genes[anno$IlmnID %in% rownames(logit_df_2a)])
log_df_2b <- subset(log_df_2b, rownames(log_df_2b) %in% anno$genes[anno$IlmnID %in% rownames(logit_df_2b)])
log_df_2c <- subset(log_df_2c, rownames(log_df_2c) %in% anno$genes[anno$IlmnID %in% rownames(logit_df_2c)])
log_df_2d <- subset(log_df_2d, rownames(log_df_2d) %in% anno$genes[anno$IlmnID %in% rownames(logit_df_2d)])


names(log_df_2a)[3] <- "temp"
names(logit_df_2a)[3] <- "temp"
names(log_df_2b)[3] <- "temp"
names(logit_df_2b)[3] <- "temp"
names(log_df_2c)[3] <- "temp"
names(logit_df_2c)[3] <- "temp"
names(log_df_2d)[3] <- "temp"
names(logit_df_2d)[3] <- "temp"


methyl_sf <- (nrow(logit_df_2a)/nrow(log_df_2a))


cond1 <- which(logit_df_2a$temp*methyl_sf < p_cutoff & (logit_df_2a[,1] <= -.2) )
cond2 <- which(logit_df_2a$temp*methyl_sf < p_cutoff & (logit_df_2a[,1] >= .2) )
length(cond1); length(cond2)
#m2b.hypo <- rownames(logit_df_2b)[cond1]; m2b.hyper <- rownames(logit_df_2b)[cond2]; save(m2b.hypo, m2b.hyper, file="~/Desktop/models_probes/m2b.meth.Rdata")

pdf(paste0("results/plots/m2a_methyl_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(logit_df_2d[,1], logit_df_2c[,1], logit_df_2b[,1], logit_df_2a[,1]), -log10(c(logit_df_2d$temp, logit_df_2c$temp, logit_df_2b$temp, logit_df_2a$temp)), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~beta)), type = "n", main = "450k Methylation Array: Grade 2", xlim = c(-1,1), yaxt = 'n')
axis(2, las = 1)

with(logit_df_2a[-c(cond1,cond2),], points(g2_beta_d, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(logit_df_2a[cond1,], points(g2_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[2] ))  
with(logit_df_2a[cond2,], points(g2_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-.2,.2), col=cols[6:7], lwd = 2)
abline(lty=5, h = -log10(p_cutoff/methyl_sf), col=c(cols[4]), lwd = 2)

dev.off()

cond1 <- which(logit_df_2b$temp*methyl_sf < p_cutoff & (logit_df_2b[,1] <= -.2) )
cond2 <- which(logit_df_2b$temp*methyl_sf < p_cutoff & (logit_df_2b[,1] >= .2) )
length(cond1); length(cond2)
#m2b.hypo <- rownames(logit_df_2b)[cond1]; m2b.hyper <- rownames(logit_df_2b)[cond2]; save(m2b.hypo, m2b.hyper, file="~/Desktop/models_probes/m2b.meth.Rdata")

pdf(paste0("results/plots/m2b_methyl_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(logit_df_2d[,1], logit_df_2c[,1], logit_df_2b[,1], logit_df_2a[,1]), -log10(c(logit_df_2d$temp, logit_df_2c$temp, logit_df_2b$temp, logit_df_2a$temp)), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~beta)), type = "n", main = "450k Methylation Array: Grade 3", xlim = c(-1,1), yaxt = 'n')
axis(2, las = 1)

with(logit_df_2b[-c(cond1,cond2),], points(g3_beta_d, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(logit_df_2b[cond1,], points(g3_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[2] ))  
with(logit_df_2b[cond2,], points(g3_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-.2,.2), col=cols[6:7], lwd = 2)
abline(lty=5, h = -log10(p_cutoff/methyl_sf), col=c(cols[4]), lwd = 2)

dev.off()


cond1 <- which(logit_df_2c$temp*methyl_sf < p_cutoff & (logit_df_2c[,1] <= -.2) )
cond2 <- which(logit_df_2c$temp*methyl_sf < p_cutoff & (logit_df_2c[,1] >= .2) )
length(cond1); length(cond2)
#m2b.hypo <- rownames(logit_df_2b)[cond1]; m2b.hyper <- rownames(logit_df_2b)[cond2]; save(m2b.hypo, m2b.hyper, file="~/Desktop/models_probes/m2b.meth.Rdata")

pdf(paste0("results/plots/m2c_methyl_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(logit_df_2d[,1], logit_df_2c[,1], logit_df_2b[,1], logit_df_2a[,1]), -log10(c(logit_df_2d$temp, logit_df_2c$temp, logit_df_2b$temp, logit_df_2a$temp)), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~beta)), type = "n", main = "450k Methylation Array: Grade 2 & 3", xlim = c(-1,1), yaxt = 'n')
axis(2, las = 1)

with(logit_df_2c[-c(cond1,cond2),], points(g23_beta_d, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(logit_df_2c[cond1,], points(g23_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[2] ))  
with(logit_df_2c[cond2,], points(g23_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-.2,.2), col=cols[6:7], lwd = 2)
abline(lty=5, h = -log10(p_cutoff/methyl_sf), col=c(cols[4]), lwd = 2)

dev.off()



cond1 <- which(logit_df_2d$temp*methyl_sf < p_cutoff & (logit_df_2d[,1] <= -.2) )
cond2 <- which(logit_df_2d$temp*methyl_sf < p_cutoff & (logit_df_2d[,1] >= .2) )
length(cond1); length(cond2)
#m2b.hypo <- rownames(logit_df_2b)[cond1]; m2b.hyper <- rownames(logit_df_2b)[cond2]; save(m2b.hypo, m2b.hyper, file="~/Desktop/models_probes/m2b.meth.Rdata")

pdf(paste0("results/plots/m2d_methyl_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(logit_df_2d[,1], logit_df_2c[,1], logit_df_2b[,1], logit_df_2a[,1]), -log10(c(logit_df_2d$temp, logit_df_2c$temp, logit_df_2b$temp, logit_df_2a$temp)), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~beta)), type = "n", main = "450k Methylation Array: Grade 4", xlim = c(-1,1), yaxt = 'n')
axis(2, las = 1)

with(logit_df_2d[-c(cond1,cond2),], points(g4_beta_d, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(logit_df_2d[cond1,], points(g4_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[2] ))  
with(logit_df_2d[cond2,], points(g4_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-.2,.2), col=cols[6:7], lwd = 2)
abline(lty=5, h = -log10(p_cutoff/methyl_sf), col=c(cols[4]), lwd = 2)

dev.off()




cond1 <- which(log_df_2a$temp < p_cutoff & log_df_2a$g2_fc_log_2 > log(2,2))
cond2 <- which(log_df_2a$temp < p_cutoff & log_df_2a$g2_fc_log_2 < log(1/2,2))
length(cond1); length(cond2)
xl <- max(abs(log_df_2d[,1]), abs(log_df_2d[,1]), abs(log_df_2b[,1]), abs(log_df_2a[,1]))

pdf(paste0("results/plots/m2a_rna_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(log_df_2d[,1], log_df_2c[,1], log_df_2b[,1], log_df_2a[,1]), -log10(c(log_df_2d[,3], log_df_2c[,3], log_df_2b[,3], log_df_2a[,3])), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~log[2]~FPKM)), type = "n", main = "RNA-seq: Grade 2", yaxt = 'n', xlim=c(-xl,xl))
axis(2, las = 1)

with(log_df_2a[-c(cond1, cond2),], points(g2_fc_log_2, -log10(temp), pch = 20, cex = .7, col = cols[1])) 
with(log_df_2a[cond1,], points(g2_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[2] )) 
with(log_df_2a[cond2,], points(g2_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[3] )) 
abline(lty=2, v = c(log(1/2,2),log(2,2)), col=cols[7:6], lwd = 2)
abline(lty=5, h = -log10(p_cutoff), col=c(cols[4]), lwd = 2)
dev.off()


cond1 <- which(log_df_2b$temp < p_cutoff & log_df_2b$g3_fc_log_2 > log(2,2))
cond2 <- which(log_df_2b$temp < p_cutoff & log_df_2b$g3_fc_log_2 < log(1/2,2))
length(cond1); length(cond2)
xl <- max(abs(log_df_2d[,1]), abs(log_df_2d[,1]), abs(log_df_2b[,1]), abs(log_df_2a[,1]))

pdf(paste0("results/plots/m2b_rna_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(log_df_2d[,1], log_df_2c[,1], log_df_2b[,1], log_df_2a[,1]), -log10(c(log_df_2d[,3], log_df_2c[,3], log_df_2b[,3], log_df_2a[,3])), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~log[2]~FPKM)), type = "n", main = "RNA-seq: Grade 3", yaxt = 'n', xlim=c(-xl,xl))
axis(2, las = 1)

with(log_df_2b[-c(cond1, cond2),], points(g3_fc_log_2, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(log_df_2b[cond1,], points(g3_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[2] )) 
with(log_df_2b[cond2,], points(g3_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[3] )) 
abline(lty=2, v = c(log(1/2,2),log(2,2)), col=cols[7:6], lwd = 2)
abline(lty=5, h = -log10(p_cutoff), col=c(cols[4]), lwd = 2)
dev.off()


cond1 <- which(log_df_2c$temp < p_cutoff & log_df_2c$g23_fc_log_2 > log(2,2))
cond2 <- which(log_df_2c$temp < p_cutoff & log_df_2c$g23_fc_log_2 < log(1/2,2))
length(cond1); length(cond2)
xl <- max(abs(log_df_2d[,1]), abs(log_df_2d[,1]), abs(log_df_2b[,1]), abs(log_df_2a[,1]))

pdf(paste0("results/plots/m2c_rna_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(log_df_2d[,1], log_df_2c[,1], log_df_2b[,1], log_df_2a[,1]), -log10(c(log_df_2d[,3], log_df_2c[,3], log_df_2b[,3], log_df_2a[,3])), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~log[2]~FPKM)), type = "n", main = "RNA-seq: Grade 2 & 3", yaxt = 'n', xlim=c(-xl,xl))
axis(2, las = 1)

with(log_df_2c[-c(cond1, cond2),], points(g23_fc_log_2, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(log_df_2c[cond1,], points(g23_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[2] )) 
with(log_df_2c[cond2,], points(g23_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[3] )) 
abline(lty=2, v = c(log(1/2,2),log(2,2)), col=cols[7:6], lwd = 2)
abline(lty=5, h = -log10(p_cutoff), col=c(cols[4]), lwd = 2)
dev.off()



cond1 <- which(log_df_2d$temp < p_cutoff & log_df_2d$g4_fc_log_2 > log(2,2))
cond2 <- which(log_df_2d$temp < p_cutoff & log_df_2d$g4_fc_log_2 < log(1/2,2))
length(cond1); length(cond2)
xl <- max(abs(log_df_2d[,1]), abs(log_df_2d[,1]), abs(log_df_2b[,1]), abs(log_df_2a[,1]))

pdf(paste0("results/plots/m2d_rna_or_", length(cond1), "_gr_", length(cond2), ".pdf"), useDingbats=FALSE)

plot(c(log_df_2d[,1], log_df_2c[,1], log_df_2b[,1], log_df_2a[,1]), -log10(c(log_df_2d[,3], log_df_2c[,3], log_df_2b[,3], log_df_2a[,3])), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(bar(Delta~log[2]~FPKM)), type = "n", main = "RNA-seq: Grade 4", yaxt = 'n', xlim=c(-xl,xl))
axis(2, las = 1)

with(log_df_2d[-c(cond1, cond2),], points(g4_fc_log_2, -log10(temp), pch = 20, cex = .7, col = cols[1]))
with(log_df_2d[cond1,], points(g4_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[2] )) 
with(log_df_2d[cond2,], points(g4_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[3] )) 
abline(lty=2, v = c(log(1/2,2),log(2,2)), col=cols[7:6], lwd = 2)
abline(lty=5, h = -log10(p_cutoff), col=c(cols[4]), lwd = 2)
dev.off()




###############   Model 3   ###############

load('results/log_limma_3.RData')
load('results/logit_limma_3.RData')
sig_genes_hypeMehtyl_upExpr <- scan('results/M3_genes.txt', what = "")

log_df <- subset(log_df, rownames(log_df) %in% anno$genes[anno$IlmnID %in% rownames(logit_df)])

methyl_sf <- (nrow(logit_df)/nrow(log_df))


cond1 <- with(logit_df,which(((temp*methyl_sf) < p_cutoff) & (g4_beta_d <= -.2) & (g4_beta_d - g2_beta_d) < -.15))
cond2 <- with(logit_df,which((temp*methyl_sf) < p_cutoff & (g4_beta_d >= .2)  & (g4_beta_d - g2_beta_d) > .15))
cond3 <- with(logit_df,which(row.names(logit_df) %in% anno$IlmnID[anno$gene_name %in% sig_genes_hypeMehtyl_upExpr] & ((temp*methyl_sf) < p_cutoff) & (g4_beta_d <= -.2) & (g4_beta_d - g2_beta_d) < -.15))
length(cond1);length(cond2);length(cond3)

cat(row.names(logit_df)[cond3], sep = '\n', file = 'results/purple_probes.txt')
#m3.hypo <- rownames(logit_df)[cond1]; m3.hyper <- rownames(logit_df)[cond2]; save(m3.hypo, m3.hyper, file="~/Desktop/models_probes/m3.meth.Rdata")

pdf(paste0("results/plots/m3_methyl_or_", length(cond1), "_gr_", length(cond2), "_purp_", length(cond3), ".pdf"), useDingbats=FALSE)

with(logit_df[-c(cond1,cond2,cond3),], plot(g4_beta_d - g2_beta_d, -log10(temp), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(Delta~bar(Delta~beta)), col = cols[1], main = "450k Methylation Array", xlim = c(-1,1), ylim=c(-1,10), yaxt = 'n')) 
axis(2, las = 1)
with(logit_df[cond1,], points(g4_beta_d - g2_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[2] ))  
with(logit_df[cond2,], points(g4_beta_d - g2_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-.15,.15), col=cols[6:7], lwd = 2)
abline(lty=5, h = -log10((p_cutoff/methyl_sf)), col=cols[4], lwd = 2)
with(logit_df[cond3,], points(g4_beta_d - g2_beta_d, -log10(temp), pch = 17, cex = .74, col =cols[5] ))  
dev.off()




cond1 <- with(log_df, which(temp < p_cutoff & (g4_fc_log_2 >= log(2,2)) & (g4_fc_log_2 - g2_fc_log_2) >1 ))
cond2 <- with(log_df, which(temp < p_cutoff & (g4_fc_log_2 <=log(1/2,2)) &  (g4_fc_log_2 - g2_fc_log_2) < -1))
cond3 <- with(log_df, which(temp < p_cutoff & (g4_fc_log_2 >= log(2,2)) & (g4_fc_log_2 - g2_fc_log_2) >1  & (row.names(log_df) %in% anno$genes[anno$gene_name %in% sig_genes_hypeMehtyl_upExpr])))
length(cond1);length(cond2);length(cond3)
xl <- max(abs(log_df$g4_fc_log_2 - log_df$g2_fc_log_2))

pdf(paste0("results/plots/m3_rna_or_", length(cond1), "_gr_", length(cond2), "_purp_", length(cond3), ".pdf"), useDingbats=FALSE)

with(log_df, plot(g4_fc_log_2 - g2_fc_log_2, -log10(temp), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(Delta~bar(Delta~log[2]~FPKM)), col = cols[1], main = "RNA-seq", yaxt = 'n',type = 'n', xlim=c(-xl,xl))) 
axis(2, las = 1)
with(log_df[-c(cond1,cond2,cond3),], points(g4_fc_log_2 - g2_fc_log_2, -log10(temp), pch = 20, cex = .7, col = cols[1])) 

with(log_df[cond1, ], points(g4_fc_log_2 - g2_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[2] )) 
with(log_df[cond2, ], points(g4_fc_log_2 - g2_fc_log_2, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-1,1), col=cols[7:6], lwd = 2)
abline(lty=5, h = -log(p_cutoff, 10), col=c(cols[4]), lwd = 2)
with(log_df[cond3, ], points(g4_fc_log_2 - g2_fc_log_2, -log10(temp), pch = 17, cex = .74, col =cols[5] )) 
dev.off()




###############   Supplement Model 3   ###############

load('results/log_limma_3.RData')
load('results/logit_limma_3.RData')
sig_genes_hypeMehtyl_upExpr <- scan('results/M3_genes.txt', what = "")

log_df <- subset(log_df, rownames(log_df) %in% anno$genes[anno$IlmnID %in% rownames(logit_df)])

methyl_sf <- (nrow(logit_df)/nrow(log_df))


cond1 <- with(logit_df,which(((temp*methyl_sf) < p_cutoff) & (g4_beta_d <= -.2) & (g4_beta_d - g2_beta_d) < -.15))
cond2 <- with(logit_df,which((temp*methyl_sf) < p_cutoff & (g4_beta_d >= .2)  & (g4_beta_d - g2_beta_d) > .15))
cond3 <- with(logit_df,which(row.names(logit_df) %in% anno$IlmnID[anno$gene_name %in% sig_genes_hypeMehtyl_upExpr] ))
length(cond1);length(cond2);length(cond3)
#m3.hypo <- rownames(logit_df)[cond1]; m3.hyper <- rownames(logit_df)[cond2]; save(m3.hypo, m3.hyper, file="~/Desktop/models_probes/m3.meth.Rdata")

pdf(paste0("results/plots/m3_sup5d_methyl_or_", length(cond1), "_gr_", length(cond2), "_purp_", length(cond3), ".pdf"), useDingbats=FALSE)

with(logit_df[-c(cond1,cond2,cond3),], plot(g4_beta_d - g2_beta_d, -log10(temp), pch = 20, cex = .7, ylab = expression(-log[10]~p), xlab = expression(Delta~bar(Delta~beta)), col = cols[1], main = "450k Methylation Array", xlim = c(-1,1), ylim=c(-1,10), yaxt = 'n')) 
axis(2, las = 1)
with(logit_df[cond1,], points(g4_beta_d - g2_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[2] ))  
with(logit_df[cond2,], points(g4_beta_d - g2_beta_d, -log10(temp), pch = 20, cex = .72, col =cols[3] ))  
abline(lty=2, v = c(-.15,.15), col=cols[6:7], lwd = 2)
abline(lty=5, h = -log10((p_cutoff/methyl_sf)), col=cols[4], lwd = 2)
with(logit_df[cond3,], points(g4_beta_d - g2_beta_d, -log10(temp), pch = 17, cex = .74, col =cols[5] ))  
dev.off()




cond3 <- with(logit_df,which(row.names(logit_df) %in% anno$IlmnID[anno$gene_name %in% sig_genes_hypeMehtyl_upExpr] & ((temp*methyl_sf) < p_cutoff) & (g4_beta_d <= -.2) & (g4_beta_d - g2_beta_d) < -.15))


pdf(paste0("results/plots/m3_sup5c_methyl_or_", length(cond1), "_gr_", length(cond2), "_purp_", length(cond3), ".pdf"), useDingbats=FALSE)


with(logit_df, plot(g2_beta_d, g4_beta_d, pch = 20, cex = .7, xlab = expression(bar(Delta~beta[non-gbm])), ylab = expression(bar(Delta~beta[gbm])), col = cols[1], main = "450k Methylation Array", xlim = c(-1,1), ylim=c(-1,1), yaxt = 'n')) 
axis(2, las = 1)
with(logit_df[cond1,], points(g2_beta_d,g4_beta_d, pch = 20, cex = .75, col =cols[2] ))  
with(logit_df[cond2,], points(g2_beta_d,g4_beta_d, pch = 20, cex = .75, col =cols[3] ))  
#abline(lty=2, h = c(-.2,.2), col=cols[6:7], lwd = 2)
#abline(lty=2, a = c(-.15), b=1, col=cols[6], lwd = 2)
#abline(lty=2, a = c(.15), b=1, col=cols[7], lwd = 2)
abline(lty=2, a = 0, b=1, col=cols[1], lwd = 1.5)

segments(-1.5,.2, .05,lty=2,col = cols[7], lwd = 2)
segments(.05,.2,1.5, 1.65,lty=2,col = cols[7], lwd = 2)
segments(-.05,-.2,-1.5, -1.65,lty=2,col = cols[6], lwd = 2)
segments(1.5,-.2, -.05,lty=2,col = cols[6], lwd = 2)

with(logit_df[cond3,], points(g2_beta_d, g4_beta_d, pch = 17, cex = .75, col =cols[5] ))  
dev.off()

cond1 <- with(log_df, which(temp < p_cutoff & (g4_fc_log_2 >= log(2,2)) & (g4_fc_log_2 - g2_fc_log_2) >1 ))
cond2 <- with(log_df, which(temp < p_cutoff & (g4_fc_log_2 <=log(1/2,2)) &  (g4_fc_log_2 - g2_fc_log_2) < -1))
cond3 <- with(log_df, which(temp < p_cutoff & (g4_fc_log_2 >= log(2,2)) & (g4_fc_log_2 - g2_fc_log_2) >1  & (row.names(log_df) %in% anno$genes[anno$gene_name %in% sig_genes_hypeMehtyl_upExpr])))

#names(log_df)[1:2] <- c("g4_beta_diff", "g23_beta_diff")

pdf(paste0("results/plots/m3_sup5c_rna_or_", length(cond1), "_gr_", length(cond2), "_purp_", length(cond3), ".pdf"), useDingbats=FALSE)

with(log_df, plot(g2_fc_log_2, g4_fc_log_2, pch = 20, cex = .7, ylab = expression(bar(Delta~log[2]~FPKM[gbm])), xlab = expression(bar(Delta~log[2]~FPKM[non-gbm])), col = cols[1], main = "RNA-seq", yaxt = 'n', ylim = c(-1,1)*7, xlim = c(-1,1)*7)) 
axis(2, las = 1)
with(log_df[cond1, ], points(g2_fc_log_2, g4_fc_log_2, pch = 20, cex = .75, col =cols[2] )) 
with(log_df[cond2, ], points(g2_fc_log_2, g4_fc_log_2, pch = 20, cex = .75, col =cols[3] ))  
#abline(lty=2, h = c(-1,1), col=cols[7:6], lwd = 2)
#abline(lty=2, a = 1, b = 1, col = cols[6], lwd = 2)
#abline(lty=2, a = -1, b = 1, col = cols[7], lwd = 2)
abline(lty=2, a = 0, b=1, col=cols[1], lwd = 1.5)

segments(-8, 1, 0,lty=2,col = cols[6], lwd = 2 )
segments(0, 1, 8,9,lty=2,col = cols[6], lwd = 2 )
segments(8, -1, 0,lty=2,col = cols[7], lwd = 2 )
segments(0, -1, -8,-9,lty=2,col = cols[7], lwd = 2 )

with(log_df[cond3, ], points(g2_fc_log_2, g4_fc_log_2, pch = 17, cex = .75, col =cols[5] )) 
dev.off()
