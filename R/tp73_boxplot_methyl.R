load('data/sample_inds.RData')
load('data/data.clean.Rdata')

tp73_probes <- c('cg03846767', 'cg04493946', 'cg16607065', 'cg18873878', 'cg24878868', 'cg25115460', 'cg26208930')

cat(tp73_probes, sep = '\n', file = 'data/tp73_probes.txt')

data.clean <- subset(data.clean, rownames(data.clean) %in% tp73_probes)

g4_p <- data.clean[,pri.v1.pairs.r1g4]
g23_p <- data.clean[,pri.v1.pairs.r1g23]

g4_r <- data.clean[,rec1.v1.pairs.r1g4]
g23_r <- data.clean[,rec1.v1.pairs.r1g23]

g4_d <- g4_r - g4_p
g23_d <- g23_r - g23_p

full_dat <- rbind(t(g23_d), t(g4_d))
gbm_groups <- rep( c("II,III", "IV"), c(ncol(g23_d), ncol(g4_d)))

sig_probes <- scan(file = 'results/sig_probes.txt', what = "")
sig_tp73 <- c("", "_sig")[as.integer(tp73_probes %in% sig_probes)+1]

for(i in seq_along(tp73_probes)){
  pdf(paste0('results/plots/tp73_', tp73_probes[i], sig_tp73[i], '.pdf'), 8, 12)
  boxplot(full_dat[,i] ~ gbm_groups, pch = 17, col = c('green4','darkorange1'), ylab = expression(Delta~beta), main = tp73_probes[i], ylim = c(-.7, .3))
  dev.off()
}


