cols <- c("gainsboro", rgb(250, 164, 26, max = 255), rgb(81, 184, 72, max = 255), "grey25", "mediumorchid4", rgb(245, 126, 32, max=255), rgb(10, 140, 68, max = 255), "blue")

xlab_cols_map <- c(rgb(55, 126, 184, max = 255), rgb(25, 174, 74, max = 255), rgb(194, 48, 28, max = 255), rgb(152, 78, 163, max = 255))

#load('data/sample_inds.RData')
load('data/data.clean.Rdata')

cc <- 0.4

tree.pats <- c("Patient17","Patient18", "Patient04", "Patient01","Patient49","Patient90")
header <- tree.pats[4]

pat_inds <- grep(header, names(data.clean))
#pat_sub <- data.clean[,pat_inds]

probe_inds <- which(apply(data.clean[,pat_inds], 1, function(x) max(x)-min(x)) > cc)
dat <- data.clean[probe_inds,c(which(names(data.clean) == 'A141_Insula'), pat_inds)]
dat <- dat[-which(apply(dat, 1, function(x) any(is.na(x)))), ]
names(dat)[1] <- "Normal Brain"

write.table(dat, file = file.path('results/phyloepi_tables/', paste0(header, '.txt')), row.names = T, sep = '\t')


d_m <- as.matrix(dat)
dat_scaled <- (d_m - mean(d_m))
tst <- svd(dat_scaled)
g <- with(tst, u%*%diag(length(d)))
h <- with(tst, diag(d) %*% t(v))

x <- g
y <- t(h)
xlabs <- names(dat)

xlab_cols <- c(cols[4], xlab_cols_map[c(1,1,1,1,2,3,2)])

p1 <- 1; t1 <- 450; p2 <- 2; t2 <- 0
#cat(row.names(dat[which(rank(-abs(x[,p1])) <= t1),]), sep = '\n', file = paste0('p01_PC', p1, '_',t1, '.txt'))
#cat(row.names(dat[which(rank(-abs(x[,p2])) <= t2),]), sep = '\n', file = paste0('p01_PC', p2, '_',t2, '.txt'))
p1 <- 3; t1 <- 0;   p2 <- 4; t2 <- 50
#cat(row.names(dat[which(rank(-abs(x[,p1])) <= t1),]), sep = '\n', file = paste0('p01_PC', p1, '_',t1, '.txt'))
#cat(row.names(dat[which(rank(-abs(x[,p2])) <= t2),]), sep = '\n', file = paste0('p01_PC', p2, '_',t2, '.txt'))
p1 <- 5; t1 <- 25;   p2 <- 6; t2 <- 0
#cat(row.names(dat[which(rank(-abs(x[,p1])) <= t1),]), sep = '\n', file = paste0('p01_PC', p1, '_',t1, '.txt'))
#cat(row.names(dat[which(rank(-abs(x[,p2])) <= t2),]), sep = '\n', file = paste0('p01_PC', p2, '_',t2, '.txt'))
p1 <- 7; t1 <- 0;   p2 <- 8; t2 <- 5
#cat(row.names(dat[which(rank(-abs(x[,p1])) <= t1),]), sep = '\n', file = paste0('p01_PC', p1, '_',t1, '.txt'))
#cat(row.names(dat[which(rank(-abs(x[,p2])) <= t2),]), sep = '\n', file = paste0('p01_PC', p2, '_',t2, '.txt'))

pal <- colorRampPalette(c("blue","gainsboro")); colsHM <- rev(pal(20))
pal <- colorRampPalette(c("red","blue")); colsHM <- rev(pal(20))
if(t1 > 0){
pdf(paste0('results/plots/svd_plots/',header, ".heat.",p1,'.',t1,".pdf"), useDingbats = F)
heatmap(as.matrix(dat[which(rank(-abs(x[,p1])) <= t1),]), Colv=NA, col=colsHM, scale = 'none')
dev.off()
}
if(t2>0){
pdf(paste0('results/plots/svd_plots/',header, ".heat.",p2,'.',t2,".pdf"), useDingbats = F)
heatmap(as.matrix(dat[which(rank(-abs(x[,p2])) <= t2),]), Colv=NA, col=colsHM, scale = 'none')
dev.off()
}

unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), abs(max(x, na.rm = TRUE)))
rangx1 <- unsigned.range(x[, p1])
rangx2 <- unsigned.range(x[, p2])
rangy1 <- unsigned.range(y[, p1])
rangy2 <- unsigned.range(y[, p2])

xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
ratio <- max(rangy1/rangx1, rangy2/rangx2)/1
pdf(paste0('results/plots/svd_plots/',header, ".SVD.",p1,".",p2,".pdf"), useDingbats = F)
plot(x[,p1], x[,p2], xlim = xlim, ylim = ylim, col = cols[1], cex = .4, pch = 20, xlab = paste0("PC",p1), ylab = paste0("PC",p2))
abline(h=0, v=0, lty=2, lwd = 1.3, col = "grey")
top_1_ind <- which(rank(-abs(x[,p1]))<=t1)
points(x[top_1_ind,p1],x[top_1_ind,p2], pch = 17, cex = .6, col = cols[4])

top_2_ind <- which(rank(-abs(x[,p2]))<=t2)
points(x[top_2_ind,p1],x[top_2_ind,p2], pch = 17, cex = .6, col = cols[4])

#text(x, xlabs, cex = 1, col = 'grey')
par(new = TRUE)
#dev.hold()
#on.exit(dev.flush(), add = TRUE)
plot(y[,p1], y[,p2], axes = FALSE, xlim = xlim * ratio, ylim = ylim * ratio, xlab = "", ylab = "", col = "blue", type = 'n')
axis(3)
axis(4)
box()
text(y[,p1]*1.1, y[,p2]*1.1, labels = xlabs, col = xlab_cols) 
arrows(0, 0, y[, p1], y[, p2], col = xlab_cols, length = .1, angle = 20)
dev.off()
