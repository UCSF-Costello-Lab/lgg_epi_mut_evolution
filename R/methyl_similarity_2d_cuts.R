if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

my_palette <- colorRampPalette(c('yellow2',"red"))(n = 10)

load('data/sample_inds.RData')
load('data/data.clean.Rdata')

sig_probes <- scan('results/sig_probes.txt', what = '')

all_p <- data.clean[,c(pri.v1.pairs.r1g23, pri.v1.pairs.r1g4)]
all_r <- data.clean[,c(rec1.v1.pairs.r1g23, rec1.v1.pairs.r1g4)]

p_sbst <- all_p[sig_probes,]
r_sbst <- all_r[sig_probes,]

names(all_p) <- gsub('_.*', '', names(all_p))
names(all_r) <- gsub('_.*', '', names(all_r))


pat <- names(all_p)[1]
p_mat <- as.matrix(table(cut(all_p[,pat], seq(0, 1, .2)), cut(all_r[,pat], seq(0, 1, .2))))


heatmap.2(log10(p_mat), cellnote = p_mat,
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
#  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  Colv="NA", scale = 'none')            # turn off column clustering



p_mat <- as.matrix(table(cut(rowMeans(p_sbst), seq(0, 1, .2)), cut(rowMeans(r_sbst), seq(0, 1, .2))))
heatmap.2(log10(p_mat+1), cellnote = p_mat,
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
#  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  Colv="NA", scale = 'none')            # turn off column clustering

