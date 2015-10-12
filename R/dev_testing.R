load('data/sample.info.RData')
load('data/data.clean.Rdata')

trans <- function(x) log(x/(1-x))

x7brain.adult <- which(sample.info$grade=="Brain" & sample.info$grade_of_rec1 == "Adult")
x31brain.fetal <- which(sample.info$grade=="Brain" & sample.info$grade_of_rec1 == "Fetal" & !(sample.info$patient_ID == "HuF02" & sample.info$sample_type!="FGM") ) 

b_adult <- data.clean[, x8brain.adult]
b_fetal <- data.clean[, x33brain.fetal]
asq_adult <- trans(b_adult)
asq_fetal <- trans(b_fetal)

library(limma)
design <- model.matrix(~factor(c(rep(1,ncol(asq_adult)), rep(2,ncol(asq_fetal))))) 
colnames(design) <- c("intercept", "slope") 
eset <- cbind.data.frame(asq_adult, asq_fetal)
asq_lmfit <- lmFit(eset, design) 
asq_ebayes <- ebayes(asq_lmfit)
pvals <- asq_ebayes$p[,2]
temp <- pvals
temp[is.na(temp)] <- 1

dev_logit_df <- cbind.data.frame(b_adult_mean = rowMeans(b_adult, na.rm=T), b_fetal_mean = rowMeans(b_fetal, na.rm=T) ,asq_ebayes$t[,2], temp)

save(dev_logit_df, file = "results/dev_logit_limma.RData")
