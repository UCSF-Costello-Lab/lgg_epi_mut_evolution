
p_cutoff <- .05
go_p_cutoff <- .001

load("data/annotations_with_gencode_1500_TSS_1000.RData", verbose=TRUE); anno <- anno_with_k_tss_p_named;

probe_bg <- scan(file = "results/probe_bg.txt", what = '')
rna_bg <- scan(file = "results/rna_bg.txt", what = '')

methyl_sf <- (length(probe_bg)/length(rna_bg))

m3_sig_probes <- scan(file = 'results/m3_hypo_probes.txt', what = "")

#### Development analysis ####
load('results/dev_logit_limma.RData')

## inds <- which((dev_logit_df$temp*methyl_sf) < p_cutoff)

inds <- with(dev_logit_df, which((temp*methyl_sf) < p_cutoff & ((b_fetal_mean - b_adult_mean) < -.2) ))

hypo_sig_in_dev <- intersect(m3_sig_probes, rownames(dev_logit_df[inds,]) )
sum(m3_sig_probes %in% rownames(dev_logit_df[inds,]) )

boot_vec <- replicate(10000, sum(m3_sig_probes %in% sample(probe_bg, length(inds))))
boot_vec2 <- replicate(10000, sum(sample(probe_bg, length(m3_sig_probes)) %in% rownames(dev_logit_df[inds,])))

p_dev1 <- sum(length(hypo_sig_in_dev) < boot_vec)
p_dev2 <- sum(length(hypo_sig_in_dev) < boot_vec2)

cat(hypo_sig_in_dev, sep='\n', file = 'results/hypo_sig_in_dev_probes.txt')


cat(setdiff(unique(subset(anno, IlmnID %in% hypo_sig_in_dev)$gene_name), NA), sep='\n', file = 'results/hypo_sig_in_dev_genes.txt')

#######  Rob's annotations  #######
load('data/annotations_with_robs_chip.RData')
##anno_with_robs_chip[,145:153]
## anno_with_robs_chip
anno_with_robs_chip <- subset(anno_with_robs_chip, IlmnID %in% probe_bg)
sig_vec <- as.integer(anno_with_robs_chip$IlmnID %in% m3_sig_probes)
for(i in levels(anno_with_robs_chip[,145])[4]){
  cnt_state <- rowSums(anno_with_robs_chip[,145:153]==i)
  inds <- which(cnt_state>1 )
  reg_probes <- anno_with_robs_chip$IlmnID[inds]
  boot_vec <- replicate(10000, sum(sample(probe_bg, length(m3_sig_probes)) %in% reg_probes))
  sum( m3_sig_probes %in% reg_probes );range(boot_vec)
sum(sum( m3_sig_probes %in% reg_probes ) < boot_vec)
  
}

cnt_state2 <- rowSums(sapply(anno_with_robs_chip[,145:153], function(x) x %in% c("E1", "E3", "E4")))

inds <- which(cnt_state2>1 )
reg_probes <- anno_with_robs_chip$IlmnID[inds]
boot_vec <- replicate(10000, sum(sample(probe_bg, length(m3_sig_probes)) %in% reg_probes))
sum( m3_sig_probes %in% reg_probes );range(boot_vec);sum(sum( m3_sig_probes %in% reg_probes ) < boot_vec)


cnt_state2 <- rowSums(sapply(anno_with_robs_chip[,145:153], function(x) x %in% c("E1", "E4")))

inds <- which(cnt_state2>1 )
reg_probes <- anno_with_robs_chip$IlmnID[inds]
boot_vec <- replicate(10000, sum(sample(probe_bg, length(m3_sig_probes)) %in% reg_probes))
sum( m3_sig_probes %in% reg_probes );range(boot_vec);sum(sum( m3_sig_probes %in% reg_probes ) < boot_vec)


cnt_state2 <- rowSums(sapply(anno_with_robs_chip[,145:153], function(x) x %in% c("E3", "E4")))

inds <- which(cnt_state2>1 )
reg_probes <- anno_with_robs_chip$IlmnID[inds]
boot_vec <- replicate(10000, sum(sample(probe_bg, length(m3_sig_probes)) %in% reg_probes))
sum( m3_sig_probes %in% reg_probes );range(boot_vec);sum(sum( m3_sig_probes %in% reg_probes ) < boot_vec)



load('data/annotations_with_nearest_genes.RData')
cnt_state <- rowSums(anno_with_robs_chip[,145:153]=='E4')
inds <- which(cnt_state>1 & sig_vec )
anno_with_nearest_genes$gene_name[inds]
reg_probes <- anno_with_robs_chip$IlmnID[inds]
boot_vec <- replicate(10000, sum(sample(probe_bg, length(m3_sig_probes)) %in% reg_probes))
