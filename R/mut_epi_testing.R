

phylogen_muts_files <- list.files('results/phylogenetic_tables/', '^Pat.*', full.names=T)

muts_per_file <- lapply(phylogen_muts_files, function(fname){
  read.table(fname, T, sep='\t',  as.is = T)
})

filt_muts_per_file <- lapply(muts_per_file, function(x) subset(x, type != 'Silent'))

filt_nonub_muts_per_file <- lapply(filt_muts_per_file, function(x){
  n <- length(unique(x$sample_type))
  uniq_mut <- apply(x, 1,  function(y) paste(y[1:5], collapse = '_') )
  x[which(uniq_mut %in% names(which(table(uniq_mut) != n))),]
})



library(reshape2)

filt_muts_per_file_wide <- lapply(filt_nonub_muts_per_file, function(muts_per_file ) dcast(muts_per_file, X.gene + contig + position + ref_allele+alt_allele ~ sample_type, fun.aggregate = length) )

names(filt_muts_per_file_wide) <- sapply(strsplit(basename(phylogen_muts_files),'\\.'), function(x) x[1])


filt_muts_per_file_wide[['Patient04']]$Primaryv5 <- rep(0, nrow(filt_muts_per_file_wide[['Patient04']]))
filt_muts_per_file_wide[['Patient04']]$Primaryv6 <- rep(0, nrow(filt_muts_per_file_wide[['Patient04']]))

filt_muts_per_file_wide_ordr <- lapply(filt_muts_per_file_wide, function(muts_per_file) muts_per_file[,c('X.gene', 'contig', 'position', 'ref_allele', 'alt_allele', sort(names(muts_per_file)[-(1:5)]))])

#names(filt_muts_per_file_wide_ordr) <- sapply(strsplit(basename(phylogen_muts_files),'\\.'), function(x) x[1])



### probes
probe_files <- list.files('results/phyloepi_tables', 'Pat.*', full.names=T)

probes_per_file <- lapply(probe_files, function(fname){
  read.table(fname, T, sep='\t', row.names = 1, as.is = T)[,-1]
})

load(file = 'data/annotations_with_gencode_1500_TSS_1000.RData'); anno <- anno_with_k_tss_p_named[, c('IlmnID', 'gene_name')]

probes_per_file_gene <- lapply(probes_per_file, function(probe_file){
  probes_per_file2 <- cbind.data.frame(IlmnID = row.names(probe_file), probe_file, stringsAsFactors = F)
  merge(probes_per_file2, anno)
})

names(probes_per_file_gene) <- gsub('.txt$','',basename(probe_files))

probes_per_file_ordr <- lapply(probes_per_file_gene, function(probes_file){ probes_file[, c('IlmnID','gene_name',sort(grep('Patient', names(probes_file), value = T)))]})



library(limma)
trans <- function(x) log(x/(1-x))

muts_prof_split <- lapply(names(filt_muts_per_file_wide_ordr), function(pat_name){
  muts <- filt_muts_per_file_wide_ordr[[pat_name]][,-c(1:5)]
  split(filt_muts_per_file_wide_ordr[[pat_name]], muts, drop = T)
})

names(muts_prof_split) <- names(filt_muts_per_file_wide_ordr)

methyl_limma_all <- lapply(names(filt_muts_per_file_wide_ordr), function(pat_name){
	cat(pat_name, '\n')
  probes <- probes_per_file_ordr[[pat_name]][,-c(1:2)]
  muts_profiles_uq <- sapply(strsplit(names(muts_prof_split[[pat_name]]),'\\.'), as.numeric)
  names(muts_profiles_uq) <- names(muts_prof_split[[pat_name]])
  apply(muts_profiles_uq, 2, function(prof){
  	design <- model.matrix(~factor(prof))
  	colnames(design) <- c("intercept", "slope") 
  	fit_model <- lmFit(trans(probes), design) 
    fit_model2 <- ebayes(fit_model)
    
    fit_model2.pvalues <- fit_model2$p[,2]
    temp <- fit_model2.pvalues
    temp[is.na(temp)] <- 1

    cbind.data.frame(t_stat = fit_model2$t[,2], pval = temp)
  })
})

names(methyl_limma_all) <- names(filt_muts_per_file_wide_ordr)

methyl_limma_gene_int <- lapply(names(filt_muts_per_file_wide_ordr), function(pat_name){
  names(methyl_limma_all[[pat_name]]) <- names(muts_prof_split[[pat_name]])
  probe_profs <- lapply(names(muts_prof_split[[pat_name]]), function(mut_prof_name){
  	int_inds <- which(probes_per_file_ordr[[pat_name]]$gene_name %in% muts_prof_split[[pat_name]][[mut_prof_name]]$X.gene)
  	cbind.data.frame(probes_per_file_ordr[[pat_name]][int_inds,], methyl_limma_all[[pat_name]][[mut_prof_name]][int_inds,])
  })
  names(probe_profs) <- names(muts_prof_split[[pat_name]])
  return(probe_profs)
})

methyl_limma_gene_int <- lapply(methyl_limma_gene_int, function(x) x[which(sapply(x, nrow)>0)])
names(methyl_limma_gene_int) <- names(filt_muts_per_file_wide_ordr)

save(methyl_limma_gene_int, file = 'results/rev1_q2_methyl_limma_gene_int.RData')
