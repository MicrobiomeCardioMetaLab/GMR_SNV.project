library(vegan)
library(dplyr)
library(parallel)


args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
pheno_id <- args[2]


pheno_df <- read.delim('PRJNA834801_phenotype_num_filter.xls', check.names = F)


snv_mtx <- read.delim(paste0(species, '/snv_matrix2base_filter.xls.gz'), check.names = F, row.names = 1)


sub_df <- pheno_df[, c('sample_id', pheno_id), drop=FALSE]

process_snv <- function(i) {
  tryCatch({
    
    print(paste0(i, " / ", nrow(snv_mtx)))
    
    sub_snv <- data.frame(t(snv_mtx[i, , drop=FALSE]), check.names = FALSE)
    sub_snv$sample_id <- sub("\\.[^.]+$", "", rownames(sub_snv))
    names(sub_snv) <- c('snv_site', 'sample_id')
    sub_snv <- sub_snv %>% filter(snv_site != '-')
    
    
    dat <- merge(sub_df, sub_snv, by='sample_id')
    names(dat) <- c('sample', 'pheno', 'snv_site')
    dat[dat==""] <- NA
    om_dat <- na.omit(dat)
    
    base_counts <- table(om_dat$snv_site)
    
    n_base1 <- base_counts[1]
    n_base2 <- base_counts[2] 
    
    ass_test <- wilcox.test(pheno ~ snv_site, data = om_dat)
    wil_p <- ass_test$p.value
    
    return(data.frame(snv=rownames(snv_mtx)[i],
                      pheno=pheno_id,
                      n_people=nrow(om_dat),
                      base1=n_base1,
                      base2=n_base2,
                      pval=wil_p))
  }, error=function(e){
    message("SNV ", rownames(snv_mtx)[i], " ERROR: ", conditionMessage(e))
    return(data.frame(snv=rownames(snv_mtx)[i],
                      pheno=pheno_id,
                      n_people=NA,
                      base1=NA,
                      base2=NA,
                      pval=NA))
  })
}


n_cores <- detectCores() - 1
results_list <- mclapply(1:nrow(snv_mtx), process_snv, mc.cores = n_cores)

results <- bind_rows(results_list) %>% arrange(pval)


write.csv(results, paste0(species, '__', pheno_id, '.csv'), quote=F, row.names=F)


