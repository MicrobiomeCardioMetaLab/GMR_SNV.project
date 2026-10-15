library(vegan)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
species <- args[1]


pheno_df <- read.delim('../data/phenotype.xls', check.names = F, row.names=1)


results <- data.frame(sp=NA, pheno=NA, N_dis=NA, n_dis_pheno=NA, r2=NA, pval=NA)

for (pheno_id in c('gender', 'BMI', 'Case_status', 'Age')) {
    
  tryCatch({
    
    print(paste("Processing:", pheno_id))
    
    sp_dis <- read.delim(paste0('../Example/', species, '/snv_matrix/snv_mtx2base_filter.mldist_pivotDis.xls'), check.names = F, row.names = 1)
    sp_mag <- sub("\\.[^.]+$", "", rownames(sp_dis))

    unique_mag <- unique(sp_mag)
    N_sample <- length(unique_mag)

    filtered_df <- sp_dis[!duplicated(sp_mag), !duplicated(sp_mag)]
    rownames(filtered_df) <- unique_mag
    colnames(filtered_df) <- unique_mag
    snv_matrix <- as.matrix(filtered_df)


    sub_pheno <- pheno_df[row.names(pheno_df) %in% unique_mag, pheno_id, drop = FALSE]
    sub_pheno_df <- sub_pheno[match(unique_mag, row.names(sub_pheno)), , drop = FALSE]


    
    sub_pheno_df <- na.omit(sub_pheno)
    names(sub_pheno_df) <- 'pheno_name'
    
    if (dim(sub_pheno_df)[1] < 3) {  
      results <- rbind(results, c(species, pheno_id, N_sample, dim(sub_pheno_df)[1], NA, NA))
      print(paste("Insufficient overlap samples:", dim(sub_pheno_df)[1], "skipping correlation analysis"))
      next  
    }
    
    dis_df <- snv_matrix[row.names(sub_pheno_df), row.names(sub_pheno_df)]

    
    test <- adonis2(dis_df ~ pheno_name, data = sub_pheno_df, permutations = 999, by = "margin")
    
    R2 <- test$R2[1]
    pvalue <- test$`Pr(>F)`[1] 

    results <- rbind(results, c(species, pheno_id, N_sample, dim(sub_pheno_df)[1], R2, pvalue))
    
  }, error = function(e) {
    warning(paste("Error processing", pheno_id, ":", e$message))
    NULL
  })
}


results <- results[-1, ]
results$r2 <- ifelse(results$r2 > 0, results$r2, 0)
results <- results %>%  arrange(pval)
results$sg_sp_fdr <- p.adjust(results$pval, method = 'BH')

dir.create("../results/association/", recursive = TRUE, showWarnings = FALSE)
write.csv(results, paste0('../results/association/', species, '_adnois.csv'), quote=F, row.names=F)


