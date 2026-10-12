library(vegan)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
species <- args[1]

set.seed(1234)

pheno_df <- read.delim('PRJNA834801_phenotype_num_filter.xls', check.names = F, row.names=1)

all_geo_matrix <- read.csv("PRJNA834801_geograph_dis.csv", header=TRUE, row.names=1, check.names=F)

result <- data.frame(sp=NA, pheno=NA, dis_sample=NA, N_pheno=NA, N_pair=NA, 
                     rho=NA, spear_p=NA, lm_beta=NA, lm_pval=NA, 
                     mantel_r=NA, mantel_p=NA, partial_mantel_r=NA, partial_mantel_p=NA)


for (pheno_id in names(pheno_df)) {
  #pheno_id <- 'Age'  
  dis_sample <- NA
  N_pheno <- NA
  N_pair <- NA
  rho <- NA
  spear_p <- NA
  beta_val <- NA
  p_val <- NA
  man_r <- NA
  man_p <- NA
  par_man_r <- NA
  par_man_p <- NA
    
  tryCatch({
    print(paste("Processing:", pheno_id))
    
    sp_dis <- read.delim(paste0(dir_path, 'snv_matrix/', species, '/snv_mtx2base_filter.mldist_pivotDis.xls.gz'), check.names = F, row.names = 1)


    sp_mag <- sub("\\.[^.]+$", "", rownames(sp_dis))
    unique_mag <- unique(sp_mag)
    N_unique_mag <- length(unique_mag)

    filtered_df <- sp_dis[!duplicated(sp_mag), !duplicated(sp_mag)]
    rownames(filtered_df) <- unique_mag
    colnames(filtered_df) <- unique_mag
    snv_matrix <- as.matrix(filtered_df)
    filtered_dis <- as.dist(filtered_df)

    sub_pheno <- pheno_df[row.names(pheno_df) %in% unique_mag, pheno_id, drop = FALSE]
    sub_pheno_df <- sub_pheno[match(unique_mag, row.names(sub_pheno)), , drop = FALSE]
    stopifnot(identical(rownames(sub_pheno_df), rownames(snv_matrix)))
    
    #sub_pheno_df <- na.omit(sub_pheno)
    N_pheno <- nrow(sub_pheno_df)

    if (N_pheno < 3) {
      print(paste("Insufficient samples:", N_pheno))
      next
    }
    
    pheno_mtr <- as.matrix(dist(sub_pheno_df[, 1, drop = FALSE], method = 'manhattan'))
    #pheno_mtr[is.na(pheno_mtr)] <- ""
    pheno_dis <- as.dist(pheno_mtr)
   

    geo_matrix <- all_geo_matrix[row.names(sub_pheno_df), row.names(sub_pheno_df)]
    geo_dis <- as.dist(geo_matrix)
      
   
    v_snv <- snv_matrix[upper.tri(snv_matrix)]
    v_pheno <- pheno_mtr[upper.tri(pheno_mtr)]
    v_geo <- geo_matrix[upper.tri(geo_matrix)]
      
    df <- data.frame(snv=v_snv, pheno=v_pheno, geo=v_geo)
    df <- na.omit(df)
    N_pair <- nrow(df)
      
    if (N_pair < 3) {
      print(paste("Insufficient pairs:", N_pair))
      next
    }
      

    scale_df <- as.data.frame(scale(df))
      

    spear_test <- cor.test(df$snv, df$pheno, method='spearman')
    rho <- spear_test$estimate
    spear_p <- spear_test$p.value
      

    lm_result <- summary(lm(snv ~ pheno + geo, data=scale_df))
    beta_val <- coef(lm_result)[2, 1]
    p_val <- coef(lm_result)[2, 4]
      

    man_test <- mantel(filtered_dis, pheno_dis, method='spearman', permutations=999, na.rm = TRUE)
    man_r <- man_test$statistic
    man_p <- man_test$signif
      
 
    par_man_test <- mantel.partial(filtered_dis, pheno_dis, geo_dis, method='spearman', permutations=999, na.rm = TRUE)
    par_man_r <- par_man_test$statistic
    par_man_p <- par_man_test$signif
      
  }, error = function(e) {
    print(paste("Error in", species, pheno_id, ":", e$message))
  })

  result <- rbind(result, data.frame(
    sp=species, 
    pheno=pheno_id, 
    dis_sample=N_unique_mag, 
    N_pheno=N_pheno, 
    N_pair=N_pair,
    rho=rho, 
    spear_p=spear_p, 
    lm_beta=beta_val, 
    lm_pval=p_val, 
    mantel_r=man_r, 
    mantel_p=man_p, 
    partial_mantel_r=par_man_r, 
    partial_mantel_p=par_man_p
  ))
}


result <- result[-1, ]
result <- result %>%  arrange(spear_p)
result$spear_fdr <- p.adjust(result$spear_p, method = 'BH')

write.csv(result, paste0(species, '_adj_geo_scale.csv'), quote=F, row.names=F)



