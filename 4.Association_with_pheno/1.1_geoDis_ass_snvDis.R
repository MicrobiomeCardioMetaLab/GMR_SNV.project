library(vegan)
library(dplyr)

set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
species <- args[1]


pheno_df <- read.delim('PRJNA834801_phenotype_num_filter.xls', check.names = F, row.names=1)
all_geo_matrix <- read.csv("PRJNA834801_geograph_dis.csv", header=TRUE, row.names=1, check.names=F)
    
dis_matrix <- read.delim(paste0(species, '/snv_mtx2base_filter.mldist_pivotDis.xls.gz'), check.names = F, row.names = 1)
N_sample <- dim(dis_matrix)[1]


sp_mag <- sub("\\.[^.]+$", "", rownames(dis_matrix))
unique_mag <- unique(sp_mag)
N_unique_mag <- length(unique_mag)

filtered_df <- dis_matrix[!duplicated(sp_mag), !duplicated(sp_mag)]
#filtered_df <- dis_matrix[unique_mag, unique_mag]
print(dim(filtered_df))
rownames(filtered_df) <- unique_mag
colnames(filtered_df) <- unique_mag
snv_matrix <- as.matrix(filtered_df)
snv_dis <- as.dist(snv_matrix)    


geo_matrix <- all_geo_matrix[unique_mag, unique_mag]
print(dim(geo_matrix))
geo_dis <- as.dist(geo_matrix)
      
   
v_snv <- snv_matrix[upper.tri(snv_matrix)]
v_geo <- geo_matrix[upper.tri(geo_matrix)]
      
df <- data.frame(snv=v_snv, geo=v_geo)
df <- na.omit(df)
N_pair <- nrow(df)
      

spear_test <- cor.test(df$snv, df$geo, method='spearman')
rho <- spear_test$estimate
spear_p <- spear_test$p.value
      

man_test <- mantel(snv_dis, geo_dis, method='spearman', permutations=999, na.rm = TRUE)
man_r <- man_test$statistic
man_p <- man_test$signif
      
 



results <- data.frame(
  sp = species,
  n_sample = N_sample,
  n_pair = N_pair,
  rho = rho,
  spear_pval = spear_p,
  man_r = man_r,
  man_p = man_p
)

write.csv(results, paste0(species, 'geoDis_ass_snvDis.csv'), quote=F, row.names=F)




