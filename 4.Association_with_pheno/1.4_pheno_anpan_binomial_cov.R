library(anpan)
library(ape)
library(loo)
library(tidyverse)


set.seed(1234)


args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
pheno <- args[2]



pheno_df <- read.delim('PRJNA834801_phenotype_num_filter.xls', check.names = F)
pheno_df <- pheno_df[, c('sample_id', pheno, 'location'), drop = FALSE]

pheno_df <- pheno_df[complete.cases(pheno_df[[pheno]]), ]

pheno_df[[pheno]] <- factor(pheno_df[[pheno]],
                            levels = c(1, 2),
                            labels = c("control", "case"))
valid_samples <- pheno_df %>%
  pull(sample_id) %>%
  unique()
print(valid_samples)

cat("Processing:", species, "\n")
 
raw_tr <- read.tree(paste0(species, '/snv_mtx2base_filter.treefile'))
raw_tr$tip.label <- sub("\\.[^.]+$", "", raw_tr$tip.label)
all_tips <- raw_tr$tip.label

keep_tips <- intersect(valid_samples, all_tips)

 
tr <- keep.tip(raw_tr, keep_tips)

metadata <- pheno_df[pheno_df$sample_id %in% keep_tips, , drop = FALSE]


metadata <- metadata[match(tr$tip.label, metadata$sample_id), ,drop = FALSE]

case_df <- table(metadata[[pheno]])


n_control <- ifelse("control" %in% names(case_df), case_df[["control"]], 0)

n_case <- ifelse("case" %in% names(case_df), case_df[["case"]], 0)

 

anpan_result <- anpan_pglmm(
      meta_file = metadata,
      tree_file = tr,
      outcome = pheno,
      covariates = "location",
      family = "binomial",
      bug_name = species,
      reg_noise = TRUE,
      omit_na = TRUE,
      loo_comparison = TRUE,
      show_plot_tree = FALSE,
      show_post  = FALSE
    )



base_k_values <- pareto_k_values(
      anpan_result$loo$base_loo
    )

pglmm_k_values <- pareto_k_values(
      anpan_result$loo$pglmm_loo
    )

comp_df <- as.data.frame(
      anpan_result$loo$comparison
    )

row <- comp_df[
      comp_df$model == "base_fit",
      ,
      drop = FALSE
    ]



results <- data.frame(
      sp = species,
      pheno_id = pheno,

      n_control = n_control,
      n_case = n_case,
      n_people = length(tr$tip.label),

      elp_diff = row$elpd_diff,
      se_diff = row$se_diff,
      p_worse = row$p_worse,

      n_base_k = length(base_k_values),
      max_base_k = max(base_k_values, na.rm = TRUE),
      pct_base_k_bad = mean(base_k_values >= 0.7, na.rm = TRUE) * 100,

      n_pglmm_k = length(pglmm_k_values),
      max_pglmm_k = max(pglmm_k_values, na.rm = TRUE),
      pct_pglmm_k_bad = mean(pglmm_k_values >= 0.7, na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )


write.csv(results, paste0(dir_path, species, '_anpan.csv'), quote=F, row.names=F)

