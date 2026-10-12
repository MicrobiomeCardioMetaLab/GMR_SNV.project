library(ggplot2)
library(dplyr)
library(tidyr) 

my_theme <- function(){
  theme_bw() +
    theme(legend.title=element_blank(), legend.text = element_text(size = 5)) +
    theme(plot.title = element_text(hjust=0.5,vjust=0.12, size = 6)) +
    theme(axis.title.x = element_text(size = 6,vjust = 0.8)) +
    theme(axis.title.y = element_text(size = 6,vjust = 0.8)) +
    theme(axis.text.x = element_blank())+
    theme(axis.text.y = element_text(size = 5,hjust = 1, color = "#666666"))+
    theme(axis.ticks.y = element_line(colour = "#666666")) +
    theme(axis.ticks.x = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank(),axis.line = element_line(size = 0.3, colour = "black"))
}


sp_list <- read.csv('base2_2637treeID.txt', header=F)

result <- data.frame(sp=NA, heap_a=NA, heap_b=NA, heap_pval=NA, new_SNV=NA, new_frac=NA, pred_N=NA, N_MAG=NA)

for (sp_name in sp_list[, 1]) {
  tryCatch({

    df <- read.delim(paste0(sp_name, '.xls'))
    df <- df[-1, ]
    df[1:5, 1:5]
    step <- df[, 1]
    snv_mat <- df[, -1]


    mean_snv <- apply(snv_mat, 1, mean)
    sd_snv   <- apply(snv_mat, 1 ,sd)

    summary_df <- data.frame(
      step = step,
      mean_snv = mean_snv,
      sd_snv = sd_snv
    )


    summary_df$slope <- c(NA, diff(summary_df$mean_snv))
    summary_df$new_fraction <- summary_df$slope / summary_df$mean_snv



    scale_factor <- max(summary_df$mean_snv) / max(summary_df$slope, na.rm=TRUE)
    p2 <- ggplot(summary_df, aes(x = step)) +
            geom_line(aes(y = mean_snv), color = "#0076B9", size = 0.9) +
            geom_ribbon(aes(ymin = mean_snv - sd_snv, ymax = mean_snv + sd_snv), fill = "#ABD9E9", alpha = 0.3) +
            geom_line(aes(y = slope * scale_factor), color = "#cc5b45", size = 0.8) +
            scale_y_continuous(name = "Total SNVs", sec.axis = sec_axis(~./scale_factor, name = "New SNVs per sample")) +
            labs(x = "Sample size", title = sp) +
            my_theme()



    fit <- lm(log(mean_snv) ~ log(step), data=summary_df)

    coef_fit <- coef(fit)

    a <- exp(coef_fit[1])
    b <- coef_fit[2]
    pval <- summary(fit)$coefficients[2, 4]
    cat("Heaps law model:\n")
    cat("S(N) =",a,"* N^",b,"\n")

    slope_fun <- function(N){
      a*b*N^(b-1)
    }

    N_pred <- seq(max(step),500000,1)

    pred_slope <- slope_fun(N_pred)

    plateau <- N_pred[pred_slope < 100][1]

    cat("Estimated plateau sample size:", plateau,"\n")

    result <- rbind(result, c(sp_name, a, b, pval, summary_df[nrow(summary_df), ncol(summary_df) -1], summary_df[nrow(summary_df), ncol(summary_df)], plateau, nrow(summary_df) + 1))
   },  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

result <- result[-1, ]
write.csv(result, 'SNV_cumul_analysis.csv', quote=F, row.names=F)





