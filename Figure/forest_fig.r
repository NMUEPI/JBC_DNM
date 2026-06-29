# ====================================================================
# Forest plot generation
#
# Notes:
# 1. Association analyses are first performed using:
#      - generate_multi_linear()   (linear regression)
#      - generate_multi_logit()    (logistic regression)
#
# 2. These functions generate the regression summary tables (Beta/OR,
#    95% CI, P values), which are subsequently reformatted into the
#    plotting dataset (prepare_forest).
#
# 3. The forest plots shown in the manuscript are generated directly
#    from these formatted regression summary tables using the
#    forestploter package.
# ====================================================================

library(dplyr)
library(stringr)
library(ggplot2)
library(forestploter)
library(grid)

prepare_forest <- function(res){
    df <- as.data.frame(res)

    beta <- ifelse(df$BetaCI=="ref", NA, as.numeric(sub(" .*","",df$BetaCI)))
    low <- ifelse(df$BetaCI=="ref", NA, as.numeric(gsub(".*\\((.*), .*","\\1",df$BetaCI)))
    high <- ifelse(df$BetaCI=="ref", NA, as.numeric(gsub(".*, (.*)\\)","\\1",df$BetaCI)))

    p <- df$P
    p[p!="-"] <- format(as.numeric(p[p!="-"]), scientific = TRUE,digits = 2)

    ci <- ifelse(df$BetaCI=="ref", "Ref", sprintf("%.2f (%.2f, %.2f)", beta, low, high))

    plot_df <- data.frame(Factors=df$Factors, Levels=df$Levels, `Mean±sd`=df$`Mean±sd`, n=df$n, plot=paste(rep(" ",nrow(df)),collapse=""), BetaCI=ci, P=p, Beta=beta, Low=low, High=high, size=0.8)

    return(plot_df)
}

plot_forest <- function(plot_df, ci_col="#e1061f"){
    p <- forest(plot_df[,c("Factors", "n", "Mean±sd", "plot", "BetaCI", "P")],
		est=plot_df$Beta,
        lower=plot_df$Low,
        upper=plot_df$High,
        sizes=plot_df$size,
        ci_column=4,
        ref_line=0,
        theme=forest_theme(ci_pch=15)
	)

    p <- edit_plot(p, col=4, which="ci", gp=gpar(col=ci_col))
    return(p)
}
