options(stringsAsFactors = FALSE)

source("JBC_DNM_func.r")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggExtra)
library(rms)
library(tidyverse)

# ============================================================
# Figure 3A. Association between ART and de novo mutations
# ============================================================

figure3a_data <- bind_rows(data.frame(
    Value = dnm_jbc_input$dSNV,
    Group = dnm_jbc_input$Conception,
    Type = "Overall dnSNV"
  ),
  data.frame(
    Value = dnm_jbc_input$dSNVP,
    Group = dnm_jbc_input$Conception,
    Type = "Paternal dnSNV"
  ),
  data.frame(
    Value = dnm_jbc_input$dSNVM,
    Group = dnm_jbc_input$Conception,
    Type = "Maternal dnSNV"
  ),
  data.frame(
    Value = dnm_jbc_input$dSNVmo,
    Group = dnm_jbc_input$Conception,
    Type = "EPZM"
  )
)

figure3a_data$Group <- factor(figure3a_data$Group, levels=c("ART","SP"), labels=c("ART","NC"))
figure3a_data$Type <- factor(figure3a_data$Type, levels=c("Overall dnSNV","Paternal dnSNV","Maternal dnSNV","EPZM"))

anno_df <- data.frame(
  Type = c("Overall dnSNV","Paternal dnSNV","Maternal dnSNV","EPZM"),
  label = c(
    "Beta = 2.38, P = 8.73×10^-20",
    "Beta = 0.74, P = 2.96×10^-3",
    "Beta = 0.62, P = 3.30×10^-4",
    "Beta = 1.02, P = 7.03×10^-34"
  )
)

ggplot(figure3a_data, aes(Group, Value, fill=Group)) +
    geom_violin(trim=FALSE, colour="black", alpha=0.8) +
    geom_boxplot(width=0.18, outlier.shape=NA, fill="white", colour="black") +
    facet_wrap(~Type, nrow=1, scales="free_y") +
    geom_text(data=anno_df, aes(x=1.5, y=Inf, label=label), inherit.aes=FALSE, vjust=1.4, size=3.2, fontface="bold") +
    scale_fill_manual(values=c("ART"="#DF6B63", "NC"="#6FA8CC")) +
    labs(x=NULL, y="Standardized number of dnSNVs") +
    theme_classic() +
    theme(
        strip.background=element_blank(),
        strip.text=element_text(size=11, face="bold"),
        axis.text=element_text(size=10, colour="black"),
        axis.title=element_text(size=11, face="bold"),
        legend.position="right",
        legend.title=element_blank()
    )
