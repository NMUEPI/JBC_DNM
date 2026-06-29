options(stringsAsFactors = FALSE)

source("JBC_DNM_func.r")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggExtra)
library(rms)

# ============================================================
# Figure 2A. Linear association between parental age and DNMs
# ============================================================
## Prepare plotting data

figure2a_data <- rbind(
  data.frame(
    Age = dnm_jbc_input$age_mat,
    Num = dnm_jbc_input$DNMM,
    Pheno = "Maternal"),
  data.frame(
    Age = dnm_jbc_input$age_pat,
    Num = dnm_jbc_input$DNMP,
    Pheno = "Paternal")
)

## Scatter plot with linear regression

p1 <- ggplot(figure2a_data,aes(Age, Num, colour = Pheno)) +
  geom_point(size = 1) +
  stat_poly_line(formula = y ~ x, show.legend = FALSE) +
  stat_poly_eq(formula = y ~ x, aes(label = ..eq.label..), label.x.npc = "right", label.y.npc = "bottom", size = 4) +
  theme_classic() +
  labs(x = "Parental age at conception", y = "Standardized number of DNMs") +
  scale_colour_manual(values = c("#FBB255","#0066BD"))

ggMarginal(p1, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.75, bins = 30)