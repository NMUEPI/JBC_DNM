options(stringsAsFactors = FALSE)

source("JBC_DNM_func.r")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggExtra)
library(rms)

# =========================================================
# Figure 2C. Restricted cubic spline analysis
# =========================================================

dd <- datadist(dnm_jbc_input)
options(datadist = "dd")

fit <- ols(dSNVM ~ rcs(age_mat,quantile(age_mat,c(0.5,0.6,0.75))) + age_pat, data = dnm_jbc_input)
pred <- Predict(fit,age_mat,age_pat = mean(dnm_jbc_input$age_pat,na.rm = TRUE),conf.int = 0.95)
pred <- as.data.frame(pred)

hist_info <- hist(dnm_jbc_input$age_mat, breaks = 20, plot = FALSE)
hist_df <- data.frame(age = hist_info$mids, count = hist_info$counts)

scale_factor <- max(hist_df$count) / max(pred$upper)
hist_df$scaled_count <- hist_df$count / scale_factor

knot1 <- quantile(dnm_jbc_input$age_mat, 0.50, na.rm = TRUE)
knot3 <- quantile(dnm_jbc_input$age_mat, 0.75, na.rm = TRUE)

ggplot() +
  geom_col(data = hist_df, aes(age, scaled_count), width = diff(hist_info$breaks)[1], fill = "#C7DDE8", alpha = 0.8) +
  geom_ribbon(data = pred, aes(x = age_mat, ymin = lower, ymax = upper), fill = "#F4A6A6", alpha = 0.35) +
  geom_line(data = pred, aes(x = age_mat, y = yhat), colour = "#E06A6A", linewidth = 0.8) +
  geom_vline(xintercept = c(knot1, knot3), linetype = 2, colour = "#3A6EA5", linewidth = 0.5) +
  annotate("text",x = 18,y = max(pred$upper)*0.95,label = "Pnonlinear = 1.05×10^-3",hjust = 0,size = 3.5) +
  scale_y_continuous(name = "Number of maternal dnSNVs",sec.axis = sec_axis( ~ . * scale_factor, name = "Number of samples")) +
  labs(x = "Maternal age at conception") +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10))