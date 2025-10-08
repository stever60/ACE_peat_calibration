# ACE PLSR log_inc using pls package 

# XRF log_inc model and Ti_ICP response variable to log_inc scaled and centered dataset

#                                                   ACE     HER42PB
# Subsamples measured by ICPMS (matched training)   265     66
# No. of measurements by XRF core scanning          14513   4009          
  
# Set up ------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()

# Libraries ---------------------------------------------------------------
packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'readr', 'ggpubr', 'ggsci', 
              'compositions', 'boot', 'broom',
              'wesanderson', 'viridis', 'itraxR','pls', 'caret')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

#-------------------------------------------------------------------------------
# Define parameters ------------------------------------------------------------

# Key XRF and ICPMS elements
xrf_icp_elements <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                      "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                      "Zr", "Zr_ICP", "Mo_coh")

xrf_icp_elements1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                       "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP", "Mo_inc",  "Mo_coh", "inc_coh", "coh_inc", "Dry_mass") # Mo_inc included

xrf_icp_elements2 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                       "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP")

# Load data --------------------------------------------------------------------

# Create a tibble of main 12 elements as predictors (ITRAX) and ICPMS equivalents as response (ICP-MS)
# Add predicted output to tibble at end

# Import the matched ICPMS-XRF log (element/inc) dataset and create a PLS-model for prediction of Ti_ICP
# Using up to 11 other elements as predictors. 

# Assumptions:
# Ti_ICPMS is the response variable.
# Other elements (e.g., Ca, Fe, Sr, Mn, etc.) measured by ITRAX are predictor variables.
# All data are numeric and cleaned (no missing values).
# The response variable is continuous.

# Load and Prepare data
ACE_matched_xrf_icp_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_xrf_icp_matched_log_inc.csv") 
is.na(ACE_matched_xrf_icp_log_inc)<-sapply(ACE_matched_xrf_icp_log_inc, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_log_inc

# Standardised and centred dataframe - Z-scores centred around 0
ACE_matched_xrf_icp_log_inc_Z <- ACE_matched_xrf_icp_log_inc 
ACE_matched_xrf_icp_log_inc_Z[, xrf_icp_elements] <- scale(ACE_matched_xrf_icp_log_inc[, xrf_icp_elements], center = TRUE, scale = TRUE)
ACE_matched_xrf_icp_log_inc_Z

# Select dataset and convert to base R dataframe for pls package ---------------
df0 <- ACE_matched_xrf_icp_log_inc %>%
  select(c(Site, depth, SH20_mean_age, SH20_mean_95CI, all_of(xrf_icp_elements1))) %>%
  drop_na()
df0

# Convert tibble to data frame for base R
df <- as.data.frame(df0)

#-------------------------------------------------------------------------------

# 1. Partial Least Squares (PLS) regression validation & optimal components-----

# ------------------------------------------------------------------------------

# Run once only
# Select dataset and convert to base R dataframe for pls package ---------------
# If the dataset is small, use eg k-fold cross-validation to evaluate the model's performance without relying on a single test set. 
# Calibrate Titanium (Ti) measured by ICP-MS against other elements measured by ITRAX using pls package.

# Define predictor ITRAX variables & ICPMS response variable for calibration model tests
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
main_elements_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP","Fe_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP","Zr_ICP")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")

# Split into predictors (X) using key elements only and response (Y) as ICPMS element e.g., Ti_ICP
X_test <- df[, key_elements_xrf]  # ITRAX elements
Y_test <- df$Ti_ICP #define response ICP-MS variable for calibration

# Write summary outputs from console to txt file in Output folder instead of console
# sink(file = "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Summary_stats.txt")
# Validation model - assessing signifcant elements & number of com --------
# LOO-PLSR (leave-one-out) test model to determine optimal number of comps -----

# Perform PLS regression with LOO - leave one out validation - useful for unbiased prediction error assessment
pls_train_model_LOO <- plsr(Y_test ~ ., data = data.frame(X_test, Y_test), validation = "LOO", jackknife = TRUE, scale = TRUE, center = TRUE) # Use scale = TRUE to perform as Z scores - do not use for calibration
RMSEP(pls_train_model_LOO)

# Get optimal number of components (lowest LOO error)
opt_comp_LOO <- which.min(RMSEP(pls_train_model_LOO)$val[1, , -1])
cat("Optimal number of components LOO:", opt_comp_LOO, "\n")

# Summary of LOO model output and jackknife t tests of regression coefficients & std errors
summary(pls_train_model_LOO)
jack.test(pls_train_model_LOO, ncomp = opt_comp_LOO, use.mean = TRUE)

# CV-PLSR (cross-validated) test model to determine optimal number of comps ----

# Perform PLS regression with 10-fold cross-validation for model tuning - useful for small datasets
pls_train_model_CV <- plsr(Y_test ~ ., data = data.frame(X_test, Y_test), validation = "CV", jackknife = TRUE, scale = TRUE, center = TRUE) # Use scale = TRUE to perform as Z scores - do not use for calibration

# Get optimal number of components (lowest CV error)
opt_comp_CV <- which.min(RMSEP(pls_train_model_CV)$val[1, , -1])
cat("Optimal number of components CV:", opt_comp_CV, "\n")

# Summary of  CV model output and jackknife t tests of regression coefficients & std errors
summary(pls_train_model_CV)
jack.test(pls_train_model_CV, ncomp = opt_comp_CV, use.mean = TRUE)

# Plot root mean squared error of prediction (RMSEP) to choose number of components
plot(RMSEP(pls_train_model_CV), legendpos = "topright")

# ------------------------------------------------------------------------------

# 1.2 PSLR training model with Ti, Ca, Sr, Zr as element matrix ----------------

# Only use signifcant predictor ITRAX variables from LOO & CV jackknife results for final model
final_elements_xrf <- c("Ti", "Ca", "Sr", "Zr")
final_elements_icp <- c("Ti_ICP", "Ca_ICP", "Sr_ICP", "Zr_ICP")
opt_comp_final = 4

df_xrf <- df[, final_elements_xrf]
df_icp <- df[, final_elements_icp]

df_xrf_matrix <- as.matrix(df[, final_elements_xrf])
df_icp_matrix <- as.matrix(df[, final_elements_icp])

# Refit training model with optimal number of components
pls_final <- plsr(df_icp_matrix ~ df_xrf_matrix, ncomp = opt_comp_final, validation = "CV", jackknife = TRUE, scale = TRUE, center = TRUE)
RMSEP(pls_train_model_LOO)

# Summary of final plsr model and jackknife t tests of regression coefficients & std errors
summary(pls_final)
jack.test(pls_final, ncomp = opt_comp_final, use.mean = TRUE)

# Predict Ti using the model
predicted_icp <- predict(pls_final, ncomp = opt_comp_final)
colnames(predicted_icp) <- paste0(colnames(predicted_icp), "_pred")
pred_icp_tbl <- as_tibble(predicted_icp) %>% 
  rename(Ti_pred = 1, Ca_pred = 2, Sr_pred = 3, Zr_pred = 4)
actuak_icp <- df_icp
actual_icp_tbl <- as_tibble(df_icp) %>% 
  rename(Ti_act = 1, Ca_act = 2, Sr_act = 3, Zr_act = 4)

# Plot predicted vs actual draft output
plot(actual_icp_tbl, predicted_icp)

# R2, RMSE predictions for final model base R ----------------------------------

# Predict Ti using the model
predicted_Ti0 <- predict(pls_final, ncomp = opt_comp_final)
Ti_predicted <- predicted_Ti0[, "Ti_ICP", 1]
Ti_actual <- df_icp$Ti_ICP
length(Ti_predicted) == length(Ti_actual)  # should be TRUE
# R-squared
r2_Ti <- cor(Ti_actual, Ti_predicted)^2
# rmse_Ti
rmse_Ti <- sqrt(mean((Ti_actual - Ti_predicted)^2))
# RMSEP (same as rmse_Ti if on test set)
rmsep_Ti <- sqrt(mean((Ti_actual - Ti_predicted)^2))
cat("Ti_R² =", r2_Ti, "\nTi_RMSE =", rmse_Ti, "\nTi_RMSEP =", rmse_Ti)

# Predict Ca using the model
predicted_Ca0 <- predict(pls_final, ncomp = opt_comp_final)
Ca_predicted <- predicted_Ca0[, "Ca_ICP", 1]
Ca_actual <- df_icp$Ca_ICP
length(Ca_predicted) == length(Ca_actual)  # should be TRUE
# R-squared
r2_Ca <- cor(Ca_actual, Ca_predicted)^2
# rmse_Ti
rmse_Ca <- sqrt(mean((Ca_actual - Ca_predicted)^2))
# RMSEP (same as rmse_Ti if on test set)
rmsep_Ca <- sqrt(mean((Ca_actual - Ca_predicted)^2))
cat("Ca_R² =", r2_Ca, "\nCa_RMSE =", rmse_Ca, "\nCa_RMSEP =", rmsep_Ca)

# Predict Sr using the model
predicted_Sr0 <- predict(pls_final, ncomp = opt_comp_final)
Sr_predicted <- predicted_Sr0[, "Sr_ICP", 1]
Sr_actual <- df_icp$Sr_ICP
length(Sr_predicted) == length(Sr_actual)  # should be TRUE
# R-squared
r2_Sr <- cor(Sr_actual, Sr_predicted)^2
# rmse_Ti
rmse_Sr <- sqrt(mean((Sr_actual - Sr_predicted)^2))
# RMSEP (same as rmse_Ti if on test set)
rmsep_Sr <- sqrt(mean((Sr_actual - Sr_predicted)^2))
cat("Sr_R² =", r2_Sr, "\nSr_RMSE =", rmse_Sr, "\nSr_RMSEP =", rmsep_Sr)

# Predict Zr using the model
predicted_Zr0 <- predict(pls_final, ncomp = opt_comp_final)
Zr_predicted <- predicted_Zr0[, "Zr_ICP", 1]
Zr_actual <- df_icp$Zr_ICP
length(Zr_predicted) == length(Zr_actual)  # should be TRUE
# R-squared
r2_Zr <- cor(Zr_actual, Zr_predicted)^2
# rmse_Ti
rmse_Zr <- sqrt(mean((Zr_actual - Zr_predicted)^2))
# RMSEP (same as rmse_Ti if on test set)
rmsep_Zr <- sqrt(mean((Zr_actual - Zr_predicted)^2))
cat("Zr_R² =", r2_Zr, "\nZr_RMSE =", rmse_Zr, "\nZr_RMSEP =", rmsep_Zr)

# Coefficients of the final model ---------------------------------------
coefficients(pls_final)

# Get coefficients and intercept from the model
coef_matrix <- coef(pls_final, ncomp = opt_comp_final)
#intercept <- attr(coef_matrix, "constant")

# Extract regression coefficients for optimal components
coefs <- coef(pls_final, ncomp = opt_comp_final)
coefs_vec <- as.vector(coefs)
names(coefs_vec) <- dimnames(coefs)[[1]]

# Show regression equation
coefs <- coef(pls_final, ncomp = opt_comp_final)[,,1]

# Ti
coefs_ti <- coefs[, "Ti_ICP", drop = FALSE]
cat("Ti[ICP-MS_pred] = ",  
    paste0(round(coefs_ti[,1], 4), "*", rownames(coefs_ti), collapse = " + "), "\n")
# Ca
coefs_Ca <- coefs[, "Ca_ICP", drop = FALSE]
cat("Ca[ICP-MS_pred] = ",  
    paste0(round(coefs_Ca[,1], 4), "*", rownames(coefs_Ca), collapse = " + "), "\n")
# Sr
coefs_Sr <- coefs[, "Sr_ICP", drop = FALSE]
cat("Sr[ICP-MS_pred] = ",  
    paste0(round(coefs_Sr[,1], 4), "*", rownames(coefs_Sr), collapse = " + "), "\n")
# Zr
coefs_Zr <- coefs[, "Zr_ICP", drop = FALSE]
cat("Zr[ICP-MS_pred] = ",  
    paste0(round(coefs_Zr[,1], 4), "*", rownames(coefs_Zr), collapse = " + "), "\n")

# Generate R2, RMSE, equation text for ggplot for all variables  ---------------

# Ti
R2_text_Ti <- paste0("Ti_logR^2: ", paste0(round(r2_Ti, 2)))
RMSE_text_Ti <- paste0("Ti_logRMSE: ", paste0(round(rmse_Ti, 4)))
eq_text_Ti <- paste0("Ti_ICPMS[pred_log]:", paste0(round(coefs_ti[,1], 4), "*", rownames(coefs_ti), collapse = " + "))
RMSEP_Ti <- "1.2634" 
RMSEP_text_Ti <- paste0("Ti_logRMSEP[bsmean]: ", paste0(RMSEP_Ti))# boostrapped RMSEP from 0.6 training:test model in Sect 2 below - same elements 
# Ca
R2_text_Ca <- paste0("Ca_logR^2: ", paste0(round(r2_Ca, 2)))
RMSE_text_Ca <- paste0("Ca_logRMSE: ", paste0(round(rmse_Ca, 4)))
eq_text_Ca <- paste0("Ca_ICPMS[pred_log]:", paste0(round(coefs_Ca[,1], 4), "*", rownames(coefs_Ca), collapse = " + "))
RMSEP_Ca <- "TBD" 
RMSEP_text_Ca <- paste0("Ca_logRMSEP[bsmean]: ", paste0(RMSEP_Ca))# boostrapped RMSEP from 0.6 training:test model in Sect 2 below - same elements 
# Sr
R2_text_Sr <- paste0("Sr_logR^2: ", paste0(round(r2_Sr, 2)))
RMSE_text_Sr <- paste0("Sr_logRMSE: ", paste0(round(rmse_Sr, 4)))
eq_text_Sr <- paste0("Sr_ICPMS[pred_log]:", paste0(round(coefs_Sr[,1], 4), "*", rownames(coefs_Sr), collapse = " + "))
RMSEP_Sr <- "TBD" 
RMSEP_text_Sr <- paste0("Sr_logRMSEP[bsmean]: ", paste0(RMSEP_Sr))# boostrapped RMSEP from 0.6 training:test model in Sect 2 below - same elements 
# Zr
R2_text_Zr <- paste0("Zr_logR^2: ", paste0(round(r2_Zr, 2)))
RMSE_text_Zr <- paste0("Zr_logRMSE: ", paste0(round(rmse_Zr, 4)))
eq_text_Zr <- paste0("Zr_ICPMS[pred_log]:", paste0(round(coefs_Zr[,1], 4), "*", rownames(coefs_Zr), collapse = " + "))
RMSEP_Zr <- "TBD" 
RMSEP_text_Zr <- paste0("Zr_logRMSEP[bsmean]: ", paste0(RMSEP_Zr))# boostrapped RMSEP from 0.6 training:test model in Sect 2 below - same elements 

# Variable Importance in Projection (VIP) scores summarize each  ---------------
# 
# # Function to calculate VIP scores
# vip_scores <- function(pls_final, opt_comp_final) {
#   W <- pls_final$loading.weights[, 1:opt_comp_final]
#   T <- pls_final$scores[, 1:opt_comp_final]
#   Q <- pls_final$Yloadings[1:opt_comp_final]
#   
#   SS <- colSums(T^2) * Q^2
#   total_SS <- sum(SS)
#   
#   vip <- sqrt(ncol(W) * rowSums((W^2) * SS) / total_SS)
#   names(vip) <- rownames(W)
#   return(vip)
# }
# 
# # Get VIP values
# vip <- vip_scores(pls_final, opt_comp_final)
# print(round(vip, 3))
# # VIP > 1: Important — strong contribution to the model
# # 0.8 < VIP < 1: Moderate contribution
# # VIP < 0.8: Unimportant — likely not useful for prediction
# 
# # Write console summary info to file 
# #sink(file = NULL)



# Convert to Tidyverse for calculations & plotting -----------------------------

# Create a tibble with site name and predicted output for ggplot
df_joined <- cbind(df, predicted_icp)
df_predicted <- as_tibble(df_joined) %>% 
  rename(Ti_ICP_pred = `Ti_ICP_pred.4 comps`) %>% 
  rename(Ca_ICP_pred = `Ca_ICP_pred.4 comps`) %>% 
  rename(Sr_ICP_pred = `Sr_ICP_pred.4 comps`) %>% 
  rename(Zr_ICP_pred = `Zr_ICP_pred.4 comps`)
df_predicted
write.csv(df_predicted,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Ti_ICP_predicted_log_inc.csv", row.names = FALSE)

# RMSE and R² and CI as log summary table for each element in Tidyverse --------

summary_df_log_inc <- map2_dfr(actual_icp_tbl, pred_icp_tbl, ~ {
  tibble(
    R2 = cor(.x, .y)^2,
    RMSE = sqrt(mean((.x - .y)^2))
  )
}, .id = "Element")

# Rename based on actual column names
summary_df_log_inc$Element <- names(actual_icp_tbl)
print(summary_df_log_inc)

ci_list <- map2(actual_icp_tbl, pred_icp_tbl, ~ {
  residuals <- .x - .y
  rmse_val <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  tibble(
    CI_lower = quantile(residuals, 0.05, na.rm = TRUE),
    CI_upper = quantile(residuals, 0.95, na.rm = TRUE)
  )
})

# RMSEP - should be the same as RMSE for training dataset
# Combine both into one tibble for comparison
comparison_df <- bind_cols(actual_icp_tbl, pred_icp_tbl)

# Calculate RMSEP for each element
rmsep_df <- comparison_df %>%
  transmute(
    Ti = (Ti_pred - Ti_act)^2,
    Car = (Ca_pred - Ca_act)^2,
    Sr = (Sr_pred - Sr_act)^2,
    Zr = (Zr_pred - Zr_act)^2
  ) %>%
  summarise(across(everything(), ~ sqrt(mean(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = "Element", values_to = "RMSEP") %>% 
  select(-Element)

# Combine into summary table
library(glue) # to make next labels for plotting
ci_df <- bind_rows(ci_list)
summary_df_log_inc_CI <- bind_cols(summary_df_log_inc, rmsep_df, ci_df) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  mutate(Element = str_replace(Element, "_act", "")) %>%
  mutate(label_R2 = paste0("R^2:", round(R2, 2))) %>%  
  mutate(label_RMSE = paste0("RMSE: ", round(RMSE, 3)))
print(summary_df_log_inc_CI)
write_csv(summary_df_log_inc_CI, "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Summary_log_inc_R2_RMSE_CI.csv")

# Plot Ti pred vs actual as log data, coloured by Site in ggplot ---------------

# Make labels
R2_text_Ti <- summary_df_log_inc_CI %>%
  filter(Element == "Ti") %>% 
  pull(label_R2)
RMSE_text_Ti <- summary_df_log_inc_CI %>%
  filter(Element == "Ti") %>% 
  pull(label_RMSE)

# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ti <- ggplot(df_predicted, aes(x = Ti_ICP, y = Ti_ICP_pred, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Ti ICP-MS log_inc model",
    x = "Measured Ln (Ti) [ICP-MS]",
    y = "Predicted Ln (Ti) [ICP-MS] (PLS log_inc Model)"
  ) +
  #theme_minimal() +
  coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 4.5, y = 9.25, label = RMSEP_text_Ti, parse = TRUE, hjust = 0) + 
  annotate("text", x = 4.5, y = 9.5, label = RMSE_text_Ti, parse = TRUE, hjust = 0) + 
  annotate("text", x = 4.5, y = 9.75, label = eq_text_Ti, parse = TRUE, hjust = 0) + 
  annotate("text", x = 4.5, y = 10, label = R2_text_Ti, parse = TRUE, hjust = 0)
print(PLS_final_Ti)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Ti_plsr_log_inc_as_log.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Ca pred vs actual as log data, colored by Site in ggplot ----------------

# Make labels
R2_text_Ca <- summary_df_log_inc_CI %>%
  filter(Element == "Ca") %>% 
  pull(label_R2)
RMSE_text_Ca <- summary_df_log_inc_CI %>%
  filter(Element == "Ca") %>% 
  pull(label_RMSE)

# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ca <- ggplot(df_predicted, aes(x = Ca_ICP, y = Ca_ICP_pred, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Ca ICP-MS log_inc model",
    x = "Measured Ln (Ca) ICP-MS",
    y = "Predicted Ln (Ca) ICP-MS (PLS log_inc Model)"
  ) +
  #theme_minimal() +
  coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  #annotate("text", x = 6, y = 10.25, label = RMSEP_text_Ca, parse = TRUE, hjust = 0) + 
  annotate("text", x = 6, y = 10.5, label = RMSE_text_Ca, parse = TRUE, hjust = 0) + 
  annotate("text", x = 6, y = 10.75, label = eq_text_Ca, parse = TRUE, hjust = 0) + 
  annotate("text", x = 6, y = 11, label = R2_text_Ca, parse = TRUE, hjust = 0)
print(PLS_final_Ca)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Ca_plsr_log_inc_as_log.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Sr pred vs actual as log data, colored by Site in ggplot ----------------

# Make labels
R2_text_Sr <- summary_df_log_inc_CI %>%
  filter(Element == "Sr") %>% 
  pull(label_R2)
RMSE_text_Sr <- summary_df_log_inc_CI %>%
  filter(Element == "Sr") %>% 
  pull(label_RMSE)

# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchiSrgo"
PLS_final_Sr <- ggplot(df_predicted, aes(x = Sr_ICP, y = Sr_ICP_pred, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Sr ICP-MS log_inc model",
    x = "Measured Ln (Sr) ICP-MS",
    y = "Predicted Ln (Sr) ICP-MS (PLS log_inc Model)"
  ) +
  coord_equal(xlim = c(3, 8), ylim = c(3, 8)) +
  #theme_minimal() +
  #coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 3, y = 7.25, label = RMSEP_text_Sr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 3, y = 7.5, label = RMSE_text_Sr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 3, y = 7.75, label = eq_text_Sr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 3, y = 8, label = R2_text_Sr, parse = TRUE, hjust = 0)
print(PLS_final_Sr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Sr_plsr_log_inc_as_log.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Zr pred vs actual as log data, colored by Site in ggplot ----------------

# Make labels
R2_text_Zr <- summary_df_log_inc_CI %>%
  filter(Element == "Zr") %>% 
  pull(label_R2)
RMSE_text_Zr <- summary_df_log_inc_CI %>%
  filter(Element == "Zr") %>% 
  pull(label_RMSE)

# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchiZrgo"
PLS_final_Zr <- ggplot(df_predicted, aes(x = Zr_ICP, y = Zr_ICP_pred, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Zr ICP-MS log_inc model",
    x = "Measured Ln (Zr) ICP-MS",
    y = "Predicted Ln (Zr) ICP-MS (PLS log_inc Model)"
  ) +
  coord_equal(xlim = c(0, 7), ylim = c(0, 7)) +
  #theme_minimal() +
  #coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 0, y = 6.25, label = RMSEP_text_Zr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 6.5, label = RMSE_text_Zr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 6.75, label = eq_text_Zr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 7, label = R2_text_Zr, parse = TRUE, hjust = 0)
print(PLS_final_Zr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Zr_plsr_log_inc_as_log.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")



# Convert to ppm & plot as measured vs actual in ppm  --------------------------

# Original data was log-transformed - therefore take exp of both sides to convert to ppm for ICP_pred
df_predicted_exp <- df_predicted %>% 
  select(Ti_ICP_pred:Zr_ICP_pred) %>% 
  exp() %>% 
  rename_with(~ paste0(., "_ppm"))
df_predicted_exp
write.csv(df_predicted,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Ti_ICP_predicted_ppm.csv", row.names = FALSE)

# Load original ICP ppm data for ICP_actual
ACE_matched_xrf_icp_cps <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") 
is.na(ACE_matched_xrf_icp_cps)<-sapply(ACE_matched_xrf_icp_cps, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_cps
df_predicted_ppm <- bind_cols(ACE_matched_xrf_icp_cps, df_predicted_exp)
write.csv(df_predicted_ppm,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/ACE_ICP_act_pred_ppm.csv", row.names = FALSE)

# RMSE and R² and CI summary table for each element  --------------------

pred_ppm <- df_predicted_ppm %>% 
  select(Ti_ICP_pred_ppm, Ca_ICP_pred_ppm, Sr_ICP_pred_ppm, Zr_ICP_pred_ppm)
actual_ppm <- df_predicted_ppm %>% 
  select(Ti_ICP, Ca_ICP, Sr_ICP, Zr_ICP)
  
summary_df_log_inc_ppm <- map2_dfr(actual_ppm, pred_ppm, ~ {
  tibble(
    R2 = cor(.x, .y)^2,
    RMSE = sqrt(mean((.x - .y)^2))
  )
}, .id = "Element")

# Rename based on actual column names
summary_df_log_inc_ppm$Element <- names(actual_ppm)
print(summary_df_log_inc_ppm)

ci_list_ppm <- map2(actual_ppm, pred_ppm, ~ {
  residuals <- .x - .y
  rmse_val <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  tibble(
    CI_lower = quantile(residuals, 0.05, na.rm = TRUE),
    CI_upper = quantile(residuals, 0.95, na.rm = TRUE)
  )
})

# RMSEP - should be the same as RMSE for training dataset
# Combine both into one tibble for comparison
comparison_df <- bind_cols(actual_ppm, pred_ppm)

# Calculate RMSEP for each element
rmsep_df_ppm <- comparison_df %>%
  transmute(
    Ti = (Ti_ICP_pred_ppm - Ti_ICP)^2,
    Car = (Ca_ICP_pred_ppm - Ca_ICP)^2,
    Sr = (Sr_ICP_pred_ppm - Sr_ICP)^2,
    Zr = (Zr_ICP_pred_ppm - Zr_ICP)^2
  ) %>%
  summarise(across(everything(), ~ sqrt(mean(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = "Element", values_to = "RMSEP") %>% 
  select(-Element)

# Combine into summary table
library(glue) # to make next labels for plotting
ci_df_ppm <- bind_rows(ci_list_ppm)
summary_df_log_inc_CI_ppm <- bind_cols(summary_df_log_inc_ppm, rmsep_df_ppm, ci_df_ppm) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  mutate(Element = str_replace(Element, "_act", "")) %>%
  mutate(label_R2 = paste0("R^2:", round(R2, 2))) %>%  
  mutate(label_RMSE = paste0("RMSE: ", round(RMSE, 0)))
print(summary_df_log_inc_CI_ppm)
write_csv(summary_df_log_inc_CI_ppm, "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Summary_log_inc_R2_RMSE_CI_ppm.csv")

# Plot Ti pred vs actual as ppm, colored by Site in ggplot -----------------------

# Make labels
R2_text_Ti_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Ti_ICP") %>% 
  pull(label_R2)
RMSE_text_Ti_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Ti_ICP") %>% 
  pull(label_RMSE)
# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_ppm_Ti <- ggplot(df_predicted_ppm, aes(x = Ti_ICP, y = Ti_ICP_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "PLSR: Measured vs Predicted Ti ppm ICP-MS log_inc model",
       x = "Measured Ti ICP-MS (ppm)",
       y = "Predicted Ti ICP-MS (ppm) [log_inc model]") +
  coord_equal(xlim = c(0, 40000), ylim = c(0, 40000)) + # coord_equal()
  #theme_minimal() +
  #coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 0, y = 36000, label = RMSE_text_Ti_ppm, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 38000, label = eq_text_Ti, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 40000, label = R2_text_Ti_ppm, parse = TRUE, hjust = 0)
print(PLS_final_ppm_Ti)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Ti_plsr_log_inc_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Ca pred vs actual as ppm, colored by Site in ggplot -----------------------

# Make labels
R2_text_Ca_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Ca_ICP") %>% 
  pull(label_R2)
RMSE_text_Ca_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Ca_ICP") %>% 
  pull(label_RMSE)
# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_ppm_Ca <- ggplot(df_predicted_ppm, aes(x = Ca_ICP, y = Ca_ICP_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "PLSR: Measured vs Predicted Ca ppm ICP-MS log_inc model",
       x = "Measured Ca ICP-MS (ppm)",
       y = "Predicted Ca ICP-MS (ppm) [log_inc model]") +
  coord_equal(xlim = c(0, 40000), ylim = c(0, 40000)) + # coord_equal()
  #theme_minimal() +
  #coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 0, y = 36000, label = RMSE_text_Ca, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 38000, label = eq_text_Ca, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 40000, label = R2_text_Ca, parse = TRUE, hjust = 0)
print(PLS_final_ppm_Ca)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Ca_plsr_log_inc_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Sr pred vs actual as ppm, colored by Site in ggplot -----------------------

# Make labels
R2_text_Sr_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Sr_ICP") %>% 
  pull(label_R2)
RMSE_text_Sr_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Sr_ICP") %>% 
  pull(label_RMSE)
# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchiSrgo"
PLS_final_ppm_Sr <- ggplot(df_predicted_ppm, aes(x = Sr_ICP, y = Sr_ICP_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "PLSR: Measured vs Predicted Sr ppm ICP-MS log_inc model",
       x = "Measured Sr ICP-MS (ppm)",
       y = "Predicted Sr ICP-MS (ppm) [log_inc model]") +
  coord_equal(xlim = c(0, 600), ylim = c(0, 600)) + # coord_equal()
  #theme_minimal() +
  #coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 0, y = 550, label = RMSE_text_Sr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 575, label = eq_text_Sr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 600, label = R2_text_Sr, parse = TRUE, hjust = 0)
print(PLS_final_ppm_Sr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Sr_plsr_log_inc_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Zr pred vs actual as ppm, colored by Site in ggplot -----------------------

# Make labels
R2_text_Zr_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Zr_ICP") %>% 
  pull(label_R2)
RMSE_text_Zr_ppm <- summary_df_log_inc_CI_ppm %>%
  filter(Element == "Zr_ICP") %>% 
  pull(label_RMSE)
# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchiZrgo"
PLS_final_ppm_Zr <- ggplot(df_predicted_ppm, aes(x = Zr_ICP, y = Zr_ICP_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "PLZr: Measured vs Predicted Zr ppm ICP-MS log_inc model",
       x = "Measured Zr ICP-MS (ppm)",
       y = "Predicted Zr ICP-MS (ppm) [log_inc model]") +
  coord_equal(xlim = c(0, 600), ylim = c(0, 600)) + # coord_equal()
  #theme_minimal() +
  #coord_equal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 0, y = 550, label = RMSE_text_Zr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 575, label = eq_text_Zr, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 600, label = R2_text_Zr, parse = TRUE, hjust = 0)
print(PLS_final_ppm_Zr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/Zr_plsr_log_inc_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# ------------------------------------------------------------------------------

# 1.3 Predicting Ti as ppm for a new xrf dataset -------------------------------

# Based on pls_final training model and Ln(Ti/inc) ITRAX dataset

# Import ACE Ti XRF-CS log_inc data and convert to ppm with +/- RMSE error
ACE_xrf_log_inc <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_log_inc.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>% 
  select(Location:MSE, all_of(main_elements_xrf), Total_scatter, inc_coh, coh_inc) %>%
  filter(qc == "TRUE") %>% 
  filter(validity == "1") %>% 
  mutate(across(all_of(main_elements_xrf), ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_log_inc

ACE_xrf_predict_PLS <- ACE_xrf_log_inc %>%
  select(Site, depth, SH20_age, Ti, Ca, Sr, Zr)

# Convert tibble to data frame for base R
X_new <- ACE_xrf_predict_PLS %>% 
  select (Ti, Ca, Sr, Zr)
df_xrf_new <- as.matrix(X_new)

# Get predicted xrf log_inc values for Ti, Ca, Sr, Zr
pred_xrf_log_inc <- predict(pls_final, newdata = df_xrf_new, ncomp = 4)
colnames(pred_xrf_log_inc) <- paste0(colnames(pred_xrf_log_inc), "_pred")

# As response was multivariate (multiple ICP columns), extract predictions as a matrix
pred_xrf_log_inc_matrix <- pred_xrf_log_inc[, , 1]  # drops 3D array to matrix (140000 x response_variables)

# Convert to tibble for easier processing
pred_log_inc_tbl <- as_tibble(pred_xrf_log_inc_matrix)

# 95% confidence intervals / error values as xrf_log_inc_new +/- RMSE values from training model

# Define R2 and RMSE values from training dataset (as above) 
RMSE_log_inc <- c(Ti_ICP_clr = rmse_Ti, Ca_ICP_clr = rmse_Ca, Sr_ICP_clr = rmse_Sr, Zr_ICP_clr = rmse_Zr)

# Simulate upper CIs for all elements in matrix by calculating prediction bounds as pred +/-RMSE from PLS model 
pred_log_inc_upper <- pred_log_inc_tbl %>%
  map2_dfc(RMSE_log_inc, ~ .x + .y)

# Simulate lower CIs for all elements in matrix by calculating prediction bounds as pred +/-RMSE from PLS model 
pred_log_inc_lower <- pred_log_inc_tbl %>%
  map2_dfc(RMSE_log_inc, ~ .x - .y)

# Rename outputs
pred_log_inc_upper <- pred_log_inc_upper %>% rename_with(~ paste0(.x, "_upper"))
pred_log_inc_lower <- pred_log_inc_lower %>% rename_with(~ paste0(.x, "_lower"))

# Original data was log-transformed - therefore take exp of both sides to convert to ppm 
pred_xrf_log_inc_exp <- exp(pred_xrf_log_inc_matrix)
colnames(pred_xrf_log_inc_exp) <- paste0(colnames(pred_xrf_log_inc_exp), "_ppm")
pred_log_inc_upper_exp <- exp(pred_log_inc_upper)
colnames(pred_log_inc_upper_exp) <- paste0(colnames(pred_log_inc_upper_exp), "_ppm")
pred_log_inc_lower_exp <- exp(pred_log_inc_lower)
colnames(pred_log_inc_lower_exp) <- paste0(colnames(pred_log_inc_lower_exp), "_ppm")

# Make matrix outputs into tibble dataframe and rename!!
df_pred_xrf_log_inc <- as_tibble(pred_xrf_log_inc) %>% 
  rename(Ti_ICP_pred = `Ti_ICP_pred.4 comps`) %>% 
  rename(Ca_ICP_pred = `Ca_ICP_pred.4 comps`) %>% 
  rename(Sr_ICP_pred = `Sr_ICP_pred.4 comps`) %>% 
  rename(Zr_ICP_pred = `Zr_ICP_pred.4 comps`)
df_pred_xrf_log_inc_exp <- as_tibble(pred_xrf_log_inc_exp)
df_pred_log_inc_lower_exp <- as_tibble(pred_log_inc_lower_exp) %>% 
  rename(Ti_low_RMSE = `Ti_ICP_pred_lower_ppm`) %>% 
  rename(Ca_low_RMSE = `Ca_ICP_pred_lower_ppm`) %>% 
  rename(Sr_low_RMSE = `Sr_ICP_pred_lower_ppm`) %>% 
  rename(Zr_low_RMSE = `Zr_ICP_pred_lower_ppm`)
df_pred_log_inc_upper_exp <- as_tibble(pred_log_inc_upper_exp) %>% 
  rename(Ti_up_RMSE = `Ti_ICP_pred_upper_ppm`) %>% 
  rename(Ca_up_RMSE = `Ca_ICP_pred_upper_ppm`) %>% 
  rename(Sr_up_RMSE = `Sr_ICP_pred_upper_ppm`) %>% 
  rename(Zr_up_RMSE = `Zr_ICP_pred_upper_ppm`)

# bind into one tibble
df_pred_combined <- tibble(df_pred_xrf_log_inc, 
                      df_pred_xrf_log_inc_exp, 
                      df_pred_log_inc_lower_exp, 
                      df_pred_log_inc_upper_exp)

# Add original dataframe and gather element data together
ACE_PLS_ppm_pred_log_inc <- bind_cols(ACE_xrf_predict_PLS, df_pred_combined) %>% 
  relocate(Ti_ICP_pred_ppm, .after = Ti_ICP_pred) %>% 
  relocate(Ti_low_RMSE, .after = Ti_ICP_pred_ppm) %>%
  relocate(Ti_up_RMSE, .after = Ti_low_RMSE) %>%
  relocate(Ca_ICP_pred_ppm, .after = Ca_ICP_pred) %>% 
  relocate(Ca_low_RMSE, .after = Ca_ICP_pred_ppm) %>%
  relocate(Ca_up_RMSE, .after = Ca_low_RMSE) %>%
  relocate(Sr_ICP_pred_ppm, .after = Sr_ICP_pred) %>% 
  relocate(Sr_low_RMSE, .after = Sr_ICP_pred_ppm) %>%
  relocate(Sr_up_RMSE, .after = Sr_low_RMSE) %>%
  relocate(Zr_ICP_pred_ppm, .after = Zr_ICP_pred) %>% 
  relocate(Zr_low_RMSE, .after = Zr_ICP_pred_ppm) %>%
  relocate(Zr_up_RMSE, .after = Zr_low_RMSE)
ACE_PLS_ppm_pred_log_inc
write.csv(ACE_PLS_ppm_pred_log_inc,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/ACE_PLS_ppm_pred_log_inc.csv", row.names = FALSE)

# ------------------------------------------------------------------------------

# 1.4 Site based plots vs depth and age ----------------------------------------

# Ti ---------------------------------------------------------------------------

# HER42PB ---------------------------------------------------------------

# Extract HER42PB log_inc predicted ppm data for Figure 4 plotting
HER42PB_PLS_ppm_pred_log_inc <- ACE_PLS_ppm_pred_log_inc %>% 
  filter(Site == "HER42PB")
HER42PB_PLS_ppm_pred_log_inc
write.csv(HER42PB_PLS_ppm_pred_log_inc,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/HER42PB_PLS_ppm_pred_log_inc.csv", row.names = FALSE)

# Import HER42PB log_inc WLS predicted ppm data for Figure 4
HER42PB_WLS_ppm_pred_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_conversion_Ti_WLS.csv") %>% 
  filter(Site == "HER42PB")
HER42PB_WLS_ppm_pred_log_inc

# Import HER42PB clr plsr predicted ppm data for Figure 4
HER42PB_PLS_ppm_pred_clr <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/HER42PB_PLS_ppm_pred_clr.csv")
HER42PB_PLS_ppm_pred_clr

# Import ICPMS Ti data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(xrf_icp_elements1)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "HER42PB")

# HER42PB - Single depth plot with RMSE errors - black and grey ----------------
ggplot(HER42PB_PLS_ppm_pred_log_inc, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_pred_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_pred_ppm)) +
  scale_y_reverse() +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/HER42PB_Ti_log_inc_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(HER42PB_PLS_ppm_pred_log_inc, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_pred_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_pred_ppm)) +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Age (cal yr BP)", y = "Ti_ICP_pred (ppm)") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/HER42PB_Ti_log_inc_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")

# HER42PB Depth - multiple plot overlaid comparing WLS and PLS outputs v --------

# Lock x and y ranges to all datasets (used later)
x_all <- range(c(HER42PB_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 HER42PB_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 HER42PB_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 HER42PB_ICP$Ti_ICP), na.rm = TRUE)
y_all <- range(c(HER42PB_PLS_ppm_pred_log_inc$depth, 
                 HER42PB_WLS_ppm_pred_log_inc$depth, 
                 HER42PB_PLS_ppm_pred_clr$depth, 
                 HER42PB_ICP$depth), na.rm = TRUE)
ggplot(HER42PB_PLS_ppm_pred_log_inc, aes(x = Ti_ICP_pred_ppm, y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_lineh(aes(x = Ti_ICP_pred_ppm, y = depth), color = "darkorange") +
  geom_point(aes(x = Ti_ICP_pred_ppm, y = depth), size = 0.5, color = "darkorange") +
  #geom_lineh(data = HER42PB_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), color = "black", alpha = 0.3) +
  #geom_point(data = HER42PB_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), size = 0.5, color = "black", alpha = 0.3) +
  geom_lineh(data = HER42PB_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), color = "brown", alpha = 0.5) +
  #geom_point(data = HER42PB_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), size = 0.5, color = "brown", alpha = 0.5) +
  geom_lineh(data = HER42PB_ICP, aes(x = Ti_ICP, y = depth), color = "blue", linewidth = 1) +
  geom_point(data = HER42PB_ICP, aes(x = Ti_ICP, y = depth), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_y_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/HER42PB_Ti_log_inc_PLS_orange_vs_clr_brown_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

# HER42PB Age - multiple plot overlaid comparing WLS and PLS outputs vs  --------

# Lock x and y ranges to all datasets (used later)
y_all <- range(c(HER42PB_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 HER42PB_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 HER42PB_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 HER42PB_ICP$Ti_ICP), na.rm = TRUE)
x_all <- range(c(HER42PB_PLS_ppm_pred_log_inc$SH20_age, 
                 HER42PB_WLS_ppm_pred_log_inc$SH20_age, 
                 HER42PB_PLS_ppm_pred_clr$SH20_age, 
                 HER42PB_ICP$depth), na.rm = TRUE)
ggplot(HER42PB_PLS_ppm_pred_log_inc, aes(y = Ti_ICP_pred_ppm, x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_line(aes(y = Ti_ICP_pred_ppm, x = SH20_age), color = "darkorange") +
  geom_point(aes(y = Ti_ICP_pred_ppm, x = SH20_age), size = 0.5, color = "darkorange") +
  #geom_line(data = HER42PB_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), color = "black", alpha = 0.3) +
  #geom_point(data = HER42PB_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), size = 0.5, color = "black", alpha = 0.3) +
  geom_line(data = HER42PB_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), color = "brown", alpha = 0.5) +
  #geom_point(data = HER42PB_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), size = 0.5, color = "brown", alpha = 0.5) +
  geom_line(data = HER42PB_ICP, aes(y = Ti_ICP, x = SH20_age), color = "blue", linewidth = 1) +
  geom_point(data = HER42PB_ICP, aes(y = Ti_ICP, x = SH20_age), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_x_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/HER42PB_Ti_log_inc_PLS_orange_vs_clr_brown.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")



# BI10 ---------------------------------------------------------------

# Extract BI10 log_inc predicted ppm data for Figure 4 plotting
BI10_PLS_ppm_pred_log_inc <- ACE_PLS_ppm_pred_log_inc %>% 
  filter(Site == "BI10")
BI10_PLS_ppm_pred_log_inc
write.csv(BI10_PLS_ppm_pred_log_inc,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/BI10_PLS_ppm_pred_log_inc.csv", row.names = FALSE)

# Import BI10 log_inc WLS predicted ppm data for Figure 4
BI10_WLS_ppm_pred_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_conversion_Ti_WLS.csv") %>% 
  filter(Site == "BI10")
BI10_WLS_ppm_pred_log_inc

# Import BI10 clr plsr predicted ppm data for Figure 4
BI10_PLS_ppm_pred_clr <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/BI10_PLS_ppm_pred_clr.csv")
BI10_PLS_ppm_pred_clr

# Import ICPMS Ti data for comparison / overlay
BI10_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(xrf_icp_elements1)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "BI10")

# BI10 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(BI10_PLS_ppm_pred_log_inc, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_pred_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_pred_ppm)) +
  scale_y_reverse() +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/BI10_Ti_log_inc_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(BI10_PLS_ppm_pred_log_inc, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_pred_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_pred_ppm)) +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Age (cal yr BP)", y = "Ti_ICP_pred (ppm)") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/BI10_Ti_log_inc_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")

# BI10 Depth - multiple plot overlaid comparing WLS and PLS outputs v --------

# Lock x and y ranges to all datasets (used later)
x_all <- range(c(BI10_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 BI10_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 BI10_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 BI10_ICP$Ti_ICP), na.rm = TRUE)
y_all <- range(c(BI10_PLS_ppm_pred_log_inc$depth, 
                 BI10_WLS_ppm_pred_log_inc$depth, 
                 BI10_PLS_ppm_pred_clr$depth, 
                 BI10_ICP$depth), na.rm = TRUE)
ggplot(BI10_PLS_ppm_pred_log_inc, aes(x = Ti_ICP_pred_ppm, y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_lineh(aes(x = Ti_ICP_pred_ppm, y = depth), color = "darkorange") +
  geom_point(aes(x = Ti_ICP_pred_ppm, y = depth), size = 0.5, color = "darkorange") +
  #geom_lineh(data = BI10_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), color = "black", alpha = 0.3) +
  #geom_point(data = BI10_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), size = 0.5, color = "black", alpha = 0.3) +
  geom_lineh(data = BI10_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), color = "brown", alpha = 0.5) +
  #geom_point(data = BI10_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), size = 0.5, color = "brown", alpha = 0.5) +
  geom_lineh(data = BI10_ICP, aes(x = Ti_ICP, y = depth), color = "blue", linewidth = 1) +
  geom_point(data = BI10_ICP, aes(x = Ti_ICP, y = depth), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_y_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/BI10_Ti_log_inc_PLS_orange_vs_clr_brown_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

# BI10 Age - multiple plot overlaid comparing WLS and PLS outputs vs  --------

# Lock x and y ranges to all datasets (used later)
y_all <- range(c(BI10_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 BI10_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 BI10_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 BI10_ICP$Ti_ICP), na.rm = TRUE)
x_all <- range(c(BI10_PLS_ppm_pred_log_inc$SH20_age, 
                 BI10_WLS_ppm_pred_log_inc$SH20_age, 
                 BI10_PLS_ppm_pred_clr$SH20_age, 
                 BI10_ICP$depth), na.rm = TRUE)
ggplot(BI10_PLS_ppm_pred_log_inc, aes(y = Ti_ICP_pred_ppm, x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_line(aes(y = Ti_ICP_pred_ppm, x = SH20_age), color = "darkorange") +
  geom_point(aes(y = Ti_ICP_pred_ppm, x = SH20_age), size = 0.5, color = "darkorange") +
  #geom_line(data = BI10_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), color = "black", alpha = 0.3) +
  #geom_point(data = BI10_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), size = 0.5, color = "black", alpha = 0.3) +
  geom_line(data = BI10_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), color = "brown", alpha = 0.5) +
  #geom_point(data = BI10_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), size = 0.5, color = "brown", alpha = 0.5) +
  geom_line(data = BI10_ICP, aes(y = Ti_ICP, x = SH20_age), color = "blue", linewidth = 1) +
  geom_point(data = BI10_ICP, aes(y = Ti_ICP, x = SH20_age), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_x_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/BI10_Ti_log_inc_PLS_orange_vs_clr_brown.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")






# KER1 ---------------------------------------------------------------

# Extract KER1 log_inc predicted ppm data for Figure 4 plotting
KER1_PLS_ppm_pred_log_inc <- ACE_PLS_ppm_pred_log_inc %>% 
  filter(Site == "KER1")
KER1_PLS_ppm_pred_log_inc
write.csv(KER1_PLS_ppm_pred_log_inc,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER1_PLS_ppm_pred_log_inc.csv", row.names = FALSE)

# Import KER1 log_inc WLS predicted ppm data for Figure 4
KER1_WLS_ppm_pred_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_conversion_Ti_WLS.csv") %>% 
  filter(Site == "KER1")
KER1_WLS_ppm_pred_log_inc

# Import KER1 clr plsr predicted ppm data for Figure 4
KER1_PLS_ppm_pred_clr <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/KER1_PLS_ppm_pred_clr.csv")
KER1_PLS_ppm_pred_clr

# Import ICPMS Ti data for comparison / overlay
KER1_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(xrf_icp_elements1)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER1")

# KER1 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(KER1_PLS_ppm_pred_log_inc, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_pred_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_pred_ppm)) +
  scale_y_reverse() +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER1_Ti_log_inc_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(KER1_PLS_ppm_pred_log_inc, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_pred_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_pred_ppm)) +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Age (cal yr BP)", y = "Ti_ICP_pred (ppm)") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER1_Ti_log_inc_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")

# KER1 Depth - multiple plot overlaid comparing WLS and PLS outputs v --------

# Lock x and y ranges to all datasets (used later)
x_all <- range(c(KER1_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 KER1_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 KER1_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 KER1_ICP$Ti_ICP), na.rm = TRUE)
y_all <- range(c(KER1_PLS_ppm_pred_log_inc$depth, 
                 KER1_WLS_ppm_pred_log_inc$depth, 
                 KER1_PLS_ppm_pred_clr$depth, 
                 KER1_ICP$depth), na.rm = TRUE)
ggplot(KER1_PLS_ppm_pred_log_inc, aes(x = Ti_ICP_pred_ppm, y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_lineh(aes(x = Ti_ICP_pred_ppm, y = depth), color = "darkorange") +
  geom_point(aes(x = Ti_ICP_pred_ppm, y = depth), size = 0.5, color = "darkorange") +
  #geom_lineh(data = KER1_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), color = "black", alpha = 0.3) +
  #geom_point(data = KER1_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), size = 0.5, color = "black", alpha = 0.3) +
  geom_lineh(data = KER1_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), color = "brown", alpha = 0.5) +
  #geom_point(data = KER1_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), size = 0.5, color = "brown", alpha = 0.5) +
  geom_lineh(data = KER1_ICP, aes(x = Ti_ICP, y = depth), color = "blue", linewidth = 1) +
  geom_point(data = KER1_ICP, aes(x = Ti_ICP, y = depth), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_y_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER1_Ti_log_inc_PLS_orange_vs_clr_brown_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

# KER1 Age - multiple plot overlaid comparing WLS and PLS outputs vs  --------

# Lock x and y ranges to all datasets (used later)
y_all <- range(c(KER1_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 KER1_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 KER1_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 KER1_ICP$Ti_ICP), na.rm = TRUE)
x_all <- range(c(KER1_PLS_ppm_pred_log_inc$SH20_age, 
                 KER1_WLS_ppm_pred_log_inc$SH20_age, 
                 KER1_PLS_ppm_pred_clr$SH20_age, 
                 KER1_ICP$depth), na.rm = TRUE)
ggplot(KER1_PLS_ppm_pred_log_inc, aes(y = Ti_ICP_pred_ppm, x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_line(aes(y = Ti_ICP_pred_ppm, x = SH20_age), color = "darkorange") +
  geom_point(aes(y = Ti_ICP_pred_ppm, x = SH20_age), size = 0.5, color = "darkorange") +
  #geom_line(data = KER1_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), color = "black", alpha = 0.3) +
  #geom_point(data = KER1_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), size = 0.5, color = "black", alpha = 0.3) +
  geom_line(data = KER1_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), color = "brown", alpha = 0.5) +
  #geom_point(data = KER1_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), size = 0.5, color = "brown", alpha = 0.5) +
  geom_line(data = KER1_ICP, aes(y = Ti_ICP, x = SH20_age), color = "blue", linewidth = 1) +
  geom_point(data = KER1_ICP, aes(y = Ti_ICP, x = SH20_age), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_x_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER1_Ti_log_inc_PLS_orange_vs_clr_brown.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")








# KER3 ---------------------------------------------------------------

# Extract KER3 log_inc predicted ppm data for Figure 4 plotting
KER3_PLS_ppm_pred_log_inc <- ACE_PLS_ppm_pred_log_inc %>% 
  filter(Site == "KER3")
KER3_PLS_ppm_pred_log_inc
write.csv(KER3_PLS_ppm_pred_log_inc,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER3_PLS_ppm_pred_log_inc.csv", row.names = FALSE)

# Import KER3 log_inc WLS predicted ppm data for Figure 4
KER3_WLS_ppm_pred_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_conversion_Ti_WLS.csv") %>% 
  filter(Site == "KER3")
KER3_WLS_ppm_pred_log_inc

# Import KER3 clr plsr predicted ppm data for Figure 4
KER3_PLS_ppm_pred_clr <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/KER3_PLS_ppm_pred_clr.csv")
KER3_PLS_ppm_pred_clr

# Import ICPMS Ti data for comparison / overlay
KER3_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(xrf_icp_elements1)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER3")

# KER3 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(KER3_PLS_ppm_pred_log_inc, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_pred_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_pred_ppm)) +
  scale_y_reverse() +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER3_Ti_log_inc_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(KER3_PLS_ppm_pred_log_inc, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_pred_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_pred_ppm)) +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Age (cal yr BP)", y = "Ti_ICP_pred (ppm)") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER3_Ti_log_inc_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")

# KER3 Depth - multiple plot overlaid comparing WLS and PLS outputs v --------

# Lock x and y ranges to all datasets (used later)
x_all <- range(c(KER3_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 KER3_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 KER3_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 KER3_ICP$Ti_ICP), na.rm = TRUE)
y_all <- range(c(KER3_PLS_ppm_pred_log_inc$depth, 
                 KER3_WLS_ppm_pred_log_inc$depth, 
                 KER3_PLS_ppm_pred_clr$depth, 
                 KER3_ICP$depth), na.rm = TRUE)
ggplot(KER3_PLS_ppm_pred_log_inc, aes(x = Ti_ICP_pred_ppm, y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_lineh(aes(x = Ti_ICP_pred_ppm, y = depth), color = "darkorange") +
  geom_point(aes(x = Ti_ICP_pred_ppm, y = depth), size = 0.5, color = "darkorange") +
  #geom_lineh(data = KER3_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), color = "black", alpha = 0.3) +
  #geom_point(data = KER3_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), size = 0.5, color = "black", alpha = 0.3) +
  geom_lineh(data = KER3_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), color = "brown", alpha = 0.5) +
  #geom_point(data = KER3_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), size = 0.5, color = "brown", alpha = 0.5) +
  geom_lineh(data = KER3_ICP, aes(x = Ti_ICP, y = depth), color = "blue", linewidth = 1) +
  geom_point(data = KER3_ICP, aes(x = Ti_ICP, y = depth), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_y_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER3_Ti_log_inc_PLS_orange_vs_clr_brown_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

# KER3 Age - multiple plot overlaid comparing WLS and PLS outputs vs  --------

# Lock x and y ranges to all datasets (used later)
y_all <- range(c(KER3_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 KER3_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 KER3_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 KER3_ICP$Ti_ICP), na.rm = TRUE)
x_all <- range(c(KER3_PLS_ppm_pred_log_inc$SH20_age, 
                 KER3_WLS_ppm_pred_log_inc$SH20_age, 
                 KER3_PLS_ppm_pred_clr$SH20_age, 
                 KER3_ICP$depth), na.rm = TRUE)
ggplot(KER3_PLS_ppm_pred_log_inc, aes(y = Ti_ICP_pred_ppm, x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_line(aes(y = Ti_ICP_pred_ppm, x = SH20_age), color = "darkorange") +
  geom_point(aes(y = Ti_ICP_pred_ppm, x = SH20_age), size = 0.5, color = "darkorange") +
  #geom_line(data = KER3_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), color = "black", alpha = 0.3) +
  #geom_point(data = KER3_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), size = 0.5, color = "black", alpha = 0.3) +
  geom_line(data = KER3_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), color = "brown", alpha = 0.5) +
  #geom_point(data = KER3_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), size = 0.5, color = "brown", alpha = 0.5) +
  geom_line(data = KER3_ICP, aes(y = Ti_ICP, x = SH20_age), color = "blue", linewidth = 1) +
  geom_point(data = KER3_ICP, aes(y = Ti_ICP, x = SH20_age), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_x_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/KER3_Ti_log_inc_PLS_orange_vs_clr_brown.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")









# PB1 ---------------------------------------------------------------

# Extract PB1 log_inc predicted ppm data for Figure 4 plotting
PB1_PLS_ppm_pred_log_inc <- ACE_PLS_ppm_pred_log_inc %>% 
  filter(Site == "PB1")
PB1_PLS_ppm_pred_log_inc
write.csv(PB1_PLS_ppm_pred_log_inc,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/PB1_PLS_ppm_pred_log_inc.csv", row.names = FALSE)

# Import PB1 log_inc WLS predicted ppm data for Figure 4
PB1_WLS_ppm_pred_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_conversion_Ti_WLS.csv") %>% 
  filter(Site == "PB1")
PB1_WLS_ppm_pred_log_inc

# Import PB1 clr plsr predicted ppm data for Figure 4
PB1_PLS_ppm_pred_clr <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/PB1_PLS_ppm_pred_clr.csv")
PB1_PLS_ppm_pred_clr

# Import ICPMS Ti data for comparison / overlay
PB1_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(xrf_icp_elements1)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "PB1")

# PB1 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(PB1_PLS_ppm_pred_log_inc, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_pred_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_pred_ppm)) +
  scale_y_reverse() +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/PB1_Ti_log_inc_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(PB1_PLS_ppm_pred_log_inc, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_pred_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_pred_ppm)) +
  labs(title = "PLSR Predictions with Confidence Interval (Ti)",
       x = "Age (cal yr BP)", y = "Ti_ICP_pred (ppm)") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/PB1_Ti_log_inc_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")

# PB1 Depth - multiple plot overlaid comparing WLS and PLS outputs v --------

# Lock x and y ranges to all datasets (used later)
x_all <- range(c(PB1_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 PB1_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 PB1_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 PB1_ICP$Ti_ICP), na.rm = TRUE)
y_all <- range(c(PB1_PLS_ppm_pred_log_inc$depth, 
                 PB1_WLS_ppm_pred_log_inc$depth, 
                 PB1_PLS_ppm_pred_clr$depth, 
                 PB1_ICP$depth), na.rm = TRUE)
ggplot(PB1_PLS_ppm_pred_log_inc, aes(x = Ti_ICP_pred_ppm, y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_lineh(aes(x = Ti_ICP_pred_ppm, y = depth), color = "darkorange") +
  geom_point(aes(x = Ti_ICP_pred_ppm, y = depth), size = 0.5, color = "darkorange") +
  #geom_lineh(data = PB1_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), color = "black", alpha = 0.3) +
  #geom_point(data = PB1_WLS_ppm_pred_log_inc, aes(x = Ti_ppm_WLS, y = depth), size = 0.5, color = "black", alpha = 0.3) +
  geom_lineh(data = PB1_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), color = "brown", alpha = 0.5) +
  #geom_point(data = PB1_PLS_ppm_pred_clr, aes(x = Ti_ICP_clr_ppm, y = depth), size = 0.5, color = "brown", alpha = 0.5) +
  geom_lineh(data = PB1_ICP, aes(x = Ti_ICP, y = depth), color = "blue", linewidth = 1) +
  geom_point(data = PB1_ICP, aes(x = Ti_ICP, y = depth), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_y_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/PB1_Ti_log_inc_PLS_orange_vs_clr_brown_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

# PB1 Age - multiple plot overlaid comparing WLS and PLS outputs vs  --------

# Lock x and y ranges to all datasets (used later)
y_all <- range(c(PB1_PLS_ppm_pred_log_inc$Ti_ICP_pred_ppm, 
                 PB1_WLS_ppm_pred_log_inc$Ti_ppm_WLS, 
                 PB1_PLS_ppm_pred_clr$Ti_ICP_clr_ppm, 
                 PB1_ICP$Ti_ICP), na.rm = TRUE)
x_all <- range(c(PB1_PLS_ppm_pred_log_inc$SH20_age, 
                 PB1_WLS_ppm_pred_log_inc$SH20_age, 
                 PB1_PLS_ppm_pred_clr$SH20_age, 
                 PB1_ICP$depth), na.rm = TRUE)
ggplot(PB1_PLS_ppm_pred_log_inc, aes(y = Ti_ICP_pred_ppm, x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_line(aes(y = Ti_ICP_pred_ppm, x = SH20_age), color = "darkorange") +
  geom_point(aes(y = Ti_ICP_pred_ppm, x = SH20_age), size = 0.5, color = "darkorange") +
  #geom_line(data = PB1_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), color = "black", alpha = 0.3) +
  #geom_point(data = PB1_WLS_ppm_pred_log_inc, aes(y = Ti_ppm_WLS, x = SH20_age), size = 0.5, color = "black", alpha = 0.3) +
  geom_line(data = PB1_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), color = "brown", alpha = 0.5) +
  #geom_point(data = PB1_PLS_ppm_pred_clr, aes(y = Ti_ICP_clr_ppm, x = SH20_age), size = 0.5, color = "brown", alpha = 0.5) +
  geom_line(data = PB1_ICP, aes(y = Ti_ICP, x = SH20_age), color = "blue", linewidth = 1) +
  geom_point(data = PB1_ICP, aes(y = Ti_ICP, x = SH20_age), shape = 21, size = 1, color = "blue", fill = "white", stroke = 1) +
  scale_x_continuous(limits = x_all) +
  scale_y_continuous(limits = y_all) +
  scale_x_reverse() +
  #abs(title = "PLSR Predictions with Confidence Interval (Ti)"),
  #    x = "Depth (cm)", y = "Ti_ICP_pred (ppm)") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/PB1_Ti_log_inc_PLS_orange_vs_clr_brown.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")










# END --------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 2. PLS Train-Test Split tidyverse method  ------------------------------------

# If the primary goal is to build a highly accurate predictive model, a larger training set might be preferred
# If the focus is on evaluating the model's generalizability, a larger testing set might be more appropriate
# As above but with 70% training vs 30% test dataset 

# Define df_train as training dataset and elements to use for Ti_ICP response model
df_Site <- df0[, (names(df0) %in% c("Site", "Ti_ICP", final_elements_xrf))] # or use 6 key_elements_xrf 

df_split <- df_Site %>%
  group_by(Site) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

train_df <- df_split %>%
  group_by(Site) %>%
  slice_sample(prop = 0.6) %>% # set to 0.6: 0.4 split as this matches r2_Ti = 0.75 and rmse_Ti from whole dataset final model
  ungroup()

test_df <- anti_join(df_split, train_df, by = c("Site", "row_id"))


#  2.1 Ti training:test & bootstrapping -----------------------------------

# Prepare matrices
X_train <- train_df %>% select(all_of(final_elements_xrf)) %>% as.matrix()
Y_train <- train_df$Ti_ICP

# Fit PLSR model
pls_model <- plsr(Y_train ~ X_train, scale = TRUE, validation = "CV", ncomp = 3)

# Optimal number of components
validationplot(pls_model, val.type = "MSEP")
opt_comp_train_CV <- which.min(RMSEP(pls_model)$val[1, , -1])
cat("Optimal number of components:", opt_comp_train_CV, "\n")
opt_comp <- opt_comp_train_CV

# Predict on Test Set
X_test <- test_df %>% select(all_of(final_elements_xrf)) %>% as.matrix()
Y_test <- test_df$Ti_ICP

Y_pred <- predict(pls_model, newdata = X_test, ncomp = opt_comp)[,,1]

# Add predictions to test set
Ti_ICP_test_df <- test_df %>%
  mutate(Ti_ICP_pred = Y_pred) %>% 
  relocate(row_id, .before = Site) %>%
  relocate(Ti_ICP, .before = Ti_ICP_pred)
Ti_ICP_test_df
write.csv(Ti_ICP_test_df,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/Ti_ICP_predicted.csv", row.names = FALSE)

# Plot PLS Scores (Training Set) Colored by Site
# Get scores (first 2 components)
scores_df <- as.data.frame(scores(pls_model)[, 1:2])
colnames(scores_df) <- c("Comp1", "Comp2")

# Combine with Site info
train_plot_df <- train_df %>%
  select(Site) %>%
  bind_cols(scores_df)

# Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
ggplot(train_plot_df, aes(x = Comp1, y = Comp2, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PLSR Score Plot (Training Data)", x = "PLS Component 1", y = "PLS Component 2") +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
  #coord_equal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/Ti_training_comp1vs2.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
  #theme_minimal()

# Compute R² and rmse_Ti on test set
SS_res <- sum((Y_test - Y_pred)^2)
SS_tot <- sum((Y_test - mean(Y_test))^2)
R2_test <- 1 - SS_res / SS_tot
RMSE_test <- sqrt(mean((Y_test - Y_pred)^2))
cat("Optimal Components:", opt_comp, "\n")
cat("Test R²:", round(R2_test, 4), "\n")
cat("Test rmse_Ti:", round(RMSE_test, 4), "\n")

# Extract coefficients for optimal number of components
coefs <- coef(pls_model, ncomp = opt_comp)
coefs_vec <- as.vector(coefs)
# Match predictor names and print to console
names(coefs_vec) <- colnames(X_train)
print(round(coefs_vec, 4))

# Predicted vs Measured Plot
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
ggplot(Ti_ICP_test_df, aes(x = Ti_ICP, y = Ti_ICP_pred, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Titanium (Ti) (30% Test Set))",
    x = "Measured Ti (ICP-MS)",
    y = "Predicted Ti (PLS Model)"
  ) +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
#coord_equal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/Ti_pred_meas.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
#theme_minimal()
# Regression Equation

# Round coefficients
coefs_rounded <- round(coefs_vec, 4)

# Build equation string
equation <- paste0("Ti_ICP_predicted = ",
                   paste(paste0(coefs_rounded, " * ", names(coefs_rounded)), collapse = " + "))
cat("PLSR Regression Equation (", opt_comp, " components):\n", equation, "\n")


# Bootstrapping to assess errors
# Set number of bootstrap iterations
n_boot <- 1000
rmse_vals <- numeric(n_boot)
train_data <-  train_df %>% select(-c(Site, row_id, Ti_ICP)) %>% as.matrix()
test_data <-  test_df %>% select(-c(Site, row_id, Ti_ICP)) %>% as.matrix()
Y <- train_df$Ti_ICP
set.seed(1)  # for reproducibility

for (i in 1:n_boot) {
  # 1. Sample training data with replacement
  boot_idx <- sample(1:nrow(train_data), replace = TRUE)
  boot_data <- train_data[boot_idx, ]
  boot_data <- as.data.frame(boot_data)
  
  # 2. Fit PLSR model
  model_boot <- plsr(Y ~ ., data = boot_data, scale = TRUE, ncomp = 3)
  
  # 3. Predict on original (unbooted) test set
  pred_scaled <- predict(model_boot, newdata = test_data, ncomp = 3)
  pred_scaled <- as.vector(pred_scaled)
  pred <- pred_scaled
  
  # 4. Calculate rmse_Ti in original units
  actual <- test_df$Ti_ICP
  rmse_vals[i] <- sqrt(mean((actual - pred)^2))
}

## Analyze the rmse_Ti distribution

# Summary statistics
summary(rmse_vals)
mean_rmse <- mean(rmse_vals)
log(mean_rmse)

# 95% confidence interval
quantile(rmse_vals, probs = c(0.025, 0.975))

# Plot the rmse_Ti distribution
rmse_df <- tibble(RMSE = rmse_vals)
#hist(rmse_vals, breaks = 30, col = "skyblue", main = "Bootstrap rmse_Ti Distribution", xlab = "rmse_Ti")

# Plot histogram in ggplot
ggplot(rmse_df, aes(x = RMSE)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Bootstrapped RMSE Values for Ti 0.6 train:test",
       x = "Ti RMSE",
       y = "Frequency") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/Ti_RMSE_hist_bootstrapping.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ------------------------------------------------------------------------------

# Predict Ti as ppm using pls_final model and Ln(Ti/inc) ITRAX dataset ---------

# Import ACE Ti XRF-CS log_inc data and convert to ppm with +/- rmse_Ti error
ACE_xrf_log_inc <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_log_inc.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>% 
  select(Location:MSE, all_of(main_elements_xrf), Total_scatter, inc_coh, coh_inc) %>%
  filter(qc == "TRUE") %>% 
  filter(validity == "1") %>% 
  mutate(across(all_of(main_elements_xrf), ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_log_inc

ACE_xrf_predict_Ti_PLS <- ACE_xrf_log_inc %>%
  select(Site, depth, SH20_age, Ti, Ca, Sr)

# Convert tibble to data frame for base R
X_new <- ACE_xrf_predict_Ti_PLS %>% 
  select (Ti, Ca, Sr) %>% 
  as.data.frame(ACE_xrf_predict_Ti_PLS)

# Predict Ti using pls_final model and bootstrap mean rmse_Ti value - for ppm reconstruction
y_pred_scaled_bs <- predict(pls_final, newdata = X_new, ncomp = 3)
y_pred_scaled_RMSE_upper_bs <- y_pred_scaled + mean_rmse # mean rmse_Ti from Bstrap
y_pred_scaled_RMSE_lower_bs <- y_pred_scaled - mean_rmse # mean rmse_Ti from Bstrap
# Y was log-transformed - therefore take exp of both sides to convert to ppm 
y_pred_scaled_exp_bs <- exp(y_pred_scaled_bs)
y_pred_scaled_RMSE_upper_exp_bs <- exp(y_pred_scaled_RMSE_upper_bs)
y_pred_scaled_RMSE_lower_exp_bs <- exp(y_pred_scaled_RMSE_lower_bs)
y_pred_scaled_exp_t <- tibble(y_pred_scaled_exp_bs, y_pred_scaled_RMSE_upper_exp_bs,y_pred_scaled_RMSE_lower_exp_bs)
y_pred_scaled_exp_t
# Just for comparison - calculate decentered and descaled predicted Ti using mean and sd scale parameters of training Y dataset
# Get original scale parameters from training Y - center = TRUE, scale = FALSE in plsr by not scaled by default
# y_mean <- mean(Y_test) # only apply if not scaled as here
# y_sd   <- sd(Y_test)
# Descale & decenter then take exponential
# y_pred_descaled <- y_pred_scaled * log(y_sd) + log(y_mean)
# rmse_descaled <- sqrt(mean((y_pred_descaled - y_pred_scaled)^2, na.rm = TRUE))
# rmse_descaled_exp <- exp(rmse_descaled)
# y_pred_descaled_RMSE_upper <- y_pred_scaled + rmse_descaled
# y_pred_descaled_RMSE_lower <- y_pred_scaled - rmse_descaled
# y_pred_descaled_exp <- exp(y_pred_descaled)
# y_pred_descaled_RMSE_upper_exp <- exp(y_pred_scaled_RMSE_upper)
# y_pred_descaled_RMSE_lower_exp <- exp(y_pred_scaled_RMSE_lower)

# Make into tibble for export and plotting & write to file
Ti_ppm_PLS1_bs <- tibble(Ti_pred_log_inc_bs = y_pred_scaled_bs,
                      Ti_ppm_PLS_log_inc_bs = y_pred_scaled_exp_bs, 
                      Ti_upper_PLS_log_inc_bs = y_pred_scaled_RMSE_upper_exp_bs,
                      Ti_lower_PLS_log_inc_bs = y_pred_scaled_RMSE_lower_exp_bs)
Ti_ppm_PLS_predicted_bs <- bind_cols(ACE_xrf_predict_Ti_PLS, Ti_ppm_PLS1_bs)
write.csv(Ti_ppm_PLS_predicted_bs,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/ACE_Ti_ppm_PLS_predicted_bstrap.csv", row.names = FALSE)

# Extract HER42PB data for Figure 4
HER42PB_Ti_ppm_PLS_predicted_bs <- Ti_ppm_PLS_predicted_bs %>% 
  filter(Site == "HER42PB")
HER42PB_Ti_ppm_PLS_predicted_bs
write.csv(HER42PB_Ti_ppm_PLS_predicted_bs,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/HER42PB_Ti_ppm_PLS_predicted_bstrap.csv", row.names = FALSE)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. PLS Train-Test Split Base R & ggplot method [OPTIONAL] --------------------

# define df_train as training dataset and elements to use for Ti_ICP response model
df_train <- df[, (names(df) %in% c("Site", "Ti_ICP", key_elements_xrf))]

# Split dataset into training (70%) and testing (30%)
sample_size <- floor(0.7 * nrow(df_train))
train_indices <- sample(seq_len(nrow(df_train)), size = sample_size)

# create traingin and test datasets with Site name include for plotting later
train_data_site <- df_train[train_indices, ]
test_data_site <- df_train[-train_indices, ]

# remove site name for plsr analysis
train_data <- train_data_site
train_data$Site <- NULL # remove Site column
test_data <- test_data_site
test_data$Site <- NULL # remove Site column

# Write summary outputs from console to txt file in Output folder
# sink(file = "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/Summary_stats.txt")

# Fit PLS model with LOO on training data (useful for small datasets)
pls_train_model_LOO <- plsr(Ti_ICP ~ ., data = train_data, validation = "CV")

# Find & plot optimal number of components (based on lowest LOO)
validationplot(pls_train_model_LOO, val.type = "MSEP")
opt_comp_train_LOO <- which.min(RMSEP(pls_train_model_LOO)$val[1, , -1])
cat("Optimal number of components:", opt_comp_train_LOO, "\n")

# Fit PLS model with 10-fold cross-validation on training data (useful for small datasets)
pls_train_model_CV <- plsr(Ti_ICP ~ ., data = train_data, validation = "CV")

# Find & plot optimal number of components (based on lowest CV)
validationplot(pls_train_model_CV, val.type = "MSEP")
opt_comp_train_CV <- which.min(RMSEP(pls_train_model_CV)$val[1, , -1])
cat("Optimal number of components:", opt_comp_train_CV, "\n")

# Extract regression coefficients for optimal components
coefs <- coef(pls_train_model_CV, ncomp = opt_comp_train_CV)
coefs_vec <- as.vector(coefs)
names(coefs_vec) <- dimnames(coefs)[[1]]

# Predict on test set
predictions <- predict(pls_train_model_CV, newdata = test_data, ncomp = opt_comp_train_CV)
actuals <- test_data$Ti_ICP

# Plot predicted vs actual on test dataset
plot(actuals, predictions, xlab = "Measured Ti (ICPMS)", ylab = "Predicted Ti (from ITRAX)", main = "PLS Calibration")
abline(0, 1, col = "red")

# Calculate r2_Ti and rmse_Ti values
r2_Ti <- cor(actuals, predictions)^2
rmse_Ti <- sqrt(mean((actuals - predictions)^2))
cat("Test R² =", round(r2_Ti, 2), "\nTest rmse_Ti =", round(rmse_Ti, 4), "\n")

# Show regression equation
coefs <- coef(pls_train_model_CV, ncomp = opt_comp_train_CV)[,,1]
cat("PLS Regression Equation:\n")
cat("Ti_ICPMS = ", paste0(round(coefs, 4), "*", names(coefs), collapse = " + "), "\n") # removed: intercept, " + ", 

# Write console summary info to file 
# sink(file = NULL)

# Create a tibble with site name and predicted output for plotting in ggplot
df_test_joined <- cbind(test_data_site, predictions)
df_test_predicted <- as_tibble(df_test_joined) %>% 
  rename(Ti_predicted = last_col())
write.csv(df_test_joined,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR_train_test/Ti/Ti_ICP_predicted_baseR.csv", row.names = FALSE)

# Plot predictions vs actual, colored by Site in ggplot

ggplot(df_test_predicted, aes(x = Ti_ICP, y = Ti_predicted, color = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Titanium (Ti) (30% Test Set)",
    x = "Measured Ti (ICP-MS)",
    y = "Predicted Ti (PLS Model)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
# ------------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Other code ---------------------------------------------------

# -------------------------------------------------------------------------
# RMSE and errors from bootstrapping  -----------------------------------



# Set bootstrap parameters
n_boot <- 1000
set.seed(42)

# convert original training dataset to matrix
df_xrf_m <- as.matrix(df_xrf)
df_icp_m <- as.matrix(df_icp)

# Placeholder to store predictions
boot_preds <- matrix(NA, nrow = nrow(df_xrf_new), ncol = n_boot)

# Bootstrap loop
for (i in 1:n_boot) {
  # Sample with replacement from training data
  boot_idx <- sample(1:nrow(df_xrf_m), replace = TRUE)
  
  # Fit PLSR model on bootstrap sample
  boot_model <- plsr(df_icp_m[boot_idx, "Ti_ICP"] ~ df_xrf_m[boot_idx, ],
                     ncomp = opt_comp_final, scale = TRUE)
  
  # Predict on new XRF dataset
  boot_preds[, i] <- predict(boot_model, newdata = df_xrf_new, ncomp = opt_comp_final)[, , 1]
}

# --- Compute summary statistics ---
# Mean prediction
pred_mean <- rowMeans(boot_preds)

# RMSE per prediction
pred_rmse <- apply(boot_preds, 1, sd)

# Confidence intervals (5th and 95th percentiles)
pred_lower <- apply(boot_preds, 1, quantile, probs = 0.05)
pred_upper <- apply(boot_preds, 1, quantile, probs = 0.95)

# --- Combine into a final results dataframe ---
results_df <- xrf_new %>%
  mutate(
    Ti_pred_mean = pred_mean,
    Ti_RMSE = pred_rmse,
    Ti_CI_lower = pred_lower,
    Ti_CI_upper = pred_upper
  )


# PLS calibration: convert XRF-CS data to ppm from log_inc input ---------------

# Apply equation to all ITRAX data to reconstruct xrf as ppm --------------

coefs_1 <- coefs_final %>% 
  pull(1)
coefs_2 <- coefs_final %>% 
  pull(2)
coefs_3 <- coefs_final %>% 
  pull(3)
#coefs_4 <- coefs_final %>% 
#  pull(4)
#coefs_5 <- coefs_final %>% 
#  pull(5)

# define coefficients and rmse_Ti values to use in equations - see ACE.PLS.R file outputs for details
coefs_1 = 0.8160493 # Ti
coefs_2 = -0.5945589 # Ca
coefs_3 = 0.3389237 # Sr
rmse_Ti = 0.6200723 # rmse_Ti value

ACE_xrf_conversion_Ti_PLS <- ACE_xrf_log_inc %>%
  select(Site, depth, SH20_age, Ti, Ca, Sr) %>%
  filter (Site == "HER42PB") %>% 
  mutate(Ti_convert_PLS = (coefs_1**Ti) + (coefs_2*Ca) + (coefs_3*Sr)) %>% #
  mutate(Ti_ppm_PLS = exp(Ti_convert_PLS)) %>% 
  mutate(Ti_upper_PLS = Ti_convert_PLS+rmse_Ti) %>% # rmse_Ti = 0.6201
  mutate(Ti_upper_RMSE_PLS = exp(Ti_upper_PLS)) %>% 
  mutate(Ti_lower_PLS = Ti_convert_PLS-rmse_Ti) %>% 
  mutate(Ti_lower_RMSE_PLS = exp(Ti_lower_PLS)) %>% 
  select(Site, depth, SH20_age, Ti, Ca, Sr, Ti_convert_PLS, Ti_ppm_PLS, Ti_lower_RMSE_PLS, Ti_upper_RMSE_PLS)
ACE_xrf_conversion_Ti_PLS
write.csv(ACE_xrf_conversion_Ti_PLS,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/ACE_plsr_xrf_Ti_converted_PLS.csv", row.names = FALSE)


# -------------------------------------------------------------------------


