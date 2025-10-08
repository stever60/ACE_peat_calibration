# ACE PLSR clr using pls package 

# XRF clr model and Ti_ICP response variable to scaled and centred clr dataset

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

# Import the matched ICPMS-XRF cps dataset and create a PLS-model for prediction of Ti_ICP
# Using up to 11 other elements as predictors. 

# Convert cps to clr - rename cps elements as element_clr

# Assumptions:
# Ti_ICPMS is the response variable.
# Other elements (e.g., Ca, Fe, Sr, Mn, etc.) measured by ITRAX are predictor variables.
# All data are numeric and cleaned (no missing values).
# The response variable is continuous.

# Load and Prepare data
ACE_matched_xrf_icp_cps <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") 
is.na(ACE_matched_xrf_icp_cps)<-sapply(ACE_matched_xrf_icp_cps, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_cps

# Select dataset and convert to base R dataframe for pls package ---------------
df0 <- ACE_matched_xrf_icp_cps %>%
  select(c(Site, depth, SH20_mean_age, SH20_mean_95CI, all_of(xrf_icp_elements1))) %>%
  drop_na()
df0

# Convert tibble to data frame for base R
df <- as.data.frame(df0)

#-------------------------------------------------------------------------------

# 1. Partial Least Squares (PLS) regression ------------------------------------

# ------------------------------------------------------------------------------

# 1.1 Validation based on log_inc - use Ti, Ca, Sr, Zr  ------------------------

# ------------------------------------------------------------------------------

# 1.2 PSLR training model with Ti, Ca, Sr, Zr as element matrix ----------------

# Define predictor ITRAX variables & ICPMS response variable for calibration model tests
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
main_elements_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP","Fe_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP","Zr_ICP")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")

# Input XRF and ICPMS ---

# Example: assume df_xrf and df_icp are your original data frames

df_key_elements_xrf <- df[, c(key_elements_xrf)]
df_key_elements_icp <- df[, c(key_elements_icp)]

df_xrf <- df[, c("Ti", "Ca", "Sr", "Zr")]
df_icp <- df[, c("Ti_ICP", "Ca_ICP", "Sr_ICP", "Zr_ICP")]

# CLR-transform predictors and response ---

# Transform to compositions and add _clr to column titles
xrf_clr <- clr(acomp(df_xrf))
colnames(xrf_clr) <- paste0(colnames(xrf_clr), "_clr")
icp_clr <- clr(acomp(df_icp))
colnames(icp_clr) <- paste0(colnames(icp_clr), "_clr")

# Refit CV-PLSR (cross-validated) training model -------------------------------
plsr_model <- plsr(icp_clr ~ xrf_clr, ncomp = 4, validation = "CV", jackknife = TRUE, scale = TRUE, center = TRUE)

# Summary of final plsr model and jackknife t tests of regression coefficients & std errors
summary(plsr_model)
jack.test(plsr_model, ncomp = 4, use.mean = TRUE) # estimates and std errors not valid

# Generate R2, RMSE, equation text for ggplot for all variables  ---------------

# Predict Ti using the model
predicted_Ti0 <- predict(plsr_model, ncomp = 4)
Ti_predicted <- predicted_Ti0[, "Ti_ICP_clr", 1]
Ti_actual <- icp_clr$Ti_ICP_clr
length(Ti_predicted) == length(Ti_actual)  # should be TRUE
# R-squared
r2_Ti <- cor(Ti_actual, Ti_predicted)^2
# RMSE
rmse_Ti <- sqrt(mean((Ti_actual - Ti_predicted)^2))
# RMSEP (same as RMSE if on test set)
rmsep_Ti <- sqrt(mean((Ti_actual - Ti_predicted)^2))
cat("Ti_R² =", r2_Ti, "\nTi_RMSE =", rmse_Ti, "\nTi_RMSEP =", rmse_Ti)

# Predict Ca using the model
predicted_Ca0 <- predict(plsr_model, ncomp = 4)
Ca_predicted <- predicted_Ca0[, "Ca_ICP_clr", 1]
Ca_actual <- icp_clr$Ca_ICP_clr
length(Ca_predicted) == length(Ca_actual)  # should be TRUE
# R-squared
r2_Ca <- cor(Ca_actual, Ca_predicted)^2
# RMSE
rmse_Ca <- sqrt(mean((Ca_actual - Ca_predicted)^2))
# RMSEP (same as RMSE if on test set)
rmsep_Ca <- sqrt(mean((Ca_actual - Ca_predicted)^2))
cat("Ca_R² =", r2_Ca, "\nCa_RMSE =", rmse_Ca, "\nCa_RMSEP =", rmsep_Ca)

# Predict Sr using the model
predicted_Sr0 <- predict(plsr_model, ncomp = 4)
Sr_predicted <- predicted_Sr0[, "Sr_ICP_clr", 1]
Sr_actual <- icp_clr$Sr_ICP_clr
length(Sr_predicted) == length(Sr_actual)  # should be TRUE
# R-squared
r2_Sr <- cor(Sr_actual, Sr_predicted)^2
# RMSE
rmse_Sr <- sqrt(mean((Sr_actual - Sr_predicted)^2))
# RMSEP (same as RMSE if on test set)
rmsep_Sr <- sqrt(mean((Sr_actual - Sr_predicted)^2))
cat("Sr_R² =", r2_Sr, "\nSr_RMSE =", rmse_Sr, "\nSr_RMSEP =", rmsep_Sr)

# Predict Zr using the model
predicted_Zr0 <- predict(plsr_model, ncomp = 4)
Zr_predicted <- predicted_Zr0[, "Zr_ICP_clr", 1]
Zr_actual <- icp_clr$Zr_ICP_clr
length(Zr_predicted) == length(Zr_actual)  # should be TRUE
# R-squared
r2_Zr <- cor(Zr_actual, Zr_predicted)^2
# RMSE
rmse_Zr <- sqrt(mean((Zr_actual - Zr_predicted)^2))
# RMSEP (same as RMSE if on test set)
rmsep_Zr <- sqrt(mean((Zr_actual - Zr_predicted)^2))
cat("Zr_R² =", r2_Zr, "\nZr_RMSE =", rmse_Zr, "\nZr_RMSEP =", rmse_Zr)

# Get coefficients and intercept from the model
coef_matrix <- coef(plsr_model, ncomp = 4)
intercept <- attr(coef_matrix, "constant")
# Extract Ti_ICP_clr column only as a named vector and retain keeping row names for equation
coef_Ti <- coef_matrix[, "Ti_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)
coef_Ca <- coef_matrix[, "Ca_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)
coef_Sr <- coef_matrix[, "Sr_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)
coef_Zr <- coef_matrix[, "Zr_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)

# Coefs fro clr are meaningless - can't convert clr coefficients into ppm directly.
# Need to apply the model to get predictions in clr
# Then back-transform predictions to composition (with clrInv()),
# Then multiply by total ICP to get back to ppm - see below and section 1.3 below for how this is done

# Show regression equation - this is meaningless for clr - conversion to ppm needs inversion of matrix first
cat("PLS clr Regression Equation:\n")
cat("Ti[ICP-MS_pred] clr = ",  paste0(round(coef_Ti, 4), "*", names(coef_Ti), collapse = " + "), "\n") #intercept removed as NULL 

# Generate text for ggplot for all model variables  ------------
#Ti
R2_text_Ti <- paste0("Ti_clrR^2: ", paste0(round(r2_Ti, 2)))
RMSE_text_Ti <- paste0("Ti_clrRMSE: ", paste0(round(rmse_Ti, 4)))
eq_text_Ti <- paste0("Ti_ICPMS[pred_clr]:", paste0(round(coef_Ti, 0), "*", names(coef_Ti), collapse = " + "))
#Ca
R2_text_Ca <- paste0("Ca_clrR^2: ", paste0(round(r2_Ca, 2)))
RMSE_text_Ca <- paste0("Ca_clrRMSE: ", paste0(round(rmse_Ca, 4)))
eq_text_Ca <- paste0("Ca_ICPMS[pred_clr]:", paste0(round(coef_Ca, 0), "*", names(coef_Ca), collapse = " + "))
#Sr
R2_text_Sr <- paste0("Sr_clrR^2: ", paste0(round(r2_Sr, 2)))
RMSE_text_Sr <- paste0("Sr_clrRMSE: ", paste0(round(rmse_Sr, 4)))
eq_text_Sr <- paste0("Sr_clrICPMS[pred_clr]:", paste0(round(coef_Sr, 0), "*", names(coef_Sr), collapse = " + "))
#Zr
R2_text_Zr <- paste0("Zr_clrR^2: ", paste0(round(r2_Zr, 2)))
RMSE_text_Zr <- paste0("Zr_clrRMSE: ", paste0(round(rmse_Zr, 4)))
eq_text_Zr <- paste0("Zr_ICPMS[pred_clr]:", paste0(round(coef_Zr, 0), "*", names(coef_Zr), collapse = " + "))

# Get CV predictions as clr, plot and then convert them to ppm values and plot
pred_clr <- predict(plsr_model, ncomp = 4, validation = "CV")[,,1]  # 3 components
colnames(pred_clr) <- paste0(colnames(pred_clr), "_pred")
pred_clr_tbl <- as_tibble(pred_clr)

# RMSE & 95% confidence intervals / error values as xrf_log_inc_new +/- --------

# Define R2 and RMSE values from training dataset (as above) 
RMSE_clr <- c(Ti_ICP_clr = rmse_Ti, Ca_ICP_clr = rmse_Ca, Sr_ICP_clr = rmse_Sr, Zr_ICP_clr = rmse_Zr)

# Simulate upper CIs for all elements in matrix by calculating p --------
pred_clr_upper <- pred_clr_tbl %>%
  map2_dfc(RMSE_clr, ~ .x + .y)

# Simulate lower CIs for all elements in matrix by calculating prediction bounds as pred +/-RMSE from PLS model 
pred_clr_lower <- pred_clr_tbl %>%
  map2_dfc(RMSE_clr, ~ .x - .y)

# Rename outputs
pred_clr_upper <- pred_clr_upper %>% rename_with(~ paste0(.x, "_upper"))
pred_clr_lower <- pred_clr_lower %>% rename_with(~ paste0(.x, "_lower"))

# Transform measured values (actual response) as clr
actual_clr <- clr(acomp(df_icp))  # same transformation used in model
colnames(actual_clr) <- paste0(colnames(actual_clr), "_clr_actual")
actual_clr_tbl <- as_tibble(actual_clr)

# RMSE and R² and CI summary table for each element  --------------------

summary_df_clr <- map2_dfr(actual_clr_tbl, pred_clr_tbl, ~ {
  tibble(
    R2 = cor(.x, .y)^2,
    RMSE = sqrt(mean((.x - .y)^2))
  )
}, .id = "Element")

# Rename based on actual column names
summary_df_clr$Element <- names(actual_clr)
print(summary_df_clr)

ci_list <- map2(actual_clr_tbl, pred_clr_tbl, ~ {
  residuals <- .x - .y
  rmse_val <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  tibble(
    CI_lower = quantile(residuals, 0.05, na.rm = TRUE),
    CI_upper = quantile(residuals, 0.95, na.rm = TRUE)
  )
})

# Combine into summary table
ci_df <- bind_rows(ci_list)
summary_df_clr_CI <- bind_cols(summary_df_clr, ci_df) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(summary_df_clr_CI)
write_csv(summary_df_clr_CI, "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/Summary_clr_R2_RMSE_CI.csv")

# Create plot df of PLSR measured vs predicted as clr  ------------------
plot_df_clr <- bind_cols(df, actual_clr, pred_clr, pred_clr_lower, pred_clr_upper)
plot_df_clr

# Ti training model - pred vs measured plot as clr ---------------------------------------
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ti_clr <- ggplot(plot_df_clr, aes(x = Ti_ICP_clr_actual, y = Ti_ICP_clr_pred, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(title = "PLSR: Measured vs Predicted Ti (ppm): clr model",
       x = "Measured Ti_ICP (clr)",
       y = "Predicted Ti_ICP (clr) [clr model]") +
  coord_equal() + # coord_equal()
  #theme_minimal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = -1, y = 1.5, label = RMSE_text_Ti, parse = TRUE, hjust = 0) +
  #annotate("text", x = -0.7, y = -0.9, label = eq_text, parse = TRUE, hjust = 0) + # clr is meaningless
  annotate("text", x = -1, y = 1.6, label = R2_text_Ti, parse = TRUE, hjust = 0)
print(PLS_final_Ti_clr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/Ti_plsr_clr_as_clr.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# Convert predicted and actual Ti clr values back to ppm values  --------

# Invert clr transform to get values in ppm
pred_acomp <- clrInv(pred_clr)
pred_upper_acomp <- clrInv(pred_clr_upper)
pred_lower_acomp <- clrInv(pred_clr_lower)
actual_acomp <- clrInv(actual_clr)
total_icp <- rowSums(df_icp) #Get total ICP per sample 

# Convert to tibble for tidyverse operations (optional)
pred_ppm <- pred_acomp %>%
  as_tibble() %>%
  mutate(total_icp = total_icp) %>%
  rowwise() %>%
  mutate(across(-total_icp, ~ .x * total_icp)) %>% #multiply each column by total_icp
  ungroup() %>%
  select(-total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm")) %>%
  mutate(across(everything(), ~ as.numeric(unlist(.x))))
pred_ppm

# Convert to ppm in tidyverse
pred_upper_ppm <- pred_upper_acomp %>%
  as_tibble() %>%
  mutate(mean_total_icp = mean_total_icp) %>%
  rowwise() %>%
  mutate(across(-mean_total_icp, ~ .x * mean_total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-mean_total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm")) %>%
  mutate(across(everything(), ~ as.numeric(unlist(.x))))
pred_upper_ppm

pred_lower_ppm <- pred_lower_acomp %>%
  as_tibble() %>%
  mutate(mean_total_icp = mean_total_icp) %>%
  rowwise() %>%
  mutate(across(-mean_total_icp, ~ .x * mean_total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-mean_total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm")) %>%
  mutate(across(everything(), ~ as.numeric(unlist(.x))))
pred_lower_ppm

# Invert actual clr back to df_icp values - check they are the same
actual_acomp <- clrInv(actual_clr)  # proportions
# Convert to tibble for tidyverse operations (optional)
actual_ppm <- actual_acomp %>%
  as_tibble() %>%
  mutate(total_icp = total_icp) %>%
  rowwise() %>%
  mutate(across(-total_icp, ~ .x * total_icp)) %>%
  ungroup() %>%
  select(-total_icp) %>% 
  rename_with(~ paste0(.x, "_actual_ppm")) %>%
  mutate(across(everything(), ~ as.numeric(unlist(.x))))
actual_ppm

# RMSE and R² and CI summary table for each element  --------------------

summary_df_ppm <- map2_dfr(actual_ppm, pred_ppm, ~ {
  tibble(
    R2 = cor(.x, .y)^2,
    RMSE = sqrt(mean((.x - .y)^2)),
  )
}, .id = "Element")

# Rename based on actual column names
summary_df_clr$Element <- names(actual_ppm)
print(summary_df_ppm)

ci_list <- map2(actual_ppm, pred_ppm, ~ {
  residuals <- .x - .y
  rmse_val <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  tibble(
    CI_lower = quantile(residuals, 0.05, na.rm = TRUE),
    CI_upper = quantile(residuals, 0.95, na.rm = TRUE)
  )
})

# Combine into summary table
ci_df <- bind_rows(ci_list)
summary_df_ppm_CI <- bind_cols(summary_df_ppm, ci_df) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(summary_df_ppm_CI)
write_csv(summary_df_ppm_CI, "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/Summary_ppm_R2_RMSE_CI.csv")

# Set up Plot df for PLSR measured vs predicted as ppm 
plot_df_ppm <- bind_cols(df, actual_ppm, pred_ppm, pred_upper_ppm, pred_lower_ppm)
plot_df_ppm

# Ti training model - pred vs measured plot as ppm ---------------------------------------
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ti_ppm <- ggplot(plot_df_ppm, aes(x = Ti_ICP_actual_ppm, y = Ti_ICP_clr_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(title = "PLSR: Measured vs Predicted Ti (ppm): clr model",
       x = "Measured Ti_ICP (ppm)",
       y = "Predicted Ti_ICP (ppm) [clr model]") +
  coord_equal(xlim = c(0, 40000), ylim = c(0, 40000)) + # coord_equal()
  #theme_minimal()
  #coord_equal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0, y = 38000, label = RMSE_text_Ti, parse = TRUE, hjust = 0) +
  #annotate("text", x = 4.5, y = 28000, label = eq_text, parse = TRUE, hjust = 0) + # clr eq is meaningless
  annotate("text", x = 0, y = 40000, label = R2_text_Ti, parse = TRUE, hjust = 0)
  print(PLS_final_Ti_ppm)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/Ti_plsr_clr_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# ------------------------------------------------------------------------------

# 1.3 Predicting Ti as ppm for a new xrf dataset -------------------------------

# Import ACE Ti XRF-CS log_inc data and convert to ppm with +/- RMSE error
ACE_xrf_cps <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_cps.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>% 
  select(Location:MSE, all_of(main_elements_xrf), Total_scatter, inc_coh, coh_inc) %>%
  filter(qc == "TRUE") %>% 
  filter(validity == "1") %>% 
  mutate(across(all_of(main_elements_xrf), ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_cps

ACE_xrf_predict_PLS <- ACE_xrf_cps %>%
  select(Site, depth, SH20_age, Ti, Ca, Sr, Zr)

# Convert tibble to data frame for base R
X_new <- ACE_xrf_predict_PLS %>% 
  select (Ti, Ca, Sr, Zr)
df_xrf_new <- as.matrix(X_new)

# Transform to clr using compositions package and add _clr to column titles
xrf_clr_new <- clr(acomp(df_xrf_new))
colnames(xrf_clr_new) <- paste0(colnames(xrf_clr_new), "_clr")
xrf_clr_new[xrf_clr_new == 0] <- NA
#xrf_clr_new  <- as.data.frame(xrf_clr_new0) %>%
#  mutate(across(all_of(c("Ti_clr", "Ca_clr", "Sr_clr")), ~ ifelse(. == 0, min(.[. != 0]) / 2, .)))
# replace zeros with half minimum value to allow clr model to work
# Recommended procedure from Bertrand et al. (submitted) - retains dataframe structure

# Input new values as xrf_clr_new
xrf_clr_new

# Get predicted xrf clr values
pred_xrf_clr <- predict(plsr_model, newdata = xrf_clr_new, ncomp = 4, validation = "CV")[,,1]  # 3 components

# Invert predicted clr back to proportions using clrInv (to est. ppm values) 
pred_xrf_acomp <- clrInv(pred_xrf_clr)  # proportions

# Calculate average total from the training ICPMS dataset (n=265) with known total row sums
mean_total_icp <- mean(rowSums(df_icp))

# Estimate a representative total (e.g., the mean row sum) clr conversion and apply that to the new predictions
# Convert to ppm in tidyverse - only Ti is valid as a predicted estimate here
pred_xrf_ppm <- pred_xrf_acomp %>%
  as_tibble() %>%
  mutate(mean_total_icp = mean_total_icp) %>%
  rowwise() %>%
  mutate(across(-mean_total_icp, ~ .x * mean_total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-mean_total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm"))
pred_xrf_ppm

# RMSE & 95% confidence intervals / error values as xrf_log_inc_new +/- --------

# Get predicted xrf clr values for Ti, Ca, Sr as tibble
pred_xrf_clr <- predict(plsr_model, newdata = xrf_clr_new, ncomp = 4, validation = "CV")[,,1]  # 3 components
pred_xrf_clr_tbl <- as_tibble(pred_xrf_clr)

# Define R2 and RMSE values from training dataset (as above) 
RMSE_clr <- c(Ti_ICP_clr = rmse_Ti, Ca_ICP_clr = rmse_Ca, Sr_ICP_clr = rmse_Sr, Zr_ICP_clr = rmse_Zr)

# Simulate upper CIs for all elements in matrix by calculating p --------
pred_xrf_clr_upper <- pred_xrf_clr_tbl %>%
  map2_dfc(RMSE_clr, ~ .x + .y)

# Simulate lower CIs for all elements in matrix by calculating prediction bounds as pred +/-RMSE from PLS model 
pred_xrf_clr_lower <- pred_xrf_clr_tbl %>%
  map2_dfc(RMSE_clr, ~ .x - .y)

# Rename outputs
pred_xrf_clr_upper <- pred_xrf_clr_upper %>% rename_with(~ paste0(.x, "_upper"))
pred_xrf_clr_lower <- pred_xrf_clr_lower %>% rename_with(~ paste0(.x, "_lower"))

# Invert clr transform to get values in ppm
pred_xrf_upper_acomp <- clrInv(pred_xrf_clr_upper)
pred_xrf_lower_acomp <- clrInv(pred_xrf_clr_lower)

# Convert to ppm in tidyverse
pred_xrf_upper_ppm <- pred_xrf_upper_acomp %>%
  as_tibble() %>%
  mutate(mean_total_icp = mean_total_icp) %>%
  rowwise() %>%
  mutate(across(-mean_total_icp, ~ .x * mean_total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-mean_total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm"))
pred_xrf_upper_ppm

pred_xrf_lower_ppm <- pred_xrf_lower_acomp %>%
  as_tibble() %>%
  mutate(mean_total_icp = mean_total_icp) %>%
  rowwise() %>%
  mutate(across(-mean_total_icp, ~ .x * mean_total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-mean_total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm"))
pred_xrf_lower_ppm

pred_xrf_ppm_df <- bind_cols(pred_xrf_ppm, pred_xrf_upper_ppm, pred_xrf_lower_ppm) %>% 
  relocate(Ti_ICP_clr_lower_ppm, .after = Ti_ICP_clr_ppm) %>% 
  relocate(Ti_ICP_clr_upper_ppm, .after = Ti_ICP_clr_lower_ppm) %>% 
  relocate(Ca_ICP_clr_lower_ppm, .after = Ca_ICP_clr_ppm) %>% 
  relocate(Ca_ICP_clr_upper_ppm, .after = Ca_ICP_clr_lower_ppm) %>% 
  relocate(Sr_ICP_clr_lower_ppm, .after = Sr_ICP_clr_ppm) %>% 
  relocate(Sr_ICP_clr_upper_ppm, .after = Sr_ICP_clr_lower_ppm) %>% 
  relocate(Zr_ICP_clr_lower_ppm, .after = Zr_ICP_clr_ppm) %>% 
  relocate(Zr_ICP_clr_upper_ppm, .after = Zr_ICP_clr_lower_ppm) %>% 
  rename(Ti_low_RMSE = Ti_ICP_clr_lower_ppm) %>% 
  rename(Ca_low_RMSE = Ca_ICP_clr_lower_ppm) %>% 
  rename(Sr_low_RMSE = Sr_ICP_clr_lower_ppm) %>% 
  rename(Zr_low_RMSE = Zr_ICP_clr_lower_ppm) %>% 
  rename(Ti_up_RMSE = Ti_ICP_clr_upper_ppm) %>% 
  rename(Ca_up_RMSE = Ca_ICP_clr_upper_ppm) %>% 
  rename(Sr_up_RMSE = Sr_ICP_clr_upper_ppm) %>% 
  rename(Zr_up_RMSE = Zr_ICP_clr_upper_ppm)
# Make into tibble for export and plotting & write to file
ACE_ppm_PLS_pred_clr <- bind_cols(ACE_xrf_predict_PLS, pred_xrf_ppm_df)
write.csv(ACE_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv", row.names = FALSE)

# Extract HER42PB data for Figure 4
HER42PB_ppm_PLS_pred_clr <- ACE_ppm_PLS_pred_clr %>% 
  filter(Site == "HER42PB")
HER42PB_ppm_PLS_pred_clr
write.csv(HER42PB_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/clr/HER42PB_PLS_ppm_pred_clr.csv", row.names = FALSE)

# END --------------------------------------------------------------------------
