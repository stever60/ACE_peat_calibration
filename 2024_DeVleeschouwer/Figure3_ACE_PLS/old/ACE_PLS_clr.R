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
              'compositions', 'wesanderson', 'viridis', 'itraxR','pls', 'caret')
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

# If the dataset is small, use eg k-fold cross-validation to evaluate the model's performance without relying on a single test set. 
# Calibrate Titanium (Ti) measured by ICP-MS against other elements measured by ITRAX using pls package.

# Define predictor ITRAX variables & ICPMS response variable for calibration model tests
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
main_elements_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP","Fe_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP","Zr_ICP")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")

# Load required libraries
library(compositions)  # clr, acomp, clrInv
library(pls)           # plsr
library(ggplot2)
library(dplyr)

# Input XRF and ICPMS ---

# Example: assume df_xrf and df_icp are your original data frames

df_key_elements_xrf <- df[, c(key_elements_xrf)]
df_key_elements_icp <- df[, c(key_elements_icp)]

df_xrf <- df[, c("Ti", "Ca", "Sr")]
df_icp <- df[, c("Ti_ICP", "Ca_ICP", "Sr_ICP")]

# CLR-transform predictors and response ---

# Transform to compositions and add _clr to column titles
xrf_clr <- clr(acomp(df_xrf))
colnames(xrf_clr) <- paste0(colnames(xrf_clr), "_clr")
icp_clr <- clr(acomp(df_icp))
colnames(icp_clr) <- paste0(colnames(icp_clr), "_clr")

# Fit CV-PLSR (cross-validated) test model to determine optimal number of comps ---

plsr_model <- plsr(icp_clr ~ xrf_clr, ncomp = 3, validation = "CV", jackknife = TRUE, scale = TRUE, center = TRUE)

# Predict Ti using the model
predicted_Ti0 <- predict(plsr_model, ncomp = 3)
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
predicted_Ca0 <- predict(plsr_model, ncomp = 3)
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
predicted_Sr0 <- predict(plsr_model, ncomp = 3)
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

# Get coefficients and intercept from the model
coef_matrix <- coef(plsr_model, ncomp = 3)
intercept <- attr(coef_matrix, "constant")
# Extract Ti_ICP_clr column only as a named vector and retain keeping row names for equation
coef_Ti <- coef_matrix[, "Ti_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)

# Coefs here are a bit pointless as you can't convert clr coefficients into ppm directly.
# Need to apply the model to get predictions in clr,
# Back-transform predictions to composition (with clrInv()),
# Multiply by total ICP to get back to ppm.

# Show regression equation
cat("PLS clr Regression Equation:\n")
cat("Ti[ICP-MS_pred] clr = ",  paste0(round(coef_Ti, 4), "*", names(coef_Ti), collapse = " + "), "\n") #intercept removed as NULL 

# Generate equation text for ggplot
R2_text <- paste0("R^2: ", paste0(round(r2_Ti, 2)))
RMSE_text <- paste0("RMSE: ", paste0(round(rmse_Ti, 4)))
eq_text <- paste0("Ti_clr[ICP-MS_predcited]:", paste0(round(coef_Ti, 0), "*", names(coef_Ti), collapse = " + "))

# Get CV predictions as clr, plot and then convert them to ppm values and plot
pred_clr <- predict(plsr_model, ncomp = 3, validation = "CV")[,,1]  # 3 components
colnames(pred_clr) <- paste0(colnames(pred_clr), "_pred")

# Transform measured values (actual response) as clr
actual_clr <- clr(acomp(df_icp))  # same transformation used in model
colnames(actual_clr) <- paste0(colnames(actual_clr), "_clr_actual")

# Plot PLSR measured vs predicted as clr 
plot_df_clr <- bind_cols(df, actual_clr, pred_clr)
plot_df_clr

# Plot PLSR measured vs predicted as ppm 
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_clr <- ggplot(plot_df_clr, aes(x = Ti_ICP_clr_actual, y = Ti_ICP_clr_pred, color = Site)) +
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
  coord_equal() + # coord_equal()
  #theme_minimal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = -1, y = 1.5, label = RMSE_text, parse = TRUE, hjust = 0) + # add equation to graph
  #annotate("text", x = -0.7, y = -0.9, label = eq_text, parse = TRUE, hjust = 0) + # add equation to graph
  annotate("text", x = -1, y = 1.6, label = R2_text, parse = TRUE, hjust = 0) # add equation to graph
print(PLS_final_clr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/clr/Ti_plsr_clr_as_clr.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# You can't convert clr coefficients into ppm directly.
# Instead:Apply the model to get predictions in clr,
# Back-transform predictions to composition (with clrInv()),
# Multiply by total ICP to get back to ppm.

# Convert predicted and actual Ti clr values back to ppm values using clrInv
pred_acomp <- clrInv(pred_clr)
actual_acomp <- clrInv(actual_clr)
total_icp <- rowSums(df_icp) #Get total ICP per sample 

# Invert pred clr to ppm values
pred_acomp <- clrInv(pred_clr)  # proportions
# Convert to tibble for tidyverse operations (optional)
pred_ppm <- pred_acomp %>%
  as_tibble() %>%
  mutate(total_icp = total_icp) %>%
  rowwise() %>%
  mutate(across(-total_icp, ~ .x * total_icp)) %>% #multiply each column by total_icp
  ungroup() %>%
  select(-total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm"))
pred_ppm

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
  rename_with(~ paste0(.x, "_actual_ppm"))
actual_ppm

# Plot PLSR measured vs predicted as ppm 
plot_df_ppm <- bind_cols(df, actual_ppm, pred_ppm)
plot_df_ppm

library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final <- ggplot(plot_df_ppm, aes(x = Ti_ICP_actual_ppm, y = Ti_ICP_clr_pred_ppm, color = Site)) +
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
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0, y = 38000, label = RMSE_text, parse = TRUE, hjust = 0) + # add equation to graph
  #annotate("text", x = 4.5, y = 28000, label = eq_text, parse = TRUE, hjust = 0) + # add equation to graph
  annotate("text", x = 0, y = 40000, label = R2_text, parse = TRUE, hjust = 0) # add equation to graph
  print(PLS_final)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/clr/Ti_plsr_clr_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Predict Ti as ppm using pls_final model and Ln(Ti/inc) ITRAX dataset ---------

# Import ACE Ti XRF-CS log_inc data and convert to ppm with +/- RMSE error
ACE_xrf_cps <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_cps.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>% 
  select(Location:MSE, all_of(main_elements_xrf), Total_scatter, inc_coh, coh_inc) %>%
  filter(qc == "TRUE") %>% 
  filter(validity == "1") %>% 
  mutate(across(all_of(main_elements_xrf), ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_cps

ACE_xrf_predict_Ti_PLS <- ACE_xrf_cps %>%
  select(Site, depth, SH20_age, Ti, Ca, Sr)

# Convert tibble to data frame for base R
X_new <- ACE_xrf_predict_Ti_PLS %>% 
  select (Ti, Ca, Sr) %>% 
  as.data.frame(ACE_xrf_predict_Ti_PLS)

# Transform to clr using compositions package and add _clr to column titles
xrf_clr_new <- clr(acomp(X_new))
colnames(xrf_clr_new) <- paste0(colnames(xrf_clr_new), "_clr")
xrf_clr_new[xrf_clr_new == 0] <- NA
#xrf_clr_new  <- as.data.frame(xrf_clr_new0) %>%
#  mutate(across(all_of(c("Ti_clr", "Ca_clr", "Sr_clr")), ~ ifelse(. == 0, min(.[. != 0]) / 2, .)))
# replace zeros with half minimum value to allow clr model to work
# Recommended procedure from Bertrand et al. (submitted) - retains dataframe structure

# Input new values as xrf_clr_new
xrf_clr_new

# Get predicted xrf clr values
pred_xrf_clr <- predict(plsr_model, newdata = xrf_clr_new, ncomp = 3, validation = "CV")[,,1]  # 3 components

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

#  RMSE error values as xrf_clr_new + RMSE of clr model

# Get predicted xrf clr values for Ti, Ca, Sr
pred_xrf_clr <- predict(plsr_model, newdata = xrf_clr_new, ncomp = 3, validation = "CV")[,,1]  # 3 components

# Define RMSE avlues rom training dataset (as above)
RMSE_clr <- c(Ti = rmse_Ti, Ca = rmse_Ca, Sr = rmse_Sr)

pred_clr_tbl <- as_tibble(pred_xrf_clr)

# Compute upper bound
pred_clr_upper <- pred_clr_tbl %>%
  map2_dfc(RMSE_clr, ~ .x + .y)

# Compute lower bound
pred_clr_lower <- pred_clr_tbl %>%
  map2_dfc(RMSE_clr, ~ .x - .y)

# Rename outputs
pred_clr_upper <- pred_clr_upper %>% rename_with(~ paste0(.x, "_upper"))
pred_clr_lower <- pred_clr_lower %>% rename_with(~ paste0(.x, "_lower"))

# Invert clr transform to get values in ppm
pred_upper_acomp <- clrInv(pred_clr_upper)
pred_lower_acomp <- clrInv(pred_clr_lower)

# Convert to ppm in tidyverse - only Ti is valid as a predicted estimate here
pred_xrf_upper_ppm <- pred_upper_acomp %>%
  as_tibble() %>%
  mutate(mean_total_icp = mean_total_icp) %>%
  rowwise() %>%
  mutate(across(-mean_total_icp, ~ .x * mean_total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-mean_total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm"))
pred_xrf_upper_ppm

pred_xrf_lower_ppm <- pred_lower_acomp %>%
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
  relocate(Sr_ICP_clr_upper_ppm, .after = Sr_ICP_clr_lower_ppm)
pred_xrf_ppm_df

# Make into tibble for export and plotting & write to file
Ti_ppm_PLS_predicted <- bind_cols(ACE_xrf_predict_Ti_PLS, pred_xrf_ppm_df)
write.csv(Ti_ppm_PLS_predicted,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/clr/ACE_Ti_ppm_PLS_predicted.csv", row.names = FALSE)

# Extract HER42PB data for Figure 4
HER42PB_Ti_ppm_PLS_predicted <- Ti_ppm_PLS_predicted %>% 
  filter(Site == "HER42PB")
HER42PB_Ti_ppm_PLS_predicted
write.csv(HER42PB_Ti_ppm_PLS_predicted,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/clr/HER42PB_Ti_ppm_PLS_predicted.csv", row.names = FALSE)

# END --------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 2. PLS Train-Test Split tidyverse method  --------------------------

# If the primary goal is to build a highly accurate predictive model, a larger training set might be preferred
# If the focus is on evaluating the model's generalizability, a larger testing set might be more appropriate
# As above but with 70% training vs 30% test dataset 

# Define df_train as training dataset and elements to use for Ti_ICP response model
df_Site <- df[, (names(df) %in% c("Site", "Ti_ICP", key_elements_xrf))]

df_split <- df_Site %>%
  group_by(Site) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

train_df <- df_split %>%
  group_by(Site) %>%
  slice_sample(prop = 0.7) %>%
  ungroup()

test_df <- anti_join(df_split, train_df, by = c("Site", "row_id"))

# Prepare matrices
X_train <- train_df %>% select(all_of(key_elements_xrf)) %>% as.matrix()
Y_train <- train_df$Ti_ICP

# Fit PLSR model
pls_model <- plsr(Y_train ~ X_train, scale = TRUE, validation = "CV", ncomp = 6)

# Optimal number of components
validationplot(pls_model, val.type = "MSEP")
opt_comp_train_LOO <- which.min(RMSEP(pls_model)$val[1, , -1])
cat("Optimal number of components:", opt_comp_train_LOO, "\n")

# Predict on Test Set
X_test <- test_df %>% select(all_of(key_elements_xrf)) %>% as.matrix()
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
ggplot(train_plot_df, aes(x = Comp1, y = Comp2, color = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PLSR Score Plot (Training Data)", x = "PLS Component 1", y = "PLS Component 2") +
  theme_minimal()

# Compute R² and RMSE on test set
SS_res <- sum((Y_test - Y_pred)^2)
SS_tot <- sum((Y_test - mean(Y_test))^2)
R2_test <- 1 - SS_res / SS_tot
RMSE_test <- sqrt(mean((Y_test - Y_pred)^2))
cat("Optimal Components:", opt_comp, "\n")
cat("Test R²:", round(R2_test, 4), "\n")
cat("Test RMSE:", round(RMSE_test, 4), "\n")

# Extract coefficients for optimal number of components
coefs <- coef(pls_model, ncomp = opt_comp)
coefs_vec <- as.vector(coefs)
# Match predictor names and print to console
names(coefs_vec) <- colnames(X_train)
print(round(coefs_vec, 4))

# Predicted vs Measured Plot
ggplot(Ti_ICP_test_df, aes(x = Ti_ICP, y = Ti_ICP_pred, color = Site)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Titanium (Ti) (30% Test Set))",
    x = "Measured Ti (ICP-MS)",
    y = "Predicted Ti (PLS Model)"
  ) +
  theme_minimal()

# Regression Equation

# Round coefficients
coefs_rounded <- round(coefs_vec, 4)

# Build equation string
equation <- paste0("Ti_ICP_predicted = ",
                   paste(paste0(coefs_rounded, " * ", names(coefs_rounded)), collapse = " + "))
cat("PLSR Regression Equation (", opt_comp, " components):\n", equation, "\n")

# Equations with intercept are for unscaled only - this doesnt work here because intecept is centred and returns NULL value
# The plsr() function centers and scales data by default, so coefficients are for scaled predictors and centered response. 
# If you want to reconstruct the full unscaled equation including the intercept, you can use the model components directly:

# Get scaled means
#X_means <- attr(X_train, "scaled:center")
#X_sds <- attr(X_train, "scaled:scale")
#Y_mean <- mean(Y_train)

# Adjust coefficients back to original scale
#coefs_raw <- coefs_vec / X_sds
#intercept <- Y_mean - sum(coefs_raw * X_means)

# Rebuild full equation with intercept
#equation_raw <- paste0("Ti_ICP_predicted = ", round(intercept, 4), " + ",
#                       paste(paste0(round(coefs_raw, 4), " * ", names(coefs_raw)), collapse = " + "))
#cat("Final regression equation (unscaled):\n", equation_raw, "\n")

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

# Calculate R2 and RMSE values
r2 <- cor(actuals, predictions)^2
rmse <- sqrt(mean((actuals - predictions)^2))
cat("Test R² =", round(r2, 2), "\nTest RMSE =", round(rmse, 4), "\n")

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

# old code\

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

# define coefficients and RMSE values to use in equations - see ACE.PLS.R file outputs for details
coefs_1 = 0.8160493 # Ti
coefs_2 = -0.5945589 # Ca
coefs_3 = 0.3389237 # Sr
rmse = 0.6200723 # RMSE value

ACE_xrf_conversion_Ti_PLS <- ACE_xrf_log_inc %>%
  select(Site, depth, SH20_age, Ti, Ca, Sr) %>%
  filter (Site == "HER42PB") %>% 
  mutate(Ti_convert_PLS = (coefs_1**Ti)+ (coefs_2*Ca) + (coefs_3*Sr)) %>% #
  mutate(Ti_ppm_PLS = exp(Ti_convert_PLS)) %>% 
  mutate(Ti_upper_PLS = Ti_convert_PLS+rmse) %>% # RMSE = 0.6201
  mutate(Ti_upper_RMSE_PLS = exp(Ti_upper_PLS)) %>% 
  mutate(Ti_lower_PLS = Ti_convert_PLS-rmse) %>% 
  mutate(Ti_lower_RMSE_PLS = exp(Ti_lower_PLS)) %>% 
  select(Site, depth, SH20_age, Ti, Ca, Sr, Ti_convert_PLS, Ti_ppm_PLS, Ti_lower_RMSE_PLS, Ti_upper_RMSE_PLS)
ACE_xrf_conversion_Ti_PLS
write.csv(ACE_xrf_conversion_Ti_PLS,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/cps/ACE_plsr_xrf_Ti_converted_PLS.csv", row.names = FALSE)

