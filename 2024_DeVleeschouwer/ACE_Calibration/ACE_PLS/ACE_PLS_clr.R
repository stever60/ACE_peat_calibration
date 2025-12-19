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
              'compositions', 'boot', 'broom', 'glue', "ggpmisc",
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

icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", 
                      "Ni_ICP", "Cu_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", 
                         "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# Define predictor ITRAX variables & ICPMS response variable for calibration model tests
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
main_elements_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP","Fe_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP","Zr_ICP")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")
final_elements_xrf <- c("Ti", "Ca", "Sr", "Zr")
final_elements_icp <- c("Ti_ICP", "Ca_ICP", "Sr_ICP", "Zr_ICP")

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

# Load and Prepare clr matched data for based on all 12 xrf_icp_Elements_min elements for validatiion Section 1.1
ACE_matched_xrf_icp_clr <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_clr.csv") 
is.na(ACE_matched_xrf_icp_clr)<-sapply(ACE_matched_xrf_icp_clr, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_clr

# Load and Prepare cps matched data for Sectuon 1.2 / 1.3 - Training and new data predictions 
ACE_matched_xrf_icp_cps <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") 
is.na(ACE_matched_xrf_icp_cps)<-sapply(ACE_matched_xrf_icp_cps, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_cps

#-------------------------------------------------------------------------------

# 1. Partial Least Squares regression (plsr) validation/optimal components Ti --

# ------------------------------------------------------------------------------

# Run once only
# Select dataset and convert to base R dataframe for pls package ---------------
df_clr <- ACE_matched_xrf_icp_cps %>%
  select(c(Site, depth, SH20_mean_age, SH20_mean_95CI, all_of(xrf_icp_elements1))) %>%
  drop_na()
df_clr
df_val <- as.data.frame(df_clr) # Convert tibble to data frame for base R
# use original cps dataset and create a need clr matrix based on [signifcant elements] only!
# from log_inc dataset, only "Ti", "Ca", "Sr", "Zr" had p<0.001 in CV & LOO validation tests

# Validation model - assessing signifcant elements & number of com --------

# If the dataset is small, use eg k-fold cross-validation to evaluate the model's performance without relying on a single test set. 
# Calibrate Titanium (Ti) measured by ICP-MS against other elements measured by ITRAX using pls package.

# Define predictor ITRAX variables & ICPMS response variable for calibration model tests
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
main_elements_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP","Fe_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP","Zr_ICP")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")

# Validation models - single element test - XRF predictors (X) and ICP-MS response
X_test <- df_val[, key_elements_xrf]  # ITRAX elements
Y_test <- df_val$Ti_ICP #define response ICP-MS variable for calibration



# LOO-PLSR (leave-one-out) test model to determine optimal number of comps -----

# Perform PLS regression with LOO - leave one out validation - useful for unbiased prediction error assessment
pls_train_model_LOO <- plsr(Y_test ~ ., data = data.frame(X_test, Y_test), validation = "LOO", 
                            jackknife = TRUE, scale = TRUE, center = TRUE) 
# Use scale = TRUE to perform as Z scores - do not use for calibration

# Get optimal number of components (lowest LOO error)
opt_comp_LOO <- which.min(RMSEP(pls_train_model_LOO)$val[1, , -1])
cat("Optimal number of components LOO:", opt_comp_LOO, "\n")

# Summary of LOO model output and jackknife t tests of regression coefficients & std errors
summary(pls_train_model_LOO)
jack.test(pls_train_model_LOO, ncomp = opt_comp_LOO, use.mean = TRUE)
RMSEP(pls_train_model_LOO) # compare RMSEP and RMSE

# Save summary output to a text file
capture.output(summary(pls_train_model_LOO), file = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PLS_model_LOO_summary_Ti.txt")

# Plot RMSEP to screen to choose number of components - where first becomes flat
plot(RMSEP(pls_train_model_LOO), legendpos = "topright")

# Open a PDF in device
pdf("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PLS_rmsep_plot_LOO.pdf", width = 4, height = 4)

# Create the plot
plot(RMSEP(pls_train_model_LOO), legendpos = "topright")

# LOO - RMSE & Jackknifing - all - write console to txt file -------------------

# Define predictor and response variables
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")

# Loop over each key ICP element
results_LOO <- list()
for (icp_element in key_elements_icp) {
  
  # Set Y_test dynamically based on current ICPMS element
  Y_test <- df_clr[[icp_element]]
  
  # Create combined data frame of predictors and response
  data_model <- data.frame(df_clr[, key_elements_xrf], Y_test = Y_test)
  
  # Run PLSR with leave-one-out validation
  pls_model <- plsr(Y_test ~ ., data = data_model, 
                    validation = "LOO", 
                    jackknife = TRUE,
                    scale = TRUE, center = TRUE)
  
  # Get RMSEP and optimal number of components
  rmsep_vals <- RMSEP(pls_model)$val[1, , -1]
  opt_ncomp <- which.min(rmsep_vals)
  
  # Get jackknife test results
  jack_res <- jack.test(pls_model, ncomp = opt_ncomp, use.mean = TRUE)
  
  # Store results in a list
  results_LOO[[icp_element]] <- list(
    RMSEP = rmsep_vals,
    optimal_components = opt_ncomp,
    jackknife_test = jack_res,
    model_summary = summary(pls_model)
  )
  
  # Optionally print status
  cat("Finished validation for", icp_element, 
      " | Optimal components:", opt_ncomp, "\n")
}

# View RMSEP values for Fe_ICP
results_LOO$Ti_ICP$RMSEP
results_LOO$Ca_ICP$RMSEP
results_LOO$Mn_ICP$RMSEP
results_LOO$Fe_ICP$RMSEP
results_LOO$Sr_ICP$RMSEP
results_LOO$Zr_ICP$RMSEP
# View jackknife test for Mn_ICP
results_LOO$Ti_ICP$jackknife_test
results_LOO$Ca_ICP$jackknife_test
results_LOO$Mn_ICP$jackknife_test
results_LOO$Fe_ICP$jackknife_test
results_LOO$Sr_ICP$jackknife_test
results_LOO$Zr_ICP$jackknife_test
# Print LOO results_LOO from console output to a txt file -------------------
sink("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/RMSEP_Jackknife_all_LOO.txt")
cat("RMSEP values for ICP elements:\n\n")
cat("Ti_ICP RMSEP:\n")
print(results_LOO$Ti_ICP$RMSEP)
cat("\nCa_ICP RMSEP:\n")
print(results_LOO$Ca_ICP$RMSEP)
cat("\nMn_ICP RMSEP:\n")
print(results_LOO$Mn_ICP$RMSEP)
cat("\nFe_ICP RMSEP:\n")
print(results_LOO$Fe_ICP$RMSEP)
cat("\nSr_ICP RMSEP:\n")
print(results_LOO$Sr_ICP$RMSEP)
cat("\nZr_ICP RMSEP:\n")
print(results_LOO$Zr_ICP$RMSEP)
cat("\n\nJackknife test results_LOO for ICP elements:\n\n")
cat("Ti_ICP jackknife_test:\n")
print(results_LOO$Ti_ICP$jackknife_test)
cat("\nCa_ICP jackknife_test:\n")
print(results_LOO$Ca_ICP$jackknife_test)
cat("\nMn_ICP jackknife_test:\n")
print(results_LOO$Mn_ICP$jackknife_test)
cat("\nFe_ICP jackknife_test:\n")
print(results_LOO$Fe_ICP$jackknife_test)
cat("\nSr_ICP jackknife_test:\n")
print(results_LOO$Sr_ICP$jackknife_test)
cat("\nZr_ICP jackknife_test:\n")
print(results_LOO$Zr_ICP$jackknife_test)
sink()





# CV-PLSR (cross-validated) test model to determine optimal number of comps ----

# Perform PLS regression with 10-fold cross-validation for model tuning - useful for small datasets
pls_train_model_CV <- plsr(Y_test ~ ., data = data.frame(X_test, Y_test), validation = "CV", 
                           jackknife = TRUE, scale = TRUE, center = TRUE)

# Get optimal number of components (lowest CV error)
opt_comp_CV <- which.min(RMSEP(pls_train_model_CV)$val[1, , -1])
cat("Optimal number of components CV:", opt_comp_CV, "\n")
RMSEP(pls_train_model_CV) # compare RMSEP and RMSE

# Summary of  CV model output and jackknife t tests of regression coefficients & std errors
summary(pls_train_model_CV)
jack.test(pls_train_model_CV, ncomp = opt_comp_CV, use.mean = TRUE)

# Save summary output to a text file
capture.output(summary(pls_train_model_CV), file = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PLS_model_CV_summary_Ti.txt")

# Plot RMSEP to screen to choose number of components - where first becomes flat
plot(RMSEP(pls_train_model_CV), legendpos = "topright")

# Open a PDF in device
pdf("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PLS_rmsep_plot_CV.pdf", width = 4, height = 4)

# Create the plot
plot(RMSEP(pls_train_model_CV), legendpos = "topright")
# CV - RMSE & Jackknifing - all - write console to txt file -------------------

# Define predictor and response variables
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")

# CVp over each key ICP element
results_CV <- list()
for (icp_element in key_elements_icp) {
  
  # Set Y_test dynamically based on current ICPMS element
  Y_test <- df_clr[[icp_element]]
  
  # Create combined data frame of predictors and response
  data_model <- data.frame(df_clr[, key_elements_xrf], Y_test = Y_test)
  
  # Run PLSR with leave-one-out validation
  pls_model <- plsr(Y_test ~ ., data = data_model, 
                    validation = "CV", 
                    jackknife = TRUE,
                    scale = TRUE, center = TRUE)
  
  # Get RMSEP and optimal number of components
  rmsep_vals <- RMSEP(pls_model)$val[1, , -1]
  opt_ncomp <- which.min(rmsep_vals)
  
  # Get jackknife test results
  jack_res <- jack.test(pls_model, ncomp = opt_ncomp, use.mean = TRUE)
  
  # Store results in a list
  results_CV[[icp_element]] <- list(
    RMSEP = rmsep_vals,
    optimal_components = opt_ncomp,
    jackknife_test = jack_res,
    model_summary = summary(pls_model)
  )
  
  # Optionally print status
  cat("Finished validation for", icp_element, 
      " | Optimal components:", opt_ncomp, "\n")
}

# View RMSEP values for Fe_ICP
results_CV$Ti_ICP$RMSEP
results_CV$Ca_ICP$RMSEP
results_CV$Mn_ICP$RMSEP
results_CV$Fe_ICP$RMSEP
results_CV$Sr_ICP$RMSEP
results_CV$Zr_ICP$RMSEP
# View jackknife test for Mn_ICP
results_CV$Ti_ICP$jackknife_test
results_CV$Ca_ICP$jackknife_test
results_CV$Mn_ICP$jackknife_test
results_CV$Fe_ICP$jackknife_test
results_CV$Sr_ICP$jackknife_test
results_CV$Zr_ICP$jackknife_test
# Print CV results_CV from console output to a txt file -------------------
sink("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/RMSEP_Jackknife_all_CV.txt")
cat("RMSEP values for ICP elements:\n\n")
cat("Ti_ICP RMSEP:\n")
print(results_CV$Ti_ICP$RMSEP)
cat("\nCa_ICP RMSEP:\n")
print(results_CV$Ca_ICP$RMSEP)
cat("\nMn_ICP RMSEP:\n")
print(results_CV$Mn_ICP$RMSEP)
cat("\nFe_ICP RMSEP:\n")
print(results_CV$Fe_ICP$RMSEP)
cat("\nSr_ICP RMSEP:\n")
print(results_CV$Sr_ICP$RMSEP)
cat("\nZr_ICP RMSEP:\n")
print(results_CV$Zr_ICP$RMSEP)
cat("\n\nJackknife test results_CV for ICP elements:\n\n")
cat("Ti_ICP jackknife_test:\n")
print(results_CV$Ti_ICP$jackknife_test)
cat("\nCa_ICP jackknife_test:\n")
print(results_CV$Ca_ICP$jackknife_test)
cat("\nMn_ICP jackknife_test:\n")
print(results_CV$Mn_ICP$jackknife_test)
cat("\nFe_ICP jackknife_test:\n")
print(results_CV$Fe_ICP$jackknife_test)
cat("\nSr_ICP jackknife_test:\n")
print(results_CV$Sr_ICP$jackknife_test)
cat("\nZr_ICP jackknife_test:\n")
print(results_CV$Zr_ICP$jackknife_test)
sink()






# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

# 1.2 PSLR training model with Ti, Ca, Sr, Zr as element matrix ----------------

# -------------------------------------------------------------------------

# Select dataset and convert to base R dataframe for pls package ---------------
df0 <- ACE_matched_xrf_icp_cps %>%
  select(c(Site, depth, SH20_mean_age, SH20_mean_95CI, all_of(xrf_icp_elements1))) %>%
  drop_na()
df0
df <- as.data.frame(df0) # Convert tibble to data frame for base R
# use original cps dataset and create a need clr matrix based on [signifcant elements] only!
# from log_inc dataset, only "Ti", "Ca", "Sr", "Zr" had p<0.001 in CV & LOO validation tests

# Define predictor ITRAX variables & ICPMS response variable for calibration model tests
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
main_elements_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP","Fe_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP","Zr_ICP")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")

# Select XRF and ICPMS input variables

# Example: assume df_xrf and df_icp are your original data frames

df_key_elements_xrf <- df[, c(key_elements_xrf)]
df_key_elements_icp <- df[, c(key_elements_icp)]
df_xrf <- df[, c("Ti", "Ca", "Sr", "Zr")]
df_icp <- df[, c("Ti_ICP", "Ca_ICP", "Sr_ICP", "Zr_ICP")]

# clr transformation of predictors and response ------------------------------

# Transform to compositions and add _clr to column titles
xrf_clr <- clr(acomp(df_xrf))
colnames(xrf_clr) <- paste0(colnames(xrf_clr), "_clr")
icp_clr <- clr(acomp(df_icp))
colnames(icp_clr) <- paste0(colnames(icp_clr), "_clr")

# Refit CV-PLSR (cross-validated) training model -------------------------------
plsr_model <- plsr(icp_clr ~ xrf_clr, ncomp = 4, validation = "CV", 
                   jackknife = TRUE, scale = TRUE, center = TRUE)

# Summary of final plsr model and jackknife t tests of regression coefficients & std errors
summary(plsr_model)
jack.test(plsr_model, ncomp = 4, use.mean = TRUE) # estimates and std errors not valid
rmsep_obj <- RMSEP(plsr_model) # check that RMSEP = RMSE for training dataset

# Extract RMSEP values as a matrix
rmsep_matrix <- rmsep_obj$val[, , ]  # usually a 3D array
# Convert to data frame (assuming 4th component response variable)
rmsep_df <- as.data.frame(t(rmsep_matrix[,,5]))  # in column index 5 due to intercept in col 1
rmsep_df <- tibble::rownames_to_column(rmsep_df, "Type")  # Add rownames as a column
rmsep_df
rmsep_labels_Ti <- rmsep_df %>% mutate(label = glue("CV = {round(CV, 3)}\nadjCV = {round(adjCV, 3)}"))
# How to interpret e.g. 4 comps Ti: CV = 0.3391, adjCV = 0.3589 - values are low and similar
# RMSE value for Ti is 0.3534 is most
# The RMSEP for 3 components ~0.331, meaning the average prediction error in the cross-validation set is 0.331 in clr space
# The adjusted estimate (0.311) confirms its reliability.

# Generate R2, RMSE, equation & labels for each element in base R ----------

# Predict Ti using the model
predicted_Ti0 <- predict(plsr_model, ncomp = 4)
Ti_predicted <- predicted_Ti0[, "Ti_ICP_clr", 1]
Ti_actual <- icp_clr$Ti_ICP_clr
length(Ti_predicted) == length(Ti_actual)  # should be TRUE
# R-squared
r2_Ti <- round(cor(Ti_actual, Ti_predicted)^2, 2)
# RMSE
rmse_Ti <- round(sqrt(mean((Ti_actual - Ti_predicted)^2)), 4)
# RMSEP from CV
rmsep_labels_Ti <- rmsep_df[rmsep_df$Type == "Ti_ICP_clr", ]
rmsep_label_text_Ti <- paste0("\n",
  "CV = ", round(rmsep_labels_Ti$CV, 4), "\n",
  "adjCV = ", round(rmsep_labels_Ti$adjCV, 4))
cat("Ti_R² =", r2_Ti, "\nTi_RMSE =", rmse_Ti, "\nTi_RMSEP:", rmsep_label_text_Ti)

# Predict Ca using the model
predicted_Ca0 <- predict(plsr_model, ncomp = 4)
Ca_predicted <- predicted_Ca0[, "Ca_ICP_clr", 1]
Ca_actual <- icp_clr$Ca_ICP_clr
length(Ca_predicted) == length(Ca_actual)  # should be TRUE
# R-squared
r2_Ca <- round(cor(Ca_actual, Ca_predicted)^2, 2)
# RMSE
rmse_Ca <- round(sqrt(mean((Ca_actual - Ca_predicted)^2)), 4)
# RMSEP from CV
rmsep_labels_Ca <- rmsep_df[rmsep_df$Type == "Ca_ICP_clr", ]
rmsep_label_text_Ca <- paste0("\n",
                              "CV = ", round(rmsep_labels_Ca$CV, 4), "\n",
                              "adjCV = ", round(rmsep_labels_Ca$adjCV, 4))
cat("Ca_R² =", r2_Ca, "\nCa_RMSE =", rmse_Ca, "\nCa_RMSEP:", rmsep_label_text_Ca)

# Predict Sr using the model
predicted_Sr0 <- predict(plsr_model, ncomp = 4)
Sr_predicted <- predicted_Sr0[, "Sr_ICP_clr", 1]
Sr_actual <- icp_clr$Sr_ICP_clr
length(Sr_predicted) == length(Sr_actual)  # should be TRUE
# R-squared
r2_Sr <- round(cor(Sr_actual, Sr_predicted)^2, 2)
# RMSE
rmse_Sr <- round(sqrt(mean((Sr_actual - Sr_predicted)^2)), 4)
# RMSEP from CV
rmsep_labels_Sr <- rmsep_df[rmsep_df$Type == "Sr_ICP_clr", ]
rmsep_label_text_Sr <- paste0("\n",
                              "CV = ", round(rmsep_labels_Sr$CV, 4), "\n",
                              "adjCV = ", round(rmsep_labels_Sr$adjCV, 4))
cat("Sr_R² =", r2_Sr, "\nSr_RMSE =", rmse_Sr, "\nSr_RMSEP:", rmsep_label_text_Sr)

# Predict Zr using the model
predicted_Zr0 <- predict(plsr_model, ncomp = 4)
Zr_predicted <- predicted_Zr0[, "Zr_ICP_clr", 1]
Zr_actual <- icp_clr$Zr_ICP_clr
length(Zr_predicted) == length(Zr_actual)  # should be TRUE
# R-squared
r2_Zr <- round(cor(Zr_actual, Zr_predicted)^2, 2)
# RMSE
rmse_Zr <- round(sqrt(mean((Zr_actual - Zr_predicted)^2)), 4)
# RMSEP from CV
rmsep_labels_Zr <- rmsep_df[rmsep_df$Type == "Zr_ICP_clr", ]
rmsep_label_text_Zr <- paste0("\n",
                              "CV = ", round(rmsep_labels_Zr$CV, 4), "\n",
                              "adjCV = ", round(rmsep_labels_Zr$adjCV, 4))
cat("Zr_R² =", r2_Zr, "\nZr_RMSE =", rmse_Zr, "\nZr_RMSEP:", rmsep_label_text_Zr)

# Get coefficients and intercept from the model
coef_matrix <- coef(plsr_model, ncomp = 4)
intercept <- attr(coef_matrix, "constant")

# Extract coefficients - Ti_ICP_clr column only as a named vector and retain keeping row names for equation
coef_Ti <- coef_matrix[, "Ti_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)
coef_Ca <- coef_matrix[, "Ca_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)
coef_Sr <- coef_matrix[, "Sr_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)
coef_Zr <- coef_matrix[, "Zr_ICP_clr", 1, drop = FALSE]
names(coef_Ti) <- rownames(coef_matrix)

# Coefs for clr are meaningless - can't convert clr coefficients into ppm directly.
# Need to apply the model to get predictions in clr
# Then back-transform predictions to composition (with clrInv()),
# Then multiply by total ICP to get back to ppm - see below and section 1.3 below for how this is done

# Show regression equation - this is meaningless for clr - conversion to ppm needs inversion of matrix first
cat("PLS clr Regression Equation:\n")
cat("Ti[ICP-MS_pred] clr = ",  paste0(round(coef_Ti, 4), "*", names(coef_Ti), collapse = " + "), "\n") #intercept removed as NULL 


# Get CV predictions as clr, plot --------------------------------------------
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

# Make RMSE and R² and CI summary table for each element in Tidyverse ---------------

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
library(glue) # to make next labels for plotting
ci_df <- bind_rows(ci_list)
summary_df_clr_CI <- bind_cols(summary_df_clr, ci_df, rmsep_df) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>% 
  mutate(Element = str_replace(Element, "_ICP_clr_actual", "")) %>%
  mutate(label_R2 = paste0("R^2:", round(R2, 2))) %>%  
  mutate(label_RMSE = paste0("RMSE: ", round(RMSE, 4))) %>% 
  mutate(label_RMSEP_CV = paste0("RMSEP_CV[4]:", round(CV, 4))) %>% 
  mutate(label_RMSEP_adjCV = paste0("RMSEP_adjCV[4]:", round(adjCV, 4)))
print(summary_df_clr_CI)
write_csv(summary_df_clr_CI, "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Summary_clr_R2_RMSE_CI.csv")

# Create plot df of PLSR measured vs predicted as clr  ------------------
plot_df_clr <- bind_cols(df, actual_clr, pred_clr, pred_clr_lower, pred_clr_upper)

# Plot Ti training model - meas vs pred as clr values --------------------------

# Make labels
R2_text_Ti <- summary_df_clr_CI %>%
  filter(Element == "Ti") %>% 
  pull(label_R2)
RMSE_text_Ti <- summary_df_clr_CI %>%
  filter(Element == "Ti") %>% 
  pull(label_RMSE)
RMSEP_CV_text_Ti <- summary_df_clr_CI %>% 
  filter(Element == "Ti") %>% 
  pull(label_RMSEP_CV)
RMSEP_adjCV_text_Ti <- summary_df_clr_CI %>% 
  filter(Element == "Ti") %>% 
  pull(label_RMSEP_adjCV)
# Plot
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ti_clr <- ggplot(plot_df_clr, aes(x = Ti_ICP_clr_pred, y = Ti_ICP_clr_actual, color = Site)) +
  scale_color_jco() +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(title = "PLSR: Measured vs Predicted Ti: clr model",
       y = "Measured Ti [ICP-MS] (clr)",
       x = "Predicted Ti [clr model]") +
  coord_equal(xlim = c(-0.5, 2.5), ylim = c(-0.5, 2.5)) +
  #coord_equal() +
  #theme_minimal()
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = -0.5, y = 2.4, label = R2_text_Ti, parse = TRUE, hjust = 0) +
  annotate("text", x = -0.5, y = 2.3, label = RMSE_text_Ti, parse = TRUE, hjust = 0) + 
  annotate("text", x = -0.5, y = 2.2, label = RMSEP_CV_text_Ti, parse = TRUE, hjust = 0) +
  annotate("text", x = -0.5, y = 2.1, label = RMSEP_adjCV_text_Ti, parse = TRUE, hjust = 0)
print(PLS_final_Ti_clr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Ti_plsr_clr_as_clr.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Ca training model - pred vs meas as clr values --------------------------

# Make labels
R2_text_Ca <- summary_df_clr_CI %>%
  filter(Element == "Ca") %>% 
  pull(label_R2)
RMSE_text_Ca <- summary_df_clr_CI %>%
  filter(Element == "Ca") %>% 
  pull(label_RMSE)
RMSEP_CV_text_Ca <- summary_df_clr_CI %>% 
  filter(Element == "Ca") %>% 
  pull(label_RMSEP_CV)
RMSEP_adjCV_text_Ca <- summary_df_clr_CI %>% 
  filter(Element == "Ca") %>% 
  pull(label_RMSEP_adjCV)
# Plot
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ca_clr <- ggplot(plot_df_clr, aes(x = Ca_ICP_clr_pred, y = Ca_ICP_clr_actual, color = Site)) +
  scale_color_jco() +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(title = "PLSR: Measured vs Predicted Ca: clr model",
       y = "Measured Ca [ICP-MS] (clr)",
       x = "Predicted Ca [clr model]") +
  coord_equal(xlim = c(1.5, 4.5), ylim = c(1.5, 4.5)) +
  #coord_equal() +
  #theme_minimal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 1.5, y = 4.4, label = R2_text_Ca, parse = TRUE, hjust = 0) +
  annotate("text", x = 1.5, y = 4.3, label = RMSE_text_Ca, parse = TRUE, hjust = 0) + 
  annotate("text", x = 1.5, y = 4.2, label = RMSEP_CV_text_Ca, parse = TRUE, hjust = 0) +
  annotate("text", x = 1.5, y = 4.1, label = RMSEP_adjCV_text_Ca, parse = TRUE, hjust = 0)
print(PLS_final_Ca_clr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Ca_plsr_clr_as_clr.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# Plot Sr training model - pred vs meas as clr values --------------------------

# Make labels
R2_text_Sr <- summary_df_clr_CI %>%
  filter(Element == "Sr") %>% 
  pull(label_R2)
RMSE_text_Sr <- summary_df_clr_CI %>%
  filter(Element == "Sr") %>% 
  pull(label_RMSE)
RMSEP_CV_text_Sr <- summary_df_clr_CI %>% 
  filter(Element == "Sr") %>% 
  pull(label_RMSEP_CV)
RMSEP_adjCV_text_Sr <- summary_df_clr_CI %>% 
  filter(Element == "Sr") %>% 
  pull(label_RMSEP_adjCV)
# Plot
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Sr_clr <- ggplot(plot_df_clr, aes(x = Sr_ICP_clr_pred, y = Sr_ICP_clr_actual, color = Site)) +
  scale_color_jco() +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(title = "PLSR: Measured vs Predicted Sr: clr model",
       y = "Measured Sr [ICP-MS] (clr)",
       x = "Predicted Sr [clr model]") +
  coord_equal(xlim = c(-3.5, 0.5), ylim = c(-3.5, 0.5)) +
  #coord_equal() + # coord_equal()
  #theme_minimal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = -3.2, y = 0.4, label = R2_text_Sr, parse = TRUE, hjust = 0) +
  annotate("text", x = -3.2, y = 0.2, label = RMSE_text_Sr, parse = TRUE, hjust = 0) + 
  annotate("text", x = -3.2, y = 0, label = RMSEP_CV_text_Sr, parse = TRUE, hjust = 0) +
  annotate("text", x = -3.2, y = -0.2, label = RMSEP_adjCV_text_Sr, parse = TRUE, hjust = 0)
print(PLS_final_Sr_clr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Sr_plsr_clr_as_clr.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Zr training model - pred vs meas as clr values --------------------------

# Make labels
R2_text_Zr <- summary_df_clr_CI %>%
  filter(Element == "Zr") %>% 
  pull(label_R2)
RMSE_text_Zr <- summary_df_clr_CI %>%
  filter(Element == "Zr") %>% 
  pull(label_RMSE)
RMSEP_CV_text_Zr <- summary_df_clr_CI %>% 
  filter(Element == "Zr") %>% 
  pull(label_RMSEP_CV)
RMSEP_adjCV_text_Zr <- summary_df_clr_CI %>% 
  filter(Element == "Zr") %>% 
  pull(label_RMSEP_adjCV)
# Plot
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Zr_clr <- ggplot(plot_df_clr, aes(x = Zr_ICP_clr_pred, y = Zr_ICP_clr_actual , color = Site)) +
  scale_color_jco() +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(title = "PLSR: Measured vs Predicted Zr: clr model",
       y = "Measured Zr [ICP-MS] (clr)",
       x = "Predicted Zr [clr model]") +
  coord_equal(xlim = c(-4.5, 1), ylim = c(-4.5, 1)) +
  #coord_equal() +
  #theme_minimal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = -4.5, y = -0, label = R2_text_Zr, parse = TRUE, hjust = 0) +
  annotate("text", x = -4.5, y = -0.2, label = RMSE_text_Zr, parse = TRUE, hjust = 0) + 
  annotate("text", x = -4.5, y = -0.4, label = RMSEP_CV_text_Zr, parse = TRUE, hjust = 0) +
  annotate("text", x = -4.5, y = -0.6, label = RMSEP_adjCV_text_Zr, parse = TRUE, hjust = 0)
print(PLS_final_Zr_clr)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Zr_plsr_clr_as_clr.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")




# -------------------------------------------------------------------------
# Convert actual and predicted Ti clr values back to ppm values  --------

# Invert clr transform to get new dataset values in ppm
pred_acomp <- clrInv(pred_clr)
pred_upper_acomp <- clrInv(pred_clr_upper)
pred_lower_acomp <- clrInv(pred_clr_lower)
actual_acomp <- clrInv(actual_clr)
total_icp <- rowSums(df_icp) #Get total ICP per sample to invert clr and convert to ppm

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
  mutate(total_icp = total_icp) %>%
  rowwise() %>%
  mutate(across(-total_icp, ~ .x * total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm")) %>%
  mutate(across(everything(), ~ as.numeric(unlist(.x)))) %>% 
  rename(Ca_ICP_clr_pred_lower_ppm = Ca_ICP_clr_pred_upper_ppm)
pred_upper_ppm

pred_lower_ppm <- pred_lower_acomp %>%
  as_tibble() %>%
  mutate(total_icp = total_icp) %>%
  rowwise() %>%
  mutate(across(-total_icp, ~ .x * total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm")) %>%
  mutate(across(everything(), ~ as.numeric(unlist(.x)))) %>% 
  rename(Ca_ICP_clr_pred_upper_ppm = Ca_ICP_clr_pred_lower_ppm)
pred_lower_ppm

# Make RMSE and R² and CI summary table for each element in Tidyverse ---------------

summary_df_ppm <- map2_dfr(actual_ppm, pred_ppm, ~ {
  tibble(
    R2 = cor(.x, .y)^2,
    RMSE = sqrt(mean((.x - .y)^2))
  )
}, .id = "Element")

# Rename based on actual column names
summary_df_ppm$Element <- names(actual_ppm)
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
library(glue) # to make next labels for plotting
ci_df <- bind_rows(ci_list)
summary_df_ppm_CI <- bind_cols(summary_df_ppm, ci_df) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  mutate(Element = str_replace(Element, "_ICP_actual_ppm", "")) %>%
  mutate(label_R2 = paste0("R^2:", round(R2, 2))) %>%  
  mutate(label_RMSE = paste0("RMSE: ", round(RMSE, 2)))
print(summary_df_ppm_CI)
write_csv(summary_df_ppm_CI, "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Summary_ppm_R2_RMSE_CI.csv")

# Set up Plot df for PLSR measured vs predicted as ppm 
plot_df_ppm <- bind_cols(df, actual_ppm, pred_ppm, pred_upper_ppm, pred_lower_ppm)
plot_df_ppm

# Plot Ti training model - pred vs measured plot as ppm ---------------------------------------

# Make labels
R2_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Ti") %>% 
  pull(label_R2)
RMSE_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Ti") %>% 
  pull(label_RMSE)
# Plot
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
  labs(
    title = expression("PLS Ti clr model: invclr*mtc to"~"mg kg"^{-1}*""),
    y = expression("Measured Ti [ICP-MS]" ~"(mg kg"^{-1}*")"),
    x = expression("Predicted Ti"~"(mg kg"^{-1}*") [invclr*mtc]")
  ) +
  coord_equal(xlim = c(0, 40000), ylim = c(0, 40000)) + # coord_equal()
  #theme_minimal()
  #coord_equal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0, y = 40000, label = R2_text_Ti_ppm, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 38000, label = RMSE_text_Ti_ppm, parse = TRUE, hjust = 0)
  print(PLS_final_Ti_ppm)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Ti_plsr_clr_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Ca training model - pred vs measured plot as ppm ---------------------------------------

# Make labels
R2_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Ca") %>% 
  pull(label_R2)
RMSE_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Ca") %>% 
  pull(label_RMSE)
# Plot
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ti_ppm <- ggplot(plot_df_ppm, aes(x = Ca_ICP_actual_ppm, y = Ca_ICP_clr_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(
    title = expression("PLS Ca clr model: invclr*mtc to"~"mg kg"^{-1}*""),
    y = expression("Measured Ca [ICP-MS]" ~"(mg kg"^{-1}*")"),
    x = expression("Predicted Ca"~"(mg kg"^{-1}*") [invclr*mtc]")
  ) +
  coord_equal(xlim = c(0, 80000), ylim = c(0, 80000)) + # coord_equal()
  #theme_minimal()
  #coord_equal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0, y = 80000, label = R2_text_Ti_ppm, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 78000, label = RMSE_text_Ti_ppm, parse = TRUE, hjust = 0)
print(PLS_final_Ti_ppm)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Ca_plsr_clr_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Sr training model - pred vs measured plot as ppm ---------------------------------------

# Make labels
R2_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Sr") %>% 
  pull(label_R2)
RMSE_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Sr") %>% 
  pull(label_RMSE)
# Plot
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ti_ppm <- ggplot(plot_df_ppm, aes(x = Sr_ICP_actual_ppm, y = Sr_ICP_clr_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(
    title = expression("PLS Sr clr model: invclr*mtc to"~"mg kg"^{-1}*""),
    y = expression("Measured Sr [ICP-MS]" ~"(mg kg"^{-1}*")"),
    x = expression("Predicted Sr"~"(mg kg"^{-1}*") [invclr*mtc]")
  ) +
  coord_equal(xlim = c(0, 700), ylim = c(0, 700)) + # coord_equal()
  #theme_minimal()
  #coord_equal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0, y = 700, label = R2_text_Ti_ppm, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 680, label = RMSE_text_Ti_ppm, parse = TRUE, hjust = 0)
print(PLS_final_Ti_ppm)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Sr_plsr_clr_as_ppm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Plot Zr training model - pred vs measured plot as ppm ---------------------------------------

# Make labels
R2_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Zr") %>% 
  pull(label_R2)
RMSE_text_Ti_ppm <- summary_df_ppm_CI %>%
  filter(Element == "Zr") %>% 
  pull(label_RMSE)
# Plot
library(ggpmisc)  # for stat_poly_eq()
theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final_Ti_ppm <- ggplot(plot_df_ppm, aes(x = Zr_ICP_actual_ppm, y = Zr_ICP_clr_pred_ppm, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
  #             formula = y ~ x, parse = TRUE,
  #             label.x = "left", label.y = "top") +
  labs(
    title = expression("PLS Zr clr model: invclr*mtc to"~"mg kg"^{-1}*""),
    y = expression("Measured Zr [ICP-MS]" ~"(mg kg"^{-1}*")"),
    x = expression("Predicted Zr"~"(mg kg"^{-1}*") [invclr*mtc]")
  ) +
  coord_equal(xlim = c(0, 500), ylim = c(0, 500)) + # coord_equal()
  #theme_minimal()
  #coord_equal()
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0, y = 500, label = R2_text_Ti_ppm, parse = TRUE, hjust = 0) + 
  annotate("text", x = 0, y = 480, label = RMSE_text_Ti_ppm, parse = TRUE, hjust = 0)
print(PLS_final_Ti_ppm)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/Zr_plsr_clr_as_ppm.pdf", 
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

# Calculate mean total of clr from the training ICPMS dataset (n=265) with known total row sums
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
RMSE_clr <- c(Ti_ICP_clr = 2*rmse_Ti, Ca_ICP_clr = 2*rmse_Ca, Sr_ICP_clr = 2*rmse_Sr, Zr_ICP_clr = 2*rmse_Zr)

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
  rename_with(~ paste0(.x, "_ppm")) %>%  
  mutate(across(everything(), ~ as.numeric(unlist(.x)))) %>% 
  rename(Ca_ICP_clr_lower_ppm = Ca_ICP_clr_upper_ppm)
pred_xrf_upper_ppm

pred_xrf_lower_ppm <- pred_xrf_lower_acomp %>%
  as_tibble() %>%
  mutate(mean_total_icp = mean_total_icp) %>%
  rowwise() %>%
  mutate(across(-mean_total_icp, ~ .x * mean_total_icp)) %>% #multiply each column by mean total_icp
  ungroup() %>%
  select(-mean_total_icp) %>% 
  rename_with(~ paste0(.x, "_ppm")) %>%  
  mutate(across(everything(), ~ as.numeric(unlist(.x)))) %>% 
  rename(Ca_ICP_clr_upper_ppm = Ca_ICP_clr_lower_ppm)
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
write.csv(ACE_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv", row.names = FALSE)

# ------------------------------------------------------------------------------

# 1.4 Site datasets & Depth, Age plots with CI errors --------------------------
 # TO FINISH

# Plotting functions - dDepth and age with ICP-MS overlay  ------------------
plot_element_depth_clr_with_ICP <- function(site_name, element, pred_data, 
                                                icp_data, output_dir, age_model, 
                                                age_breaks = NULL) {
  
  # make sure tidypaleo is loaded for age secondary axis
  library(tidypaleo)
  
  # Define age breaks automatically if not provided
  if (is.null(age_breaks)) {
    age_range <- range(age_model$age, na.rm = TRUE)
    age_breaks <- pretty(age_range, n = 6)
  }
  
  # Dynamically construct column names
  element_pred <- sym(paste0(element, "_ICP_clr_ppm"))
  element_low  <- sym(paste0(element, "_low_RMSE"))
  element_up   <- sym(paste0(element, "_up_RMSE"))
  element_icp  <- sym(paste0(element, "_ICP"))
  
  # Build plot
  p <- ggplot(pred_data, aes(y = depth)) +
    geom_ribbon(aes(xmin = !!element_low, xmax = !!element_up), color = "grey80", fill = "grey80", na.rm = TRUE) + 
    #geom_point(aes(x = !!element_pred), color = "blue", size = 0.3) +
    geom_lineh(aes(x = !!element_pred), linewidth = 0.5, color = "blue", na.rm = TRUE) + 
    geom_lineh(data = icp_data, aes(x = !!element_icp, y = depth), color = "red", linewidth = 0.75, na.rm = TRUE) +
    geom_errorbarh(data = icp_data,aes_string(y = "depth", 
                                              xmin = paste0(element, "_ICP - ", element, "_ICP_sd"),
                                              xmax = paste0(element, "_ICP + ", element, "_ICP_sd")),
                   height = 0, color = "red", linewidth = 0.75) +
    geom_point(data = icp_data, aes(x = !!element_icp, y = depth), shape = 21, size = 0.75, color = "red", 
               fill = "white", stroke = 0.75, na.rm = TRUE) +  
    #scale_y_reverse() +
    scale_y_depth_age(
      age_model,
      name = "Depth (cm)",
      age_name = "Age (cal a BP)",
      age_breaks = age_breaks) +
    labs(
      title = paste0(site_name, ": PLSR ", element, " Predictions with Confidence Interval (log_inc model)"),
      x = bquote(.(element) ~ "(mg kg"^{-1}*")") #,
      #y = "Depth (cm)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(color = "black", size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.ticks.length = unit(0.3, "cm"),
      axis.line = element_line(size = 0.75),
      axis.ticks = element_line(size = 0.75)
    )
  
  # Save plot
  if (!is.null(output_dir)) {
    file_path <- file.path(output_dir, paste0(site_name, "_", element, "_log_inc_PLS_depth.pdf"))
    ggsave(file_path, plot = p, height = 24, width = 12, dpi = 600, units = "cm")
  }
  return(p)
}
plot_element_age_clr_with_ICP <- function(site_name, element, pred_data, icp_data, output_dir) {
  # Dynamically construct column names
  element_pred <- sym(paste0(element, "_ICP_clr_ppm"))
  element_low  <- sym(paste0(element, "_low_RMSE"))
  element_up   <- sym(paste0(element, "_up_RMSE"))
  element_icp  <- sym(paste0(element, "_ICP"))
  # Build plot
  p <- ggplot(pred_data, aes(x = SH20_age)) +
    geom_ribbon(aes(ymin = !!element_low, ymax = !!element_up), color = "grey80", fill = "grey80", na.rm = TRUE) + 
    #geom_point(aes(y = !!element_pred), color = "blue", size = 0.3) +
    geom_line(aes(y = !!element_pred), linewidth = 0.5, color = "blue", na.rm = TRUE) +
    geom_line(data = icp_data, aes(x = SH20_age, y = !!element_icp), color = "red", linewidth = 1, na.rm = TRUE) +
    geom_errorbar(data = icp_data,aes_string(x = "SH20_age",
                                             ymin = paste0(element, "_ICP", " - ", element, "_ICP_sd"),
                                             ymax = paste0(element, "_ICP", " + ", element, "_ICP_sd")),
                  width = 0, color = "red", linewidth = 1) +
    geom_point(data = icp_data, aes(x = SH20_age, y = !!element_icp), shape = 21, size = 1, color = "red", 
               fill = "white", na.rm = TRUE) + #stroke = 0.5, 
    scale_x_reverse() +
    labs(
      title = paste0(site_name, ": PLSR ", element, " Predictions with Confidence Interval (log_inc model)"),
      x = "Age (cal yr BP)",
      y = bquote("Predicted" ~ .(element) ~ "(mg kg"^{-1}*") [log_inc]")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(color = "black", size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.ticks.length = unit(0.3, "cm"),
      axis.line = element_line(size = 0.75),
      axis.ticks = element_line(size = 0.75)
    )
  # Save plot
  if (!is.null(output_dir)) {
    file_path <- file.path(output_dir, paste0(site_name, "_", element, "_log_inc_PLS_age.pdf"))
    ggsave(file_path, plot = p, height = 12, width = 24, dpi = 600, units = "cm")
  }
  return(p)
}

# -------------------------------------------------------------------------
# Ti ---------------------------------------------------------------------------

# -------------------------------------------------------------------------
# HER42PB Import ---------------------------------------------------------------

# Import HER42PB clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "HER42PB")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ti data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "HER42PB")
HER42PB_ICP

# HER42PB Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# HER42PB Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "HER42PB",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("HER42PB", "Ti", HER42PB_ppm_PLS_pred_clr, 
                                                    HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# BI10 Import ---------------------------------------------------------------

# Import BI10 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "BI10")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ti data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "BI10")
HER42PB_ICP

# BI10 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# BI10 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "BI10",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("BI10", "Ti", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER1 Import ---------------------------------------------------------------

# Import KER1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ti data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER1")
HER42PB_ICP

# KER1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER1",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER1", "Ti", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER3 Import ---------------------------------------------------------------

# Import KER3 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER3")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ti data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER3")
HER42PB_ICP

# KER3 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER3 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER3",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER3", "Ti", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------
# PB1 Import ---------------------------------------------------------------

# Import PB1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "PB1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ti data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "PB1")
HER42PB_ICP

# PB1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# PB1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "PB1",
  element = "Ti",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("PB1", "Ti", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Ca ---------------------------------------------------------------------------

# -------------------------------------------------------------------------
# HER42PB Import ---------------------------------------------------------------

# Import HER42PB clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "HER42PB")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ca data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "HER42PB")
HER42PB_ICP

# HER42PB Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# HER42PB Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "HER42PB",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("HER42PB", "Ca", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# BI10 Import ---------------------------------------------------------------

# Import BI10 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "BI10")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ca data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "BI10")
HER42PB_ICP

# BI10 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# BI10 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "BI10",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("BI10", "Ca", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER1 Import ---------------------------------------------------------------

# Import KER1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ca data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER1")
HER42PB_ICP

# KER1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER1",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER1", "Ca", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER3 Import ---------------------------------------------------------------

# Import KER3 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER3")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ca data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER3")
HER42PB_ICP

# KER3 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER3 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER3",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER3", "Ca", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------
# PB1 Import ---------------------------------------------------------------

# Import PB1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "PB1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Ca data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "PB1")
HER42PB_ICP

# PB1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# PB1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "PB1",
  element = "Ca",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("PB1", "Ca", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Sr ---------------------------------------------------------------------------

# -------------------------------------------------------------------------
# HER42PB Import ---------------------------------------------------------------

# Import HER42PB clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "HER42PB")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Sr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "HER42PB")
HER42PB_ICP

# HER42PB Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# HER42PB Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "HER42PB",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("HER42PB", "Sr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# BI10 Import ---------------------------------------------------------------

# Import BI10 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "BI10")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Sr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "BI10")
HER42PB_ICP

# BI10 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# BI10 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "BI10",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("BI10", "Sr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER1 Import ---------------------------------------------------------------

# Import KER1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Sr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER1")
HER42PB_ICP

# KER1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER1",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER1", "Sr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER3 Import ---------------------------------------------------------------

# Import KER3 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER3")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Sr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER3")
HER42PB_ICP

# KER3 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER3 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER3",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER3", "Sr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------
# PB1 Import ---------------------------------------------------------------

# Import PB1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "PB1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Sr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "PB1")
HER42PB_ICP

# PB1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# PB1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "PB1",
  element = "Sr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("PB1", "Sr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Zr ---------------------------------------------------------------------------

# -------------------------------------------------------------------------
# HER42PB Import ---------------------------------------------------------------

# Import HER42PB clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "HER42PB")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Zr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "HER42PB")
HER42PB_ICP

# HER42PB Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "HER42PB",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# HER42PB Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "HER42PB",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("HER42PB", "Zr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# BI10 Import ---------------------------------------------------------------

# Import BI10 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "BI10")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Zr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "BI10")
HER42PB_ICP

# BI10 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "BI10",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# BI10 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "BI10",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("BI10", "Zr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER1 Import ---------------------------------------------------------------

# Import KER1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Zr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER1")
HER42PB_ICP

# KER1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER1",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER1",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER1", "Zr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)

# -------------------------------------------------------------------------
# KER3 Import ---------------------------------------------------------------

# Import KER3 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "KER3")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Zr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "KER3")
HER42PB_ICP

# KER3 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "KER3",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# KER3 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "KER3",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("KER3", "Zr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------
# PB1 Import ---------------------------------------------------------------

# Import PB1 clr predicted ppm data 
HER42PB_ppm_PLS_pred_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/ACE_PLS_ppm_pred_clr.csv") %>% 
  filter(Site == "PB1")
HER42PB_ppm_PLS_pred_clr

# Import ICPMS Zr data for comparison / overlay
HER42PB_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, depth, SH20_mean_age, all_of(icp_Elements_min), all_of(icp_Elements_min_sd)) %>% 
  rename(SH20_age = SH20_mean_age) %>% 
  filter(Site == "PB1")
HER42PB_ICP

# PB1 Depth - ICPMS overlay ------------------------------------------------
# Inputs for depth function - save to file

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_ppm_PLS_pred_clr, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_clr_with_ICP(
  site_name = "PB1",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)

# PB1 Age - ICPMS overlay ---------------------------------------------------
# Inputs for depth function - save to file
plot_element_age_clr_with_ICP(
  site_name = "PB1",
  element = "Zr",
  pred_data = HER42PB_ppm_PLS_pred_clr,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1"
)
# Create site_depth for plot matrix
HER42PB_age_Ti <- plot_element_age_clr_with_ICP("PB1", "Zr", HER42PB_ppm_PLS_pred_clr, 
                                                HER42PB_ICP, output_dir = NULL)


# -------------------------------------------------------------------------



# END --------------------------------------------------------------------------

# old code
# Extract HER42PB data ----------------------------------------------------

HER42PB_ppm_PLS_pred_clr <- ACE_ppm_PLS_pred_clr %>% 
  filter(Site == "HER42PB")
HER42PB_ppm_PLS_pred_clr
write.csv(HER42PB_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB_PLS_ppm_pred_clr.csv", row.names = FALSE)

# HER42PB - Single depth plot with RMSE errors - black and grey ----------------
ggplot(HER42PB_ppm_PLS_pred_clr, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_clr_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_clr_ppm)) +
  scale_y_reverse() +
  labs(title = "HER42PB: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Depth (cm)", y = "Predicted Ti (ppm) [clr model]") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/log_inc/HER42PB_Ti_clr_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(HER42PB_ppm_PLS_pred_clr, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_clr_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_clr_ppm)) +
  labs(title = "HER42PB: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Age (cal yr BP)", y = "Predicted Ti (ppm) [clr model]") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/HER42PB_Ti_clr_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")


# Extract BI10 data -----------------------------------------------------
BI10_ppm_PLS_pred_clr <- ACE_ppm_PLS_pred_clr %>% 
  filter(Site == "BI10")
BI10_ppm_PLS_pred_clr
write.csv(BI10_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10_PLS_ppm_pred_clr.csv", row.names = FALSE)

# BI10 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(BI10_ppm_PLS_pred_clr, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_clr_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_clr_ppm)) +
  scale_y_reverse() +
  labs(title = "BI10: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Depth (cm)", y = "Predicted Ti (ppm) [clr model]") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/log_inc/BI10_Ti_clr_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(BI10_ppm_PLS_pred_clr, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_clr_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_clr_ppm)) +
  labs(title = "BI10: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Age (cal yr BP)", y = "Predicted Ti (ppm) [clr model]") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/BI10_Ti_clr_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")



# Extract KER1 data -------------------------------------------------------

KER1_ppm_PLS_pred_clr <- ACE_ppm_PLS_pred_clr %>% 
  filter(Site == "KER1")
KER1_ppm_PLS_pred_clr
write.csv(KER1_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1_PLS_ppm_pred_clr.csv", row.names = FALSE)

# KER1 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(KER1_ppm_PLS_pred_clr, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_clr_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_clr_ppm)) +
  scale_y_reverse() +
  labs(title = "KER1: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Depth (cm)", y = "Predicted Ti (ppm) [clr model]") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/log_inc/KER1_Ti_clr_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(KER1_ppm_PLS_pred_clr, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_clr_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_clr_ppm)) +
  labs(title = "KER1: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Age (cal yr BP)", y = "Predicted Ti (ppm) [clr model]") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER1_Ti_clr_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")



# Extract KER3 data -------------------------------------------------------

KER3_ppm_PLS_pred_clr <- ACE_ppm_PLS_pred_clr %>% 
  filter(Site == "KER3")
KER3_ppm_PLS_pred_clr
write.csv(KER3_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3_PLS_ppm_pred_clr.csv", row.names = FALSE)

# KER3 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(KER3_ppm_PLS_pred_clr, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_clr_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_clr_ppm)) +
  scale_y_reverse() +
  labs(title = "KER3: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Depth (cm)", y = "Predicted Ti (ppm) [clr model]") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/log_inc/KER3_Ti_clr_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(KER3_ppm_PLS_pred_clr, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_clr_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_clr_ppm)) +
  labs(title = "KER3: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Age (cal yr BP)", y = "Predicted Ti (ppm) [clr model]") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/KER3_Ti_clr_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")



# Extract PB1 data --------------------------------------------------------


PB1_ppm_PLS_pred_clr <- ACE_ppm_PLS_pred_clr %>% 
  filter(Site == "PB1")
PB1_ppm_PLS_pred_clr
write.csv(PB1_ppm_PLS_pred_clr,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1_PLS_ppm_pred_clr.csv", row.names = FALSE)

# PB1 - Single depth plot with RMSE errors - black and grey ----------------
ggplot(PB1_ppm_PLS_pred_clr, aes(y = depth)) +
  geom_ribbon(aes(xmin = Ti_low_RMSE, xmax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(x = Ti_ICP_clr_ppm), size = 0.3) +
  geom_lineh(aes(x = Ti_ICP_clr_ppm)) +
  scale_y_reverse() +
  labs(title = "PB1: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Depth (cm)", y = "Predicted Ti (ppm) [clr model]") +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/log_inc/PB1_Ti_clr_PLS_depth.pdf",
       height = c(24), width = c(12), dpi = 600, units = "cm")

ggplot(PB1_ppm_PLS_pred_clr, aes(x = SH20_age)) +
  geom_ribbon(aes(ymin = Ti_low_RMSE, ymax = Ti_up_RMSE), alpha = 0.5, color = "lightgrey") +
  geom_point(aes(y = Ti_ICP_clr_ppm), size = 0.3) +
  geom_line(aes(y = Ti_ICP_clr_ppm)) +
  labs(title = "PB1: PLSR Ti Predictions with Confidence Interval (clr model)",
       x = "Age (cal yr BP)", y = "Predicted Ti (ppm) [clr model]") +
  scale_x_reverse() +
  theme_minimal()
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/clr/PB1_Ti_clr_PLS_age.pdf",
       height = c(12), width = c(24), dpi = 600, units = "cm")







