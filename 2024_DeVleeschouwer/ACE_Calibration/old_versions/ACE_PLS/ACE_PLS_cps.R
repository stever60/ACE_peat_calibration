# ACE PLSR cps using pls package 

# XRF cps model and Ti_ICP response variable to scaled and centred cps dataset

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

# Assumptions:
# Ti_ICPMS is the response variable.
# Other elements (e.g., Ca, Fe, Sr, Mn, etc.) measured by ITRAX are predictor variables.
# All data are numeric and cleaned (no missing values).
# The response variable is continuous.

# Load and Prepare data
ACE_matched_xrf_icp_cps <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") 
is.na(ACE_matched_xrf_icp_cps)<-sapply(ACE_matched_xrf_icp_cps, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_cps

# Standardised and centred dataframe - Z-scores centred around 0
ACE_matched_xrf_icp_Z <- ACE_matched_xrf_icp_cps
ACE_matched_xrf_icp_Z[, xrf_icp_elements] <- scale(ACE_matched_xrf_icp_cps[, xrf_icp_elements], center = TRUE, scale = TRUE)
ACE_matched_xrf_icp_Z

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
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")

# Split into predictors (X) using key elements only and response (Y) as ICPMS element e.g., Ti_ICP
X_test <- df[, key_elements_xrf]  # ITRAX elements
Y_test <- df$Ti_ICP #define response ICP-MS variable for calibration

# Write summary outputs from console to txt file in Output folder instead of console
#sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/cps/Summary_stats.txt")

# LOO-PLSR (leave-one-out) test model to determine optimal number of comps -----

# Perform PLS regression with LOO - leave one out validation - useful for unbiased prediction error assessment
pls_train_model_LOO <- plsr(Y_test ~ ., data = data.frame(X_test, Y_test), validation = "LOO", jackknife = TRUE, scale = TRUE, center = TRUE) # Use scale = TRUE to perform as Z scores - do not use for calibration

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

# Final PSLR model --------------------------------------------------------

# Only use signifcant predictor ITRAX variables from LOO & CV jackknife results for final model
final_elements_xrf <- c("Ti", "Ca", "Sr")
opt_comp_final = 3

# Split into predictors (X) using key elements only and response (Y) as Ti_ICP
X_final <- df[, final_elements_xrf]  # ITRAX elements
Y_final <- df$Ti_ICP  # Define response ICP-MS variable for calibration

# Refit model with optimal number of components
pls_final <- plsr(Y_final ~ ., data = X_final, ncomp = opt_comp_final, scale = TRUE, center = TRUE, validation = "CV", jackknife = TRUE,)

# Summary of final plsr model and jackknife t tests of regression coefficients & std errors
summary(pls_final)
jack.test(pls_final, ncomp = opt_comp_final, use.mean = TRUE)

# Predict Ti using the model
predicted_Ti <- predict(pls_final, ncomp = opt_comp_final)

# Plot predicted vs actual draft output
plot(Y_final, predicted_Ti, xlab = "Measured Ti (ICPMS)", ylab = "Predicted Ti (from ITRAX)", main = "PLS Calibration cps")
abline(0, 1, col = "red")

# Calculate R2, RMSE, RMSEP
r2 <- cor(Y_final, predicted_Ti)^2
rmse <- sqrt(mean((Y_final - predicted_Ti)^2))
rmsep_result <- RMSEP(pls_final)
print(rmsep_result)
cat("R² =", r2, "\nRMSE =", rmse, "\n")

# Generate equation text for ggplot
R2_text <- paste0("R^2: ", paste0(round(r2, 2)))
RMSE_text <- paste0("RMSE: ", paste0(round(rmse, 4)))                  

# Coefficients of the final model
coefficients(pls_final)

# Get coefficients and intercept from the model
coef_matrix <- coef(pls_final, ncomp = opt_comp_final)
intercept <- attr(coef_matrix, "constant")

# Print the intercept - check = NULL
cat("Intercept for", opt_comp_final, "components:", intercept, "\n")

# Extract regression coefficients for optimal components
coefs <- coef(pls_final, ncomp = opt_comp_final)
coefs_vec <- as.vector(coefs)
names(coefs_vec) <- dimnames(coefs)[[1]]

# Show regression equation
#intercept <- coef(pls_final, ncomp = opt_comp_final)[,,1][1]
coefs <- coef(pls_final, ncomp = opt_comp_final)[,,1]
cat("PLS Regression Equation:\n")
cat("Ti[ICP-MS_pred] = ",  paste0(round(coefs, 4), "*", names(coefs), collapse = " + "), "\n") #intercept removed as NULL 

# Variable Importance in Projection (VIP) scores summarize each variable’s contribution to the PLS model.
# Function to calculate VIP scores
vip_scores <- function(pls_final, opt_comp_final) {
  W <- pls_final$loading.weights[, 1:opt_comp_final]
  T <- pls_final$scores[, 1:opt_comp_final]
  Q <- pls_final$Yloadings[1:opt_comp_final]
  
  SS <- colSums(T^2) * Q^2
  total_SS <- sum(SS)
  
  vip <- sqrt(ncol(W) * rowSums((W^2) * SS) / total_SS)
  names(vip) <- rownames(W)
  return(vip)
}

# Get VIP values
vip <- vip_scores(pls_final, opt_comp_final)
print(round(vip, 3))
# VIP > 1: Important — strong contribution to the model
# 0.8 < VIP < 1: Moderate contribution
# VIP < 0.8: Unimportant — likely not useful for prediction

# Write console summary info to file 
#sink(file = NULL)

# Convert to Tidyverse for calculations later ----------------------------------

# Generate coefficients table
coefs_final <- tibble(term = names(coefs),estimate = as.numeric(coefs)) %>% 
  pivot_wider(names_from = term, values_from = estimate)
coefs_final

# Generate equation text for ggplot
eq_text <- paste0("Ti[ICP-MS_predcited]:", paste0(round(coefs, 4), "*", names(coefs), collapse = " + "))
eq_text

# Create a tibble with site name and predicted output for ggplot
df_joined <- cbind(df, predicted_Ti)
df_predicted <- as_tibble(df_joined) %>% 
  rename(Ti_predicted = last_col())
df_predicted
write.csv(df_predicted,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/cps/Ti_ICP_predicted.csv", row.names = FALSE)

# Plot predictions vs actual, coloured by Site in ggplot -----------------------

theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final <- ggplot(df_predicted, aes(x = Ti_ICP, y = Ti_predicted, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted cps model",
    x = "Measured Ti ppm [ICP-MS]",
    y = "Predicted Ti ppm ICP-MS] (PLS cps Model)"
  ) +
  #theme_minimal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 4.5, y = 25000, label = RMSE_text, parse = TRUE, hjust = 0) + # add equation to graph
  annotate("text", x = 4.5, y = 26000, label = eq_text, parse = TRUE, hjust = 0) + # add equation to graph
  annotate("text", x = 4.5, y = 27000, label = R2_text, parse = TRUE, hjust = 0) + # add equation to graph
  coord_equal()
print(PLS_final)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/cps/Ti_plsr_final.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Ti as ppm conversion from original plsr model  -------------------------------

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

# Predict Ti using pls_final model and ITRAX ACE dataset as X_new value - use for ppm reconstruction
y_pred_scaled <- predict(pls_final, newdata = X_new, ncomp = 3)
y_pred_scaled_RMSE_upper <- y_pred_scaled + rmse
y_pred_scaled_RMSE_lower <- y_pred_scaled - rmse

# # Just for comparison - calculate decentered and descaled predicted Ti using mean and sd scale parameters of training Y dataset
# # Get original scale parameters from training Y - center = TRUE, scale = FALSE in plsr by not scaled by default
# y_mean <- mean(Y_test) # only apply if not scaled as here
# y_sd   <- sd(Y_test)
# # Descale & decenter then take exponential 
# y_pred_descaled <- y_pred_scaled * y_sd + y_mean # 
# y_pred_descaled_RMSE_upper <- y_pred_scaled + rmse
# y_pred_descaled_RMSE_lower <- y_pred_scaled - rmse

# Make into tibble for export and plotting & write to file
Ti_ppm_PLS1 <- tibble(Ti_pred_cps = y_pred_scaled,
                      Ti_ppm_PLS_cps = y_pred_scaled, 
                      Ti_upper_PLS_cps = y_pred_scaled_RMSE_upper,
                      Ti_lower_PLS_cps = y_pred_scaled_RMSE_lower)
Ti_ppm_PLS_predicted <- bind_cols(ACE_xrf_predict_Ti_PLS, Ti_ppm_PLS1)
write.csv(Ti_ppm_PLS_predicted,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/cps/ACE_Ti_ppm_PLS_predicted.csv", row.names = FALSE)

# Extract HER42PB data for Figure 4
HER42PB_Ti_ppm_PLS_predicted <- Ti_ppm_PLS_predicted %>% 
  filter(Site == "HER42PB")
HER42PB_Ti_ppm_PLS_predicted
write.csv(HER42PB_Ti_ppm_PLS_predicted,"Papers_R/2024_DeVleeschouwer/ACE_Calibration/ACE_PLS/Output/PLSR/cps/HER42PB_Ti_ppm_PLS_predicted_cps.csv", row.names = FALSE)

# END --------------------------------------------------------------------------

