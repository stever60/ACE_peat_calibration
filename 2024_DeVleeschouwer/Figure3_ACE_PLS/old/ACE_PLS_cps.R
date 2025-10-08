# ACE PLS using pls package 

# XRF clr model and Ti_ICP reponse variable to clr dataset

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

# Import the matched ICPMS-XRF log (element/inc) dataset and create a PLS-model for prediction of Ti_ICP
# Using up to 11 other elements as predictors. 

# Import the matched xrf-icp dataset, and convert cps to clr - rename cps elements as element_clr

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

# 1. Partial Least Squares (PLS) regression ------------------------------------

# If the dataset is small, use eg k-fold cross-validation to evaluate the model's performance without relying on a single test set. 
# Calibrate Titanium (Ti) measured by ICP-MS against other elements measured by ITRAX using pls package.
# Create a tibble of main 12 elements as predictors (ITRAX) and ICPMS equivalents as response (ICP-MS)
# Add predicted output to tibble at end

# Define predictor ITRAX variables & ICPMS response variable for calibration model tests
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")

# Split into predictors (X) using key elements only and response (Y) as ICPMS element e.g., Ti_ICP
X_test <- df[, key_elements_xrf]  # ITRAX elements
Y_test <- df$Ti_ICP #define response ICP-MS variable for calibration

# Write summary outputs from console to txt file in Output folder instead of console
#sink(file = "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/Summary_stats.txt")

# LOO-PLSR (leave-one-out) test model to determine optimal number of comps -----

# Perform PLS regression with LOO - leave one out validation - useful for unbiased prediction error assessment
pls_train_model_LOO <- plsr(Y_test ~ ., data = data.frame(X_test, Y_test), validation = "LOO", jackknife = TRUE) # Use scale = TRUE to perform as Z scores - do not use for calibration

# Get optimal number of components (lowest LOO error)
opt_comp_LOO <- which.min(RMSEP(pls_train_model_LOO)$val[1, , -1])
cat("Optimal number of components LOO:", opt_comp_LOO, "\n")

# Summary of LOO model output and jackknife t tests of regression coefficients & std errors
summary(pls_train_model_LOO)
jack.test(pls_train_model_LOO, ncomp = opt_comp_LOO, use.mean = TRUE)

# CV-PLSR (cross-validated) test model to determine optimal number of comps ----

# Perform PLS regression with 10-fold cross-validation for model tuning - useful for small datasets
pls_train_model_CV <- plsr(Y_test ~ ., data = data.frame(X_test, Y_test), validation = "CV", jackknife = TRUE) # Use scale = TRUE to perform as Z scores - do not use for calibration

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
pls_final <- plsr(Y_final ~ ., data = X_final, ncomp = opt_comp_final, validation = "CV", jackknife = TRUE, center = TRUE)

# Summary of final plsr model and jackknife t tests of regression coefficients & std errors
summary(pls_final)
jack.test(pls_final, ncomp = opt_comp_final, use.mean = TRUE)

# Predict Ti using the model
predicted_Ti <- predict(pls_final, ncomp = opt_comp_final)

# Plot predicted vs actual draft output
plot(Y_final, predicted_Ti, xlab = "Measured Ti (ICPMS)", ylab = "Predicted Ti (from ITRAX)", main = "PLS Calibration")
abline(0, 1, col = "red")

# Calculate R2 and RMSE
r2 <- cor(Y_final, predicted_Ti)^2
rmse <- sqrt(mean((Y_final - predicted_Ti)^2))
cat("R² =", r2, "\nRMSE =", rmse, "\n")

# Generate equation text for ggplot
R2_text <- paste0("R^2: ", paste0(round(r2, 2)))
RMSE_text <- paste0("RMSEP: ", paste0(round(rmse, 4)))                  

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

# Write console summary info to file 
#sink(file = NULL)

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
write.csv(df_predicted,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/Ti_ICP_predicted.csv", row.names = FALSE)

# Plot predictions vs actual, coloured by Site in ggplot -----------------------

theme_set(theme_classic(base_size=12))
palette_set <- "jco" # define cb-friendly palette used in other plots or use "npg", "uchicago"
PLS_final <- ggplot(df_predicted, aes(x = Ti_ICP, y = Ti_predicted, color = Site)) +
  ggpubr::color_palette(palette_set) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Ln[Ti/inc)",
    x = "Measured Ln [Ti / inc] (ICP-MS)",
    y = "Predicted Ln [Ti / inc] (ICP-MS) (PLS Model)"
  ) +
  #theme_minimal() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  annotate("text", x = 4.5, y = 9.25, label = RMSE_text, parse = TRUE, hjust = 0) + # add equation to graph
  annotate("text", x = 4.5, y = 9.5, label = eq_text, parse = TRUE, hjust = 0) + # add equation to graph
  annotate("text", x = 4.5, y = 9.75, label = R2_text, parse = TRUE, hjust = 0) + # add equation to graph
  coord_equal()
print(PLS_final)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/Ti_plsr_final.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

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

# Ti ppm & flux conversions ----------------------------------------------------

# Import ACE Ti XRF-CS log_inc data and convert to ppm with +/- RMSE error
ACE_xrf_log_inc <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_log_inc.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>% 
  select(Location:MSE, all_of(main_elements_xrf), Total_scatter, inc_coh, coh_inc) %>%
  filter(qc == "TRUE") %>% 
  filter(validity == "1") %>% 
  mutate(across(all_of(main_elements_xrf), ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_log_inc

# PLS calibration: convert XRF-CS data to ppm from log_inc input ---------------

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
write.csv(ACE_xrf_conversion_Ti_PLS,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/ACE_plsr_xrf_Ti_converted_PLS.csv", row.names = FALSE)

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

