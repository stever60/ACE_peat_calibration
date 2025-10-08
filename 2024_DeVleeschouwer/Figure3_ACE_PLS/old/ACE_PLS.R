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
packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'readr', 'ggpubr', 'patchwork',
              'gridExtra', 'cowplot', 'vegan', 'rioja', 'ellipse', 'factoextra',
              'reshape2', 'GGally', 'ggsci', 'ggdendro', 'dendextend', 'dynamicTreeCut',
              'colorspace', 'cluster', 'magrittr', 'mgcv', 'gtable', 'repr',
              'bestNormalize','sjmisc', 'chemometrics', 'compositions', 
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis',
              'ggrepel', 'itraxR', 'PeriodicTable', 'errors', 'forecast', 'broom',
              'directlabels', 'performance', 'lmtest', 'ggpmisc', 'cowplot', 'Hmisc', 
              'mdatools', 'pls', 'caret')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

# Define parameters ------------------------------------------------------------

# XRF-CS elements defined by ITRAX acf and matched to Francois ICPMS element list

acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")
# Scatter parameters
scatter_param <- c("Mo_inc", "Mo_coh", "inc_coh", "Ln_inc_coh", "coh_inc", "Ln_coh_inc", "Total_scatter",
                   "Mo_inc_sd", "Mo_coh_sd", "inc_coh_sd", "Ln_inc_coh_sd", "coh_inc_sd", "Ln_coh_inc_sd", "Total_scatter_sd")

acf_icp_Elements_key <- c("K", "Ca", "Ti", "Mn", "Fe", "Zn", "Rb", "Sr", "Zr", "Mo_coh")
acf_icp_Elements_key_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Zn_sd", 
                             "Rb_sd", "Sr_sd", "Zr_sd", "Mo_coh_sd")
acf_icp_Elements_key1 <- c("K", "Ca", "Ti", "Mn", "Fe", "Zn", "Rb", "Sr", "Zr", 
                           "Mo_inc", "Mo_coh") # Mo_inc included
acf_icp_Elements_key_sd1 <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Zn_sd", 
                              "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd","Mo_coh_sd")

# ICPMS - elements defined by Francois & by ITRAX acf
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", 
                      "Co_ICP", "Ni_ICP", "Cu_ICP", "Zn_ICP", "As_ICP", "Rb_ICP", 
                      "Sr_ICP", "Zr_ICP", "Pb_ICP", "Dry_mass")
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", 
                      "Ni_ICP", "Cu_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", 
                         "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# Key XRF and ICPMS elements
xrf_icp_elements <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                      "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                      "Zr", "Zr_ICP", "Mo_coh")
xrf_icp_elements1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                       "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP", "Mo_inc",  "Mo_coh") # Mo_inc included

xrf_icp_elements2 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                      "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                      "Zr", "Zr_ICP")

# MSCL
mscl_param  <- c("Den1_SAT", "MS1_SAT", "DCMS1_SAT", "Impedance_SAT", 
                 "Fract_Porosity_SAT", "Resistivity_SAT")

# Subsample parameters
subsample_param <- c("Water_Content", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density","Dry_density",	"Dry_density_err", 
                     "DMAR", "DMAR_err", "DMAR_err_pc")

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density")


# Set up PLS regression & load data -------------------------------------------------

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

# Select dataset to use
data <- ACE_matched_xrf_icp_log_inc
data

# Select predictors (ITRAX elements) and response (ICP-MS Ti)
df0 <- data %>%
  select(c(Site, depth, SH20_mean_age, SH20_mean_95CI, all_of(xrf_icp_elements1), inc_coh, coh_inc, Dry_mass)) %>%
  drop_na()
df0

# Convert tibble to data frame for base R
df <- as.data.frame(df0)

# Define predictor (ITRAX - multivariate) variables for calibration
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")



#-------------------------------------------------------------------------------

# 1. Partial Least Squares (PLS) regression with cross validation --------------

# If the dataset is small, consider using techniques like k-fold cross-validation to evaluate the model's performance without relying on a single test set. 
# Calibrate Titanium (Ti) measured by ICP-MS against other elements measured by ITRAX using pls package.
# Create a tibble with main 12 elements as predictors (ITRAX) and ICPMS equivalents as response (ICP-MS)

# Write summary outputs from console to txt file in Output folder
#sink(file = "Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/Summary_stats.txt")

# Split into predictors (X) using key elements only and response (Y) as Ti_ICP
X <- df[, key_elements_xrf]  # ITRAX elements
Y <- df$Ti_ICP  # Define response ICP-MS variable for calibration

# Perform PLS regression with LOO - leave one out validation - useful for unbiased prediction error assessment
pls_train_model_LOO <- plsr(Y ~ ., data = data.frame(X, Y), validation = "LOO") # Use scale = TRUE to perform as Z scores - do not use for calibration

# Get optimal number of components (lowest LOO error)
opt_comp_LOO <- which.min(RMSEP(pls_train_model_LOO)$val[1, , -1])
cat("Optimal number of components LOO:", opt_comp_LOO, "\n")

# Perform PLS regression with 10-fold cross-validation for model tuning - useful for small datasets
pls_train_model_CV <- plsr(Y ~ ., data = data.frame(X, Y), validation = "CV") # Use scale = TRUE to perform as Z scores - do not use for calibration

# Get optimal number of components (lowest CV error)
opt_comp_CV <- which.min(RMSEP(pls_train_model_CV)$val[1, , -1])
cat("Optimal number of components CV:", opt_comp_LOO, "\n")
# check CV vs LOO match (or not)

# Summary of the CV model output
summary(pls_train_model_CV)

# Plot root mean squared error of prediction (RMSEP) to choose number of components
plot(RMSEP(pls_train_model_CV), legendpos = "topright")

# Refit model with optimal number of components
pls_final <- plsr(Y ~ ., data = X, ncomp = opt_comp_CV)

# Predict Ti using the model
predicted_Ti <- predict(pls_final, ncomp = opt_comp_CV)

# Plot predicted vs actual draft output
plot(Y, predicted_Ti, xlab = "Measured Ti (ICPMS)", ylab = "Predicted Ti (from ITRAX)", main = "PLS Calibration")
abline(0, 1, col = "red")

# Calculate R2 and RMSE
r2 <- cor(Y, predicted_Ti)^2
rmse <- sqrt(mean((Y - predicted_Ti)^2))
cat("R² =", r2, "\nRMSE =", rmse, "\n")

# Coefficients of the final model
coefficients(pls_final)

# Extract regression coefficients for optimal components
coefs <- coef(pls_final, ncomp = opt_comp_CV)
coefs_vec <- as.vector(coefs)
names(coefs_vec) <- dimnames(coefs)[[1]]

# Show regression equation
#intercept <- coef(pls_final, ncomp = opt_comp_CV)[,,1][1]
coefs <- coef(pls_final, ncomp = opt_comp_CV)[,,1]
cat("PLS Regression Equation:\n")
cat("Ti_ICPMS = ",  paste0(round(coefs, 4), "*", names(coefs), collapse = " + "), "\n") #removed intercept, " + ",

# Write console summary info to file 
#sink(file = NULL)

# Create a tibble with site name and predicted output for ggplot
df_joined <- cbind(df, predicted_Ti)
df_predicted <- as_tibble(df_joined) %>% 
  rename(Ti_predicted = last_col())
write.csv(df_joined,"Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/Ti/Ti_ICP_predicted.csv", row.names = FALSE)

# Plot predictions vs actual, colored by Site in ggplot


ggplot(df_predicted, aes(x = Ti_ICP, y = Ti_predicted, color = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "PLS Regression: Measured vs Predicted Titanium (Ti)",
    x = "Measured Ti (ICP-MS)",
    y = "Predicted Ti (PLS Model)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()


# 2. PLS Train-Test Split tidyverse method - as above but with 70% training vs 30% test dataset ---------------------------------

# If the primary goal is to build a highly accurate predictive model, a larger training set might be preferred
# If the focus is on evaluating the model's generalizability, a larger testing set might be more appropriate. 

# define df_train as training dataset and elements to use for Ti_ICP response model
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

# Equations with intercept - unscaled only - doesnt work here because some data return NULL value
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


# 3. PLS Train-Test Split Base R & ggplot method - as above ---------------------------------

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


# ########## END ########## -----------------------------------------------












################################################################################
# Tidyverse solution - doesnt work  -------------------------------------

# Split into Training and Testing Sets (Optional)
set.seed(123)  # for random reproducibility if used elsewhere
split <- initial_split(df, prop = 0.8)
train <- training(split)
test  <- testing(split)

# Perform PLS Regression ------------------------------------------------

# Fit PLS model
pls_model <- plsr(Ti_ICPMS ~ ., data = train, scale = TRUE, validation = "CV") # or "LOO" Leave-one-out

# Summary of model
summary(pls_model)

# Plot RMSEP to choose number of components
plot(RMSEP(pls_model), legendpos = "topright")

# Evaluate Model Performance -------------------------------------------- 

# Predict on test set
predictions <- predict(pls_model, newdata = test, ncomp = 3)  # choose optimal components based on RMSEP

# Combine predictions and actual
results <- tibble(
  actual = test$Ti_ICPMS,
  predicted = as.vector(predictions)
)

# Plot predicted vs actual
ggplot(results, aes(x = actual, y = predicted)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "PLS Regression: Predicted vs Actual Ti (ICP-MS)",
       x = "Actual Ti (ICP-MS)", y = "Predicted Ti") +
  theme_minimal()

# Compute R-squared
cor(results$actual, results$predicted)^2

# Options 
# Use scale = TRUE to normalize predictors.
# Adjust ncomp based on cross-validation RMSEP plots.
# Include more elements from ITRAX as predictors if available.
# For many elements, consider variable selection or regularized PLS variants.






# Single element only 
# Below is an example of R code that performs Partial Least Squares Regression (PLSR) between elemental data for Titanium (Ti) measured by ICP-MS and XRF scanning.

# Assume your dataset contains two columns:

# Ti_ICPMS: Titanium measured by ICP-MS
# Ti_XRF: Titanium measured by XRF

  # Load necessary package
  install.packages("pls")  # Run only if 'pls' isn't already installed
library(pls)

# Simulated or actual dataset
# Replace this with your real data
# Example:
# data <- read.csv("your_data.csv")
# head(data)

# Example mock data
set.seed(123)
Ti_ICPMS <- rnorm(50, mean = 100, sd = 10)
Ti_XRF <- Ti_ICPMS + rnorm(50, sd = 5)
data <- data.frame(Ti_ICPMS, Ti_XRF)

# Fit PLS regression model
# We use Ti_XRF as predictor (X), and Ti_ICPMS as response (Y)
pls_model <- plsr(Ti_ICPMS ~ Ti_XRF, data = data, validation = "CV")

data(oliveoil)
mod <- plsr(sensory ~ chemical, ncomp = 4, data = oliveoil, validation = "LOO")
RMSEP(mod)
## Not run: plot(R2(mod))


# Summary of the model
summary(pls_model)

# Plot Root Mean Square Error of Prediction (RMSEP)
plot(RMSEP(pls_model), main = "PLSR - RMSEP")

# Predict using the model
predictions <- predict(pls_model, ncomp = 1)

# Plot actual vs predicted
plot(data$Ti_ICPMS, predictions, 
     xlab = "Observed Ti_ICPMS", 
     ylab = "Predicted Ti_ICPMS", 
     main = "Observed vs Predicted Ti (PLSR)")
abline(0, 1, col = "red")
# END