# -------------------------------------------------------------------------
# Regression Comparison Multivariate function
# Matched XRF-CS log_inc & log ICP-MS matched dataset
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# ReadMe description README.txt ----------------------------------------

# This function performs a complete calibration, validation, and reconstruction 
# workflow for multiple geochemical elements using log transformed ICP-MS 
# calibration data and log transformed XRF datasets. For each element, it fits 
# eight univariate and, where applicable, multivariate statistical models 
# (OLS, WLS, weighted OLS, weighted WLS, Bayesian Generlised Linear Model (BGLM) 
# regression, Random Forest, PLS with LOO, and PLS with k-fold), evaluates model 
# performance using R², RMSE, 10-fold RMSEP, AIC/BIC (where applicable), and runs 
# diagnostic tests for normality and heteroscedasticity. It produces predicted 
# vs. observed plots, residual and influence diagnostics, and saves model 
# rankings based on R² → RMSE → RMSEP. The function then generates ppm-scale 
# predictions with 95% confidence intervals for new XRF data, produces site-level 
# and multi-site and multi-element comparison plots for depth and age, and outputs 
# both 8 different (unweighted, weighted OLS, WLS) and multivariate (2 PLS models, 
# RF, BGLM) models and best model-multivariate-only visualisations. Finally, it 
# compiles global summary tables (in log space and ppm space), saves per-element 
# prediction files, and organises all outputs into structured folders for 
# streamlined interpretation and reporting.

# It nows includes Random Forest OOB (Out-of-Bag) - a method for estimating a 
# random forest model's performance without using a separate test set. RF works
# by having each tree in the forest evaluate the data points that were not 
# included in its own training sample (the "out-of-bag" data). By averaging these 
# predictions across all trees and their respective OOB samples, the out-of-bag 
# error provides a robust and unbiased estimate of the model's generalization error.

# This function fits: 
# • OLS, WLS (var), weighted OLS (sd), 2-step WLS, Bayes, RF, PLS LOO & k-fold
# • R2 (true) as Sum of Squares = 1-(SSres / SS_tot)
# • writes per-element folders into your chosen directory
# Saves:
# • Residual distribution plots (PDF)
# • Predicted vs Observed plots (PDF)
# • PLS RMSEP validation plots (PDF)
# • Influence plots (Cook’s distance + residuals vs leverage) as PDF
# Writes a per-element diagnostics to a .csv file with:
# • AIC & BIC tests
# • Normality test (Shapiro–Wilk)
# • Heteroscedasticity (Breusch–Pagan) for LM/Bayes
# • Unified 10-fold CV RMSE for every model (where possible)
# Native CV metrics for:
# • RF OOB MSE & R² 
# • PLS RMSEP at best ncomp
# •	Explanations for each test in plain English
# Makes a Best Model performance summary
# • ranking based on R2, RMSE, RMSEP where AIC & BIC dont exist (PLS, RF)
# • where R2 higher, RMSE/RMSEP lower, AIC, BIC lower is better
# — sequential R that works on a Mac.

# Key features for each element:
# • OLS, WLS, weighted OLS, 2-step weighted WLS - based on existing code (Steve Roberts)
# • Bayesian GLM - using the bayesglmglm function in the arm.R package
# • Random Forest - using the randomForest.R package
# • PLS (LOO and 10-fold CV) - using the pls.r package 
# • Residual plots, normality tests, heteroscedasticity tests, influence plots
# • Predicted vs observed plots colored by Site using a JCO palette
# • RF variable importance, PLS RMSEP validation plots
# • Per-model/per-element diagnostics (PDF + CSV)
# • Performance summary CSV (R², RMSE, RMSEP, AIC, BIC)
# • Best-model table
# • README.txt file explaining criteria and rankings

# Multivariate run set up
# 
# Single-predictor models (based on existing manual code):
# •	OLS
# •	WLS (variance weighted OLS)
# •	OLS_weighted
# •	WLS_weighted
# 
# Multi-predictor models (updated from univariate function):
# •	Bayes → uses ALL 6 elements as predictors
# •	PLS_LOO → uses ALL 6 elements as predictors
# •	PLS_kfold → uses ALL 6 elements as predictors
# •	Random Forest → uses ALL 6 elements as predictors - NEW in v. 1.1
# 
# Ranking Logic
# 
# Models ranked by:
# 1.	Highest R² (best explanatory fit)
# 2.	Lowest RMSEP (best predictive CV performance)
# 3.	Lowest RMSE
# 
# No PLS priority. No AIC/BIC in ranking.
# 
# Directory Output Structure Per Element
# *_Model_Ranking.csv (R² → RMSEP → RMSE ranking)
# *_Diagnostics_Summary.csv
# *_PredVsObs.pdf
# *_Residuals.pdf
# *_Influence.pdf (LM models only)
# *_PLS_RMSEP.pdf
# README.txt
# 
# Global outputs in All_models & Best_models folders:
# •	AllElements_ModelSummary.csv
# •	BestModels_PerElement.csv

# Acknowledgements
# Based on original for OLS, WLS, PLS regression code by Steve Roberts
# Refined and made into a function with assistance from ChatGPT v5.1
# We acknowledge the use of ChatGPT (OpenAI, version 5.1) for assistance in 
# constructing functions based on manual code written by the authors and 
# generating some preliminary code The final implementations of the code and 
# functions were produced and verified independently by the authors.

# Improvements in v1.2-1.4  ---------------------------------------------------

# v. 1.4 
# Optional sections 13-18 applying the calibration to the XRF data in ACE_dataset
# to give calibration predictions and pred vs obs plots in ppm 
# repeats section 1-12 structure and diagnostic tests - but no age or depth plots

# v. 1.3-1.3.2
# - adding details to Site age and depth plots eg median, n. = values etc.

# v. 1.2.2
# - Correcting Section 12 age plot to one column only
# - Making new age plot with global min-max y axis ppm values 

# v. 1.2.3
# - chnaged pred vs obs back to jco site colours
# - made best model age comparison equal y axis log scale
# - added median ICP-MS line to age axis plots for each element

# v. 1.2.1
# Structural Improvements
# •	Function now requires the user to specify save_dir
# → No default hard-coded path; users choose their output directory in the run call.
# •	Ensures folder hierarchy is created inside a top-level folder
# → e.g., /Users/.../Regression_multivariate/ automatically gets:
#   
#   All_models_ppm/
#   Best_models_ppm/
#   Sites/
#   elements/
#   <element>/
#   
# •	Added library(rlang) as requested.
# •	Model predictions never cause models to be skipped anymore
# → If prediction fails, NA is inserted into output columns, but the model is still included in:
#     •	rankings
#     •	summaries
#     •	diagnostics
#     •	plots (except where no numeric prediction exists)
# 
# ⸻
# 
# Plotting Improvements
# 
# 1. Dynamic Layouts
# •	Plots no longer assume 6 elements; they now choose grid layout based on length(elements).
# 
# 2. Site-level plots restored
# •	All site-level depth-based and age-based plots now correctly appear in:
#   •	All_models
# •	Best_models
# •	Under Sites/ folder
# (This had been broken earlier but is now fixed.)
# 
# 3. Element-level comparison pages redesigned
# •	Each element now generates multi-site age–concentration comparison pages:
#   •	All models
# •	Multivariate-only models
# 
# 4. Forced single-column layout
# •	These element-level comparison plots now use 1-column vertical stacking so:
#     •	All sites share the same horizontal age axis
#     •	Between-site comparison becomes visually direct
# 
# 5. Added second set of plots with matched Y-ranges
# 
# For each element, we now generate:
#   
# Page Type	Y-axis	Filename suffix
# All models	free	_all.pdf
# All models	equalised	_all_y_equal.pdf
# Multivariate only	free	_multivariate.pdf
# Multivariate only	equalised	_multivariate_y_equal.pdf
# 
# The equalised versions compute global min–max ppm ranges across:
# •	All predicted values
# •	All ICP-MS observations
# 
# Then apply consistent ylim() across all site panels.
# 
# ⸻
# 
# Modelling Improvements
# •	Prediction intervals now computed for all models when possible, and fall back to RMSE-based bounds when not.
# •	Cross-validation RMSE no longer stops execution if one fold errors; it simply returns NA gracefully.
# •	All models are evaluated even when predictions contain NA.
# •	Diagnostic tables improved in clarity and formatting.
# 
# ⸻
# 
# Output Improvements
# •	AllElements_ModelSummary_ppm.csv
# •	PPM-scale error metrics (absolute error, multiplicative factors)
# •	Mean ppm per element
# •	Unified global model ranking table
# •	Per-element README files updated.
# 
# ⸻
# 
# Bug Fixes
# •	Fixed critical select() masking error (MASS::select vs dplyr::select)
#     → All now explicitly use dplyr::select.
# •	Fixed plotting errors and blank pages caused by empty objects.
# •	Fixed missing “Sites” folder and missing site-based plots.
# •	Fixed earlier crash due to invalid -suppressPackageStartupMessages.
# •	Repaired the collapsed plots caused by patchwork misuse.
# •	Corrected all spelling mistakes (e.g., ACE_dataset).
# 
# ⸻
# 
# Overall Result
# 
# A robust, flexible, user-friendly, fail-safe multivariate regression engine with:
# •	Full OLS/WLS/Bayes/RF/PLS support
# •	Clean dynamic layouts
# •	Multi-site + multi-model comparison pages
# •	Depth-based & age-based analysis
# •	Back-transformed predictions with intervals
# •	Extensive QA diagnostics
# •	Fully parallelisable architecture
# •	No crashes, no dropped models

# Sections ----------------------------------------------------------------

# Section 1–7   (setup, checks, helpers)
# Section 8     (calibration element loop)
# Section 9     (site-level depth plots)
# Section 10    (element-level age–conc comparisons)
# Section 11    (global calibration summary)
# Section 11.1  (ppm-scale calibration summary)
# Section 12    (finish message for calibration)
# Section 13    (ACE prediction)
# Section 14    (ACE diagnostic comparison vs ICP)
# Section 15    (ACE diagnostic plots)
# Section 17    (ACE global summary)
# Section 18    (ACE ppm-scale summary)
# Section 21 (optional correlation diagnostics)

# FLOW DIAGRAM OF FULL PROCESSING PIPELINE -------------------------------------

# Initial set up 
# 
# Start
# │
# ├─ Load packages
# │
# ├─ Create root output folder:  Regression_Multivariate
# │
# ├─ Create subfolders:
#   │   ├── All_models_ppm
# │   └── Best_models_ppm
# │
# └─ Prepare model labels, colours, plotting themes
# 
# Input data 
# 
# Inputs:
# │   • Calibration dataset (log-space ICP + XRF)
# │   • ICP_ppm (ppm-scale ICP observational data)
# │   • XRF_new (optional)
# │   • ACE_dataset (log-space XRF for prediction)
# │
# ├─ Check required columns
# │
# └─ Factorise Site variable
# 
# PER-ELEMENT LOOP (Main calibration engine)
# 
# For each element:
#   │
# ├─ Create element folder under All_models_ppm/<element>
#   │
# ├─ Build formulas:
# │     • Single-input: y ~ x
# │     • Multi-input : y ~ all elements
# │
# ├─ Fit 8 models:
# │     1. OLS
# │     2. WLS
# │     3. OLS_weighted
# │     4. WLS_weighted
# │     5. Bayes
# │     6. Random Forest
# │     7. PLS (LOO)
# │     8. PLS (k-fold = 10)
# │
# ├─ Compute validation metrics:
# │     • R²
# │     • RMSE (log-space)
# │     • RMSEP (10-fold CV)
# │     • AIC / BIC (for lm/glm)
# │
# ├─ Diagnostics:
# │     • Shapiro–Wilk
# │     • Breusch–Pagan
# │
# ├─ Generate plots:
# │     • Predicted vs Observed
# │     • Residual histograms
# │     • Influence plots (lm only)
# │     • PLS RMSEP curves
# │
# ├─ Save per-model diagnostics + CSV summary
# │
# └─ Collect global summary stats
# 
# SITE-LEVEL VISUALISATION (Depth & Age Profiles) - from XRF_new inputted as XRF_pred
# 
# For each Site:
#   │
# ├─ Generate depth profiles:
# │     • All models
# │     • Best multivariate models
# │     • ICP overlay (no transparency)
# │     • No grey ribbon for *_age_comparison_all
# │     • Grey ribbon kept for *_age_comparison_all_y_equal
# │
# ├─ Generate age profiles:
# │     • All models (with/without grey CI depending on Y-equal)
# │     • Best multivariate models
# │     • Add red “.” if PLS_k CI exceeds 2×ICP_sd (upper only)
# │
# └─ Save A4 multi-panel PDFs
# 
# For each Site:
#   │
# ├─ Generate depth profiles:
# │     • All models
# │     • Best multivariate models
# │     • ICP overlay (no transparency)
# │     • No grey ribbon for *_age_comparison_all
# │     • Grey ribbon kept for *_age_comparison_all_y_equal
# │
# ├─ Generate age profiles:
# │     • All models (with/without grey CI depending on Y-equal)
# │     • Best multivariate models
# │     • Add red “.” if PLS_k CI exceeds 2×ICP_sd (upper only)
# │
# └─ Save A4 multi-panel PDFs
# 
# ELEMENT-LEVEL MULTI-SITE AGE PROFILES
# 
# For each element:
#   │
# ├─ Combine all site profiles vertically (1 column)
# ├─ Y-equal version with shared log10 scale
# ├─ Add exceedance red “.” markers for PLS_k CI > 2×ICP_SD
# ├─ Remove all pink error bars (Section 10 patch)
# │
# └─ Save PDFs to:
# • All_models_ppm/elements
# • Best_models_ppm/elements
# 
# GLOBAL SUMMARY TABLES (Calibration)
# 
# Compute for each element:
# │    • Rank of all 8 models (R² → RMSEP → RMSE)
# │    • Best model per element
# │    • ppm-scale metrics (RMSE_abs_ppm, RMSE_factor …)
# │
# Save:
# │    • AllElements_ModelSummary.csv
# │    • BestModels_PerElement.csv
# │    • AllElements_ModelSummary_ppm.csv
# 
# CALIBRATION_pred FOLDER GENERATION
# 
# Create:
#   All_models_ppm/Calibration_pred/<element>/
#   │
# For each element:
#   │    ├─ Extract ICP_obs data (ppm)
# │    ├─ Extract ACE_dataset XRF predictors (log)
# │    ├─ Join metadata (Site, depth, SH20_age)
# │    └─ Save base file: ICP_pred_<element>.csv
# 
# APPLY MODELS TO ACE_dataset Calibration dataset as XRF_new (ppm predictions)
# 
# For each element:
#   │
# ├─ Load ICP_pred_<element>.csv
# │
# ├─ For each of 8 models:
# │     • Predict log-space ACE_dataset
# │     • Compute 95% CI (lm intervals or RMSE-based)
# │     • Backtransform to ppm
# │     • Add columns:
# │           <model>_<element>_Pred_ppm
# │           <model>_<element>_L95_ppm
# │           <model>_<element>_U95_ppm
# │
# └─ Save updated ICP_pred_<element>.csv
# 
# DIAGNOSTIC PLOTS FOR ACE PREDICTIONS (Calibration_pred)
# 
# For each element:
#   │
# ├─ Pred vs Obs (ACE vs ICP)
# ├─ Residual histograms
# ├─ Influence plots (lm only)
# ├─ PLS RMSEP
# |
#   └─ Save into:
#   All_models_ppm/Calibration_pred/<element>/diagnostics/
#   
#   GLOBAL SUMMARY TABLES (for Calibration_pred)
# 
# Compute:
# │    • R², RMSE, RMSEP rankings for ACE predictions
# │    • Best ACE model per element
# │
# Save:
# │    • AllElements_ModelSummary_pred.csv
# │    • BestModels_PerElement_pred.csv
# 
# ppm-SCALE SUMMARY  (for Calibration_pred)
# 
# Compute ppm-scale RMSE:
#   │    RMSE_abs_ppm = RMSE_log × mean_ppm
# │
# Save:
#   │    AllElements_ModelSummary_ppm_pred.csv
# 
# MULTI-ELEMENT CORRELATION DIAGNOSTICS
# 
# Compute & save:
#   │
# ├─ Calibration ICP log-space correlations
# │
# ├─ ACE_dataset XRF predictor correlations
# │
# ├─ ACE predicted ppm correlations (best model per element)
# │
# ├─ Outputs:
#   │     Pearson / Spearman CSVs
# │     Heatmaps (PDF)
# │     Scatter matrix (PDF)
# │
# └─ Folder:
#   All_models_ppm/Correlation_Diagnostics/

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Clear all ------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Install libraries needed - do once ---------------------------------------------------------------
install.packages(c(
  "broom", "performance", "see", "ggplot2", "dplyr", "purrr",
  "boot", "psych", "lmtest", "arm", "randomForest", "pls", "ggpubr", 
  "qqplotr", "future", "progressr"
))


# Load libraries used  --------------------------------------------------
packages <-c(
  "broom", "performance", "see", "ggplot2", "dplyr", "purrr",
  "boot", "psych", "lmtest", "arm", "randomForest", "pls", "ggpubr", 
  "qqplotr", "ggsci", "car", "performance", "ggpmisc", "future", "progressr"
)
lapply(packages, library, character.only=TRUE)

# Define elements to use - not really needed but matches other definitions -------------------------------------------------------

# ICP elements of interest defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "Dry_mass")

# XRF-CS acf elements matched to Francois ICPMS elements
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_only_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                               "Zn", "Rb", "Sr", "Zr")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")
acf_icp_Elements_only_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                                  "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd")

# Scatter parameters
scatter_param <- c("Mo_inc", "Mo_coh", "inc_coh", "Ln_inc_coh", "coh_inc", "Ln_coh_inc", "Total_scatter",
                   "Mo_inc_sd", "Mo_coh_sd", "inc_coh_sd", "Ln_inc_coh_sd", "coh_inc_sd", "Ln_coh_inc_sd", "Total_scatter_sd")

# ICPMS elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# key elements to simplify plotting 
xrf_icp_Elements_key <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP",
                          "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "Dry_mass")

# key elements_reduced for more simplified plots 
xrf_icp_Elements_key_reduced <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")

# key elements for Figure 2a
xrf_icp_Elements_key_Fig2 <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP", 
                               "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "Dry_mass")

# Subsample parameters
subsample_param <- c("Water_Content", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density","Dry_density",	"Dry_density_err", 
                     "DMAR", "DMAR_err", "DMAR_err_pc")

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density")

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Import datasets
# • ACE matched XRF-CS log (element/inc) file for calibration ------------------
ACE_dataset <- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/ACE_Regression/Input/ACE_subsample_xrf_icp_matched_log_inc.csv",
  name_repair = "minimal") %>%   # equivalent to check.names = FALSE 
  rename (SH20_age = SH20_mean_age) %>% 
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>% # Replace infinite values with NA
  #mutate(Site = as.factor(Site)) %>%  # Convert Site to factor
  print()
str(ACE_dataset)
names(ACE_dataset)
# Check dataset is loaded OK & working before running function
c("Ti", "Ti_ICP") %in% names(ACE_dataset)
c("Ca", "Ca_ICP") %in% names(ACE_dataset)
c("Fe", "Fe_ICP") %in% names(ACE_dataset)
summary(ACE_dataset$Ti)
var(ACE_dataset$Ti)
summary(ACE_dataset$Ti_ICP)
var(ACE_dataset$Ti_ICP)
cor(ACE_dataset$Ti, ACE_dataset$Ti_ICP)

# • XRF_pred for new ppm predictions from calibration ----------------------

main_elements_xrf <- c("Ca", "Fe", "Mn", "Sr", "Ti", "Zr")

XRF_pred <- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/ACE_Regression/Input/ACE_ITRAX_qc_acf_log_inc.csv"
) %>%
  dplyr::select(Site, depth, SH20_age, dplyr::all_of(main_elements_xrf)) %>%
  dplyr::filter(Site %in% c("BI10","HER42PB","KER1","KER3","PB1")) %>%
  dplyr::mutate(across(all_of(main_elements_xrf), ~ ifelse(. <= -10, NA, .))) %>% 
  print()
str(XRF_pred)
names(ACE_dataset)

# • ICP_obs for plotting measured ICP-MS and error data  ---------------
ICP_obs<- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/ACE_Regression/Input/ACE_subsample_icp_xrf_matched_cps.csv",
  name_repair = "minimal") %>%   # equivalent to check.names = FALSE 
  rename (SH20_age = SH20_mean_age) %>% 
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>% # Replace infinite values with NA
  #mutate(Site = as.factor(Site)) %>%  # Convert Site to factor
  print()
str(ICP_obs)
names(ICP_obs)

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Regression comparison function & application to new XRF data --------------------

run_full_regressions <- function(
    elements,           # e.g. c("Ti","Fe","Rb","Sr","Zr","Ca")
    data,               # calibration dataset (log-space ICP etc.)
    XRF_new = NULL,     # new XRF dataset (predictors in log-space)
    ICP_ppm,            # observed ICP in ppm
    save_dir,           # base directory
    ACE_predict = TRUE, # for sections 13-21 to run 
    var_suffix      = "_var",
    sd_suffix       = "_ICP_sd",
    summary_all_csv = "AllElements_ModelSummary.csv",
    best_models_csv  = "BestModels_PerElement.csv"
) {
  
  # ===============================================================
  # 1. LOAD PACKAGES
  # ===============================================================
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(randomForest)
    library(pls)
    library(arm)
    library(ggsci)
    library(progress)
    library(lmtest)
    library(tidyr)
    library(patchwork)
    library(rlang)
    library(scales)
  })
  
  # ===============================================================
  # 2. OUTPUT FOLDERS (ALWAYS UNDER "Regression_Multivariate")
  # ===============================================================
  out_dir <- file.path(save_dir, "Regression_Multivariate")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_base  <- file.path(out_dir, "All_models")
  best_base <- file.path(out_dir, "Best_models")
  dir.create(all_base,  recursive = TRUE, showWarnings = FALSE)
  dir.create(best_base, recursive = TRUE, showWarnings = FALSE)
  
  # NEW: calibration outputs moved into subfolder
  calib_base <- file.path(all_base, "Calibration_outputs")
  dir.create(calib_base, recursive = TRUE, showWarnings = FALSE)
  
  # ===============================================================
  # 3. SITE FACTOR
  # ===============================================================
  if (!"Site" %in% names(data))
    stop("Calibration dataset must contain a 'Site' column.")
  
  data$Site <- factor(
    data$Site,
    levels = c("BI10","HER42PB","KER1","KER3","PB1")
  )
  
  # ===============================================================
  # 4. THEME & MODEL LABELS
  # ===============================================================
  theme_small <- theme(
    text       = element_text(size = 8),
    axis.text  = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.title   = element_text(size = 8, face = "bold")
  )
  
  model_levels <- c(
    "OLS","WLS",
    "OLS_weighted","WLS_weighted",
    "Bayes","RF","PLS_LOO","PLS_kfold"
  )
  
  model_labels <- c(
    "OLS","WLS",
    "OLS(wt)","WLS(wt)",
    "Bayes","RF","PLS(LOO)","PLS(k)"
  )
  
  cols_models <- ggsci::pal_npg("nrc")(length(model_labels))
  names(cols_models) <- model_labels
  
  # Force specific colours
  cols_models["PLS(k)"] <- "salmon3"
  cols_models["Bayes"]  <- "gold2"
  
  cols_all <- c(cols_models, "ICP-MS (Observed)" = "blue")
  
  # ===============================================================
  # 5. DATA CHECK
  # ===============================================================
  check_dataset_integrity <- function(data, el, y_col, x_col, var_col, sd_col) {
    errs <- c()
    
    req <- c(y_col, x_col)
    for (v in req)
      if (!v %in% names(data))
        errs <- c(errs, paste("Missing:", v))
    
    if (any(is.na(data[[y_col]]))) errs <- c(errs, paste("NA in", y_col))
    if (any(is.na(data[[x_col]]))) errs <- c(errs, paste("NA in", x_col))
    
    if (!is.numeric(data[[y_col]])) errs <- c(errs, paste(y_col, "not numeric"))
    if (!is.numeric(data[[x_col]])) errs <- c(errs, paste(x_col, "not numeric"))
    
    if (var(data[[x_col]], na.rm = TRUE) == 0)
      errs <- c(errs, "Predictor has zero variance")
    
    if (sd_col %in% names(data))
      if (any(data[[sd_col]] <= 0, na.rm = TRUE))
        errs <- c(errs, paste(sd_col, "contains non-positive values"))
    
    if (var_col %in% names(data))
      if (any(data[[var_col]] <= 0, na.rm = TRUE))
        errs <- c(errs, paste(var_col, "contains non-positive values"))
    
    errs
  }
  
  # ===============================================================
  # 6. 10-FOLD CV RMSE HELPER
  # ===============================================================
  cv_rmse_for_model <- function(
    mname, formula, data, y_col,
    var_col, sd_col, best_ncomp = NA,
    k = 10
  ) {
    set.seed(123)
    n <- nrow(data)
    if (n < k) return(NA_real_)
    
    folds <- sample(rep(1:k, length.out = n))
    cv_sse <- 0
    cv_n   <- 0
    
    for (fold in 1:k) {
      train <- data[folds != fold, , drop = FALSE]
      test  <- data[folds == fold, , drop = FALSE]
      
      fit <- tryCatch({
        if (mname == "OLS") {
          lm(formula, data = train)
        } else if (mname == "OLS_weighted" && sd_col %in% names(train)) {
          lm(formula, data = train, weights = 1/(train[[sd_col]]^2))
        } else if (mname == "WLS_weighted" && sd_col %in% names(train)) {
          init <- lm(formula, data = train, weights = 1/(train[[sd_col]]^2))
          rm   <- lm(abs(init$residuals) ~ init$fitted.values)
          w2   <- 1/(rm$fitted.values^2)
          lm(formula, data = train, weights = w2)
        } else if (mname == "WLS") {
          ols <- lm(formula, data = train)
          rm  <- lm(abs(ols$residuals) ~ ols$fitted.values)
          w   <- 1/(rm$fitted.values^2)
          lm(formula, data = train, weights = w)
        } else if (mname == "Bayes") {
          bayesglm(formula, data = train)
        } else if (mname == "RF") {
          randomForest(formula, data = train, ntree = 300)
        } else if (mname %in% c("PLS_LOO","PLS_kfold")) {
          plsr(formula, data = train, scale = TRUE, validation = "none")
        } else NULL
      }, error = function(e) NULL)
      
      if (is.null(fit)) next
      
      pred <- tryCatch({
        if (inherits(fit, "mvr"))
          as.numeric(predict(fit, ncomp = best_ncomp, newdata = test))
        else
          predict(fit, newdata = test)
      }, error = function(e) NULL)
      
      if (is.null(pred)) next
      
      y_test <- test[[y_col]]
      ok <- !is.na(y_test) & !is.na(pred)
      if (!any(ok)) next
      
      cv_sse <- cv_sse + sum((y_test[ok] - pred[ok])^2)
      cv_n   <- cv_n   + sum(ok)
    }
    
    if (cv_n == 0) return(NA_real_)
    sqrt(cv_sse / cv_n)
  }
  
  # ===============================================================
  # 7. EXPLANATION TEXT & GLOBAL STORAGE
  # ===============================================================
  expl_rank   <- "Ranking: R² highest → lowest; RMSEP lowest → highest; RMSE lowest → highest."
  expl_shapiro <- "Shapiro–Wilk: p > 0.05 desirable."
  expl_bp      <- "Breusch–Pagan: p < 0.05 indicates heteroscedasticity."
  expl_cv10    <- "10-fold CV RMSE (log-space)."
  expl_r2      <- "R² = 1 − SS_res/SS_tot (log-space)."
  expl_rmse    <- "RMSE = calibration error (log-space)."
  expl_rmsep   <- "RMSEP = 10-fold CV error (log-space)."
  
  global_summary <- tibble(
    Element = character(),
    Model   = character(),
    R2      = numeric(),
    RMSE    = numeric(),
    RMSEP   = numeric(),
    AIC     = numeric(),
    BIC     = numeric()
  )
  
  preds_store <- list()
  
  # ===============================================================
  # 8. ELEMENT LOOP
  # ===============================================================
  pb_elements <- progress_bar$new(
    format = "Processing Elements [:bar] :current/:total (:percent) eta::eta",
    total  = length(elements),
    clear  = FALSE,
    width  = 80
  )
  
  # GLOBAL STORAGE FOR LATER SECTIONS (ACE predictions, depth/age plots)
  model_store      <- list()
  calibration_rmse <- list()
  preds_store      <- list()
  
  for (el in elements) {
    
    pb_elements$tick()
    
    y_col  <- paste0(el, "_ICP")
    x_col  <- el
    var_col <- paste0(el, var_suffix)
    sd_col  <- paste0(el, sd_suffix)
    
    # OLD LOCATION:
    # element_dir <- file.path(all_base, el)
    # NEW LOCATION (requested change):
    element_dir <- file.path(calib_base, el)
    dir.create(element_dir, recursive = TRUE, showWarnings = FALSE)
    
    errors <- check_dataset_integrity(data, el, y_col, x_col, var_col, sd_col)
    if (length(errors) > 0)
      warning(paste("Issues in", el, ":", paste(errors, collapse = "; ")))
    
    # Formulas
    fmla_single <- as.formula(paste0(y_col, " ~ ", x_col))
    fmla_multi  <- as.formula(paste0(y_col, " ~ ", paste(elements, collapse = " + ")))
    
    y <- data[[y_col]]
    
    # ---------------------------------------------------------------
    # 8.1 FIT MODELS (never skip; failures -> NULL)
    # ---------------------------------------------------------------
    model_list <- list(
      OLS = tryCatch(lm(fmla_single, data = data), error = function(e) NULL),
      
      WLS = tryCatch({
        ols <- lm(fmla_single, data = data)
        rm  <- lm(abs(ols$residuals) ~ ols$fitted.values)
        w   <- 1/(rm$fitted.values^2)
        lm(fmla_single, data = data, weights = w)
      }, error = function(e) NULL),
      
      OLS_weighted = tryCatch({
        if (sd_col %in% names(data))
          lm(fmla_single, data = data, weights = 1/(data[[sd_col]]^2))
        else NULL
      }, error = function(e) NULL),
      
      WLS_weighted = tryCatch({
        if (sd_col %in% names(data)) {
          init <- lm(fmla_single, data = data, weights = 1/(data[[sd_col]]^2))
          rm   <- lm(abs(init$residuals) ~ init$fitted.values)
          w2   <- 1/(rm$fitted.values^2)
          lm(fmla_single, data = data, weights = w2)
        } else NULL
      }, error = function(e) NULL),
      
      Bayes = tryCatch(bayesglm(fmla_multi, data = data), error = function(e) NULL),
      
      RF = tryCatch(
        randomForest(fmla_multi, data = data, ntree = 500, importance = TRUE),
        error = function(e) NULL
      ),
      
      PLS_LOO = tryCatch({
        pls <- plsr(fmla_multi, data = data, scale = TRUE, validation = "LOO")
        best <- selectNcomp(pls, method = "onesigma", plot = FALSE)
        attr(pls, "best_ncomp") <- best
        pls
      }, error = function(e) NULL),
      
      PLS_kfold = tryCatch({
        pls <- plsr(fmla_multi, data = data, scale = TRUE,
                    validation = "CV", segments = 10)
        best <- selectNcomp(pls, method = "onesigma", plot = FALSE)
        attr(pls, "best_ncomp") <- best
        pls
      }, error = function(e) NULL)
    )
    
    # ---------------------------------------------------------------
    # 8.2 PER-ELEMENT TABLES & PRED DF
    # ---------------------------------------------------------------
    ranking_df <- tibble(
      Model = character(),
      R2    = numeric(),
      RMSE  = numeric(),
      RMSEP = numeric(),
      AIC   = numeric(),
      BIC   = numeric()
    )
    
    diag_df <- tibble(
      Model             = character(),
      Shapiro_W         = numeric(),
      Shapiro_p         = numeric(),
      Shapiro_expl      = character(),
      BreuschPagan_stat = numeric(),
      BreuschPagan_p    = numeric(),
      BreuschPagan_expl = character(),
      CV10_RMSE         = numeric(),
      CV10_expl         = character(),
      AIC               = numeric(),
      AIC_expl          = character(),
      BIC               = numeric(),
      BIC_expl          = character()
    )
    
    do_predictions <- !is.null(XRF_new)
    preds_df <- if (do_predictions) XRF_new else NULL
    
    pb_models <- progress_bar$new(
      format = paste0("  Models for ", el, " [:bar] :current/:total (:percent)"),
      total  = length(model_list),
      clear  = FALSE,
      width  = 70
    )
    
    # ===============================================================
    # 8.3 MODEL LOOP
    # ===============================================================
    for (mname in names(model_list)) {
      
      pb_models$tick()
      m <- model_list[[mname]]
      
      # Predictions on calibration data
      pred <- tryCatch({
        if (inherits(m, "mvr"))
          as.numeric(predict(m, ncomp = attr(m, "best_ncomp"), newdata = data))
        else if (!is.null(m))
          as.numeric(predict(m))
        else
          rep(NA_real_, nrow(data))
      }, error = function(e) rep(NA_real_, nrow(data)))
      
      resid_vals <- y - pred
      
      # Metrics
      R2 <- tryCatch({
        SS_res <- sum((y - pred)^2, na.rm = TRUE)
        SS_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
        1 - SS_res / SS_tot
      }, error = function(e) NA_real_)
      
      RMSE <- tryCatch(
        sqrt(mean((y - pred)^2, na.rm = TRUE)),
        error = function(e) NA_real_
      )
      
      cv_formula <- if (mname %in% c("OLS","WLS","OLS_weighted","WLS_weighted"))
        fmla_single else fmla_multi
      
      RMSEP_val <- tryCatch(
        cv_rmse_for_model(
          mname, cv_formula, data, y_col,
          var_col, sd_col,
          best_ncomp = if (inherits(m, "mvr")) attr(m, "best_ncomp") else NA
        ),
        error = function(e) NA_real_
      )
      
      AIC_val <- if (inherits(m, c("lm","glm","bayesglm")))
        tryCatch(AIC(m), error = function(e) NA_real_) else NA_real_
      
      BIC_val <- if (inherits(m, c("lm","glm","bayesglm")))
        tryCatch(BIC(m), error = function(e) NA_real_) else NA_real_
      
      ranking_df <- dplyr::add_row(
        ranking_df,
        Model = mname,
        R2    = R2,
        RMSE  = RMSE,
        RMSEP = RMSEP_val,
        AIC   = AIC_val,
        BIC   = BIC_val
      )
      
      global_summary <- dplyr::add_row(
        global_summary,
        Element = el,
        Model   = mname,
        R2      = R2,
        RMSE    = RMSE,
        RMSEP   = RMSEP_val,
        AIC     = AIC_val,
        BIC     = BIC_val
      )
      
      # Diagnostics
      sh_W <- sh_p <- NA_real_
      if (length(resid_vals) >= 3) {
        sh <- tryCatch(shapiro.test(resid_vals), error = function(e) NULL)
        if (!is.null(sh)) {
          sh_W <- unname(sh$statistic)
          sh_p <- sh$p.value
        }
      }
      
      bp_stat <- bp_p <- NA_real_
      if (inherits(m, c("lm","bayesglm"))) {
        bp <- tryCatch(bptest(m), error = function(e) NULL)
        if (!is.null(bp)) {
          bp_stat <- bp$statistic
          bp_p    <- bp$p.value
        }
      }
      
      diag_df <- dplyr::add_row(
        diag_df,
        Model             = mname,
        Shapiro_W         = sh_W,
        Shapiro_p         = sh_p,
        Shapiro_expl      = expl_shapiro,
        BreuschPagan_stat = bp_stat,
        BreuschPagan_p    = bp_p,
        BreuschPagan_expl = expl_bp,
        CV10_RMSE         = RMSEP_val,
        CV10_expl         = expl_cv10,
        AIC               = AIC_val,
        AIC_expl          = "Not used for ranking.",
        BIC               = BIC_val,
        BIC_expl          = "Not used for ranking."
      )
      
      # Predictions for XRF_new
      if (do_predictions) {
        pred_new_log <- tryCatch({
          if (inherits(m, "mvr"))
            as.numeric(predict(m, ncomp = attr(m, "best_ncomp"), newdata = XRF_new))
          else if (!is.null(m))
            as.numeric(predict(m, newdata = XRF_new))
          else
            rep(NA_real_, nrow(XRF_new))
        }, error = function(e) rep(NA_real_, nrow(XRF_new)))
        
        lower_log <- pred_new_log - 1.96 * RMSE
        upper_log <- pred_new_log + 1.96 * RMSE
        
        if (mname %in% c("OLS","WLS","OLS_weighted","WLS_weighted") && !is.null(m)) {
          pi_obj <- tryCatch(
            predict(m, newdata = XRF_new, interval = "prediction", level = 0.95),
            error = function(e) NULL
          )
          if (!is.null(pi_obj)) {
            lower_log <- pi_obj[, "lwr"]
            upper_log <- pi_obj[, "upr"]
          }
        }
        
        preds_df[[paste0(mname, "_", el, "_Pred_ppm")]] <- exp(pred_new_log)
        preds_df[[paste0(mname, "_", el, "_L95_ppm")]]  <- exp(lower_log)
        preds_df[[paste0(mname, "_", el, "_U95_ppm")]]  <- exp(upper_log)
      }
      
      # PLOTS: Pred vs Obs, Residuals, Influence, PLS RMSEP
      stats_label <- paste0(
        "R² = ", signif(R2, 3),
        "\nRMSEP = ", signif(RMSEP_val, 3),
        "\nRMSE = ", signif(RMSE, 3)
      )
      
      df_pred <- data.frame(
        Observed  = y,
        Predicted = pred,
        Site      = data$Site
      )
      
      # No transparency on points
      p_po <- ggplot(df_pred, aes(Observed, Predicted, colour = Site)) +
        geom_point(size = 1.25) +
        ggsci::scale_color_jco(name = "Site") +
        geom_abline(slope = 1, intercept = 0) +
        annotate(
          "text",
          x = min(df_pred$Observed,  na.rm = TRUE),
          y = max(df_pred$Predicted, na.rm = TRUE),
          hjust = 0, vjust = 1,
          label = stats_label, size = 2.7
        ) +
        theme_bw() + theme_small +
        ggtitle(paste(el, mname, "- Predicted vs Observed"))
      
      ggsave(
        file.path(element_dir, paste0(el, "_", mname, "_PredVsObs.pdf")),
        p_po, width = 13, height = 9, units = "cm"
      )
      
      res_df   <- data.frame(resid = resid_vals)
      hist_obj <- hist(resid_vals, plot = FALSE)
      
      p_res <- ggplot(res_df, aes(resid)) +
        geom_histogram(bins = 30, fill = "grey80") +
        geom_density(colour = "red") +
        annotate(
          "text",
          x = min(resid_vals, na.rm = TRUE),
          y = max(hist_obj$counts, na.rm = TRUE),
          hjust = 0, vjust = 1,
          label = stats_label, size = 2.7
        ) +
        theme_bw() + theme_small +
        ggtitle(paste(el, mname, "- Residuals"))
      
      ggsave(
        file.path(element_dir, paste0(el, "_", mname, "_Residuals.pdf")),
        p_res, width = 13, height = 9, units = "cm"
      )
      
      # Influence (lm only)
      if (inherits(m, "lm") && !inherits(m, "bayesglm")) {
        cooks <- cooks.distance(m)
        lev   <- hatvalues(m)
        rstd  <- rstandard(m)
        
        inf_df <- data.frame(
          obs      = seq_along(cooks),
          Cook     = cooks,
          Leverage = lev,
          StdResid = rstd
        )
        
        pdf(file.path(element_dir, paste0(el, "_", mname, "_Influence.pdf")),
            width = 13/2.54, height = 9/2.54)
        
        print(
          ggplot(inf_df, aes(obs, Cook)) +
            geom_col() +
            theme_bw() + theme_small +
            ggtitle(paste(el, mname, "- Cook's Distance"))
        )
        
        print(
          ggplot(inf_df, aes(Leverage, StdResid)) +
            geom_point(size = 1.25) +
            theme_bw() + theme_small +
            ggtitle(paste(el, mname, "- Leverage"))
        )
        
        dev.off()
      }
      
      # PLS RMSEP
      if (inherits(m, "mvr")) {
        pdf(file.path(element_dir, paste0(el, "_", mname, "_PLS_RMSEP.pdf")),
            width = 13/2.54, height = 9/2.54)
        
        validationplot(
          m,
          val.type = "RMSEP",
          main = paste(el, mname, "- RMSEP vs Components")
        )
        
        dev.off()
      }
      
    } # end model loop
    
    # ===============================================================
    # STORE MODELS + RMSE FOR LATER ACE SECTIONS
    # ===============================================================
    model_store[[el]] <- model_list
    
    calibration_rmse[[el]] <- ranking_df %>%
      dplyr::select(Model, RMSE) %>%
      tibble::deframe()
    
    # preds_store filled in Section 8.5
    
    # ===============================================================
    # 8.4 MODEL RANKING & DIAG OUTPUT
    # ===============================================================
    if (nrow(ranking_df) > 0) {
      ranking_df <- ranking_df %>%
        dplyr::arrange(dplyr::desc(R2), RMSEP, RMSE) %>%
        dplyr::mutate(Rank = dplyr::row_number())
      
      ranking_expl_df <- ranking_df %>%
        dplyr::mutate(
          Rank_explanation  = expl_rank,
          R2_explanation    = expl_r2,
          RMSE_explanation  = expl_rmse,
          RMSEP_explanation = expl_rmsep
        )
      
      write.csv(
        ranking_expl_df,
        file.path(element_dir, paste0(el, "_Model_Ranking.csv")),
        row.names = FALSE
      )
    }
    
    write.csv(
      diag_df,
      file.path(element_dir, paste0(el, "_Diagnostics_Summary.csv")),
      row.names = FALSE
    )
    
    # ===============================================================
    # 8.5 SAVE PREDICTIONS
    # ===============================================================
    if (!is.null(preds_df)) {
      write.csv(
        preds_df,
        file.path(element_dir, paste0(el, "_preds.csv")),
        row.names = FALSE
      )
      preds_store[[el]] <- preds_df
    }
    # ===============================================================
    # 8.6 README
    # ===============================================================
    readme_txt <- c(
      paste("Element:", el),
      "",
      "Calibration responses are ln(mg kg^-1).",
      "Predictions stored as *_Pred_ppm, *_L95_ppm, *_U95_ppm.",
      "",
      "Models: OLS, WLS, OLS(wt), WLS(wt), Bayes, RF, PLS(LOO), PLS(k).",
      "Ranking: highest R², then lowest RMSEP, then lowest RMSE."
    )
    
    writeLines(readme_txt, file.path(element_dir, "README.txt"))
    
  } # end element loop
  
  # ===============================================================
  # 9. SITE-LEVEL A4 PLOTS (DEPTH + AGE) WITH ICP_ppm ERRORBARS
  # ===============================================================
  
  dynamic_layout <- function(n_panels, max_cols = 3) {
    ncol <- min(max_cols, max(1, ceiling(sqrt(n_panels))))
    nrow <- ceiling(n_panels / ncol)
    list(nrow = nrow, ncol = ncol)
  }
  
  if (!is.null(XRF_new) && length(preds_store) > 0) {
    
    has_site    <- "Site"     %in% names(XRF_new)
    has_depth   <- "depth"    %in% names(XRF_new)
    has_age_new <- "SH20_age" %in% names(XRF_new)
    
    sites_dir_all  <- file.path(all_base,  "Sites")
    sites_dir_best <- file.path(best_base, "Sites")
    dir.create(sites_dir_all,  recursive = TRUE, showWarnings = FALSE)
    dir.create(sites_dir_best, recursive = TRUE, showWarnings = FALSE)
    
    if (has_site) {
      site_levels <- sort(unique(XRF_new$Site))
      multivar_labels <- c("Bayes","RF","PLS(LOO)","PLS(k)")
      
      for (s in site_levels) {
        
        # ------------------------------------------
        # 9.1 DEPTH PROFILES PER SITE
        # ------------------------------------------
        if (has_depth) {
          
          depth_plots_all  <- list()
          depth_plots_best <- list()
          
          for (el in elements) {
            
            df_el <- preds_store[[el]]
            if (is.null(df_el) ||
                !all(c("Site","depth") %in% names(df_el))) next
            
            df_site <- df_el %>% dplyr::filter(Site == s)
            if (nrow(df_site) == 0) next
            
            pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                              names(df_site), value = TRUE)
            if (length(pred_cols) == 0) next
            
            df_long <- df_site %>%
              dplyr::select(depth, Site, dplyr::all_of(pred_cols)) %>%
              tidyr::pivot_longer(
                cols = dplyr::all_of(pred_cols),
                names_to  = "Model",
                values_to = "Pred_ppm"
              )
            
            df_long$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long$Model)
            df_long$Model <- factor(df_long$Model,
                                    levels = model_levels,
                                    labels = model_labels)
            
            df_long_best <- df_long %>%
              dplyr::filter(Model %in% multivar_labels)
            
            # ICP overlay from ICP_ppm (depth domain)
            obs_col <- paste0(el, "_ICP")
            sd_col  <- paste0(el, "_ICP_sd")
            has_sd_el <- sd_col %in% names(ICP_ppm)
            
            ace_overlay <- NULL
            if (all(c("Site","depth",obs_col) %in% names(ICP_ppm))) {
              ace_overlay <- ICP_ppm %>%
                dplyr::filter(Site == s,
                              !is.na(.data[[obs_col]]),
                              !is.na(depth)) %>%
                dplyr::mutate(
                  ICP_ppm    = .data[[obs_col]],
                  ICP_sd_ppm = if (has_sd_el) .data[[sd_col]] else NA_real_
                )
              if (nrow(ace_overlay) == 0) ace_overlay <- NULL
            }
            
            # No transparency lines
            p_depth_all <- ggplot(df_long,
                                  aes(x = Pred_ppm, y = depth,
                                      colour = Model, group = Model)) +
              geom_path(linewidth = 0.4) +
              scale_y_reverse() +
              labs(
                x = bquote(.(el)~"(mg"~kg^{-1}*")"),
                y = "Depth (cm)",
                title = paste("Site", s, "-", el, "(depth, all models)")
              ) +
              scale_color_manual(values = cols_all, breaks = names(cols_all)) +
              theme_bw() + theme_small
            
            p_depth_best <- ggplot(df_long_best,
                                   aes(x = Pred_ppm, y = depth,
                                       colour = Model, group = Model)) +
              geom_path(linewidth = 0.4) +
              scale_y_reverse() +
              labs(
                x = bquote(.(el)~"(mg"~kg^{-1}*")"),
                y = "Depth (cm)",
                title = paste("Site", s, "-", el, "(depth, multivariate)")
              ) +
              scale_color_manual(
                values = cols_all,
                breaks = c(multivar_labels, "ICP-MS (Observed)")
              ) +
              theme_bw() + theme_small
            
            if (!is.null(ace_overlay)) {
              # ICP error bars (horizontal)
              if (has_sd_el) {
                p_depth_all <- p_depth_all +
                  geom_segment(
                    data = ace_overlay,
                    aes(
                      y = depth,
                      x = ICP_ppm - ICP_sd_ppm,
                      xend = ICP_ppm + ICP_sd_ppm
                    ),
                    inherit.aes = FALSE,
                    colour = "pink"
                  )
                
                p_depth_best <- p_depth_best +
                  geom_segment(
                    data = ace_overlay,
                    aes(
                      y = depth,
                      x = ICP_ppm - ICP_sd_ppm,
                      xend = ICP_ppm + ICP_sd_ppm
                    ),
                    inherit.aes = FALSE,
                    colour = "pink"
                  )
              }
              
              # ICP path + points
              p_depth_all <- p_depth_all +
                geom_path(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth,
                      colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth,
                      colour = "ICP-MS (Observed)"),
                  shape = 21,
                  fill  = "white",
                  stroke = 0.4,
                  size   = 1.4,
                  inherit.aes = FALSE
                )
              
              p_depth_best <- p_depth_best +
                geom_path(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth,
                      colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth,
                      colour = "ICP-MS (Observed)"),
                  shape = 21,
                  fill  = "white",
                  stroke = 0.4,
                  size   = 1.4,
                  inherit.aes = FALSE
                )
            }
            
            depth_plots_all[[el]]  <- p_depth_all
            depth_plots_best[[el]] <- p_depth_best
          }
          
          if (length(depth_plots_all) > 0) {
            panels <- lapply(elements, function(e) {
              if (is.null(depth_plots_all[[e]])) ggplot() + theme_void() else depth_plots_all[[e]]
            })
            used <- sum(!vapply(panels, function(p) inherits(p$theme, "theme_void"), logical(1)))
            if (used > 0) {
              lay <- dynamic_layout(used)
              combined <- wrap_plots(
                panels,
                ncol   = lay$ncol,
                nrow   = lay$nrow,
                guides = "collect"
              ) & theme(legend.position = "bottom")
              
              ggsave(
                file.path(sites_dir_all,
                          paste0("Site_", s, "_Profiles_depth_all.pdf")),
                combined, width = 21, height = 29.7, units = "cm"
              )
            }
          }
          
          if (length(depth_plots_best) > 0) {
            panels <- lapply(elements, function(e) {
              if (is.null(depth_plots_best[[e]])) ggplot() + theme_void() else depth_plots_best[[e]]
            })
            used <- sum(!vapply(panels, function(p) inherits(p$theme, "theme_void"), logical(1)))
            if (used > 0) {
              lay <- dynamic_layout(used)
              combined <- wrap_plots(
                panels,
                ncol   = lay$ncol,
                nrow   = lay$nrow,
                guides = "collect"
              ) & theme(legend.position = "bottom")
              
              ggsave(
                file.path(sites_dir_best,
                          paste0("Site_", s, "_Profiles_depth_multivariate.pdf")),
                combined, width = 21, height = 29.7, units = "cm"
              )
            }
          }
        } # end depth
        
        # ------------------------------------------
        # 9.2 AGE PROFILES PER SITE (ppm vs age)
        # ------------------------------------------
        if (has_age_new) {
          
          age_plots_all  <- list()
          age_plots_best <- list()
          
          for (el in elements) {
            
            df_el <- preds_store[[el]]
            if (is.null(df_el) ||
                !all(c("Site","SH20_age") %in% names(df_el))) next
            
            df_site <- df_el %>% dplyr::filter(Site == s)
            if (nrow(df_site) == 0) next
            
            pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                              names(df_site), value = TRUE)
            if (length(pred_cols) == 0) next
            
            df_long <- df_site %>%
              dplyr::select(SH20_age, Site, dplyr::all_of(pred_cols)) %>%
              tidyr::pivot_longer(
                cols = dplyr::all_of(pred_cols),
                names_to  = "Model",
                values_to = "Pred_ppm"
              )
            
            df_long$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long$Model)
            df_long$Model <- factor(df_long$Model,
                                    levels = model_levels,
                                    labels = model_labels)
            
            df_long_best <- df_long %>%
              dplyr::filter(Model %in% multivar_labels)
            
            # ICP overlay from ICP_ppm (age domain)
            obs_col <- paste0(el, "_ICP")
            sd_col  <- paste0(el, "_ICP_sd")
            has_sd_el <- sd_col %in% names(ICP_ppm)
            
            ace_overlay_age <- NULL
            if (all(c("Site","SH20_age",obs_col) %in% names(ICP_ppm))) {
              ace_overlay_age <- ICP_ppm %>%
                dplyr::filter(Site == s,
                              !is.na(.data[[obs_col]]),
                              !is.na(SH20_age)) %>%
                dplyr::mutate(
                  ICP_ppm    = .data[[obs_col]],
                  ICP_sd_ppm = if (has_sd_el) .data[[sd_col]] else NA_real_
                )
              if (nrow(ace_overlay_age) == 0) ace_overlay_age <- NULL
            }
            
            p_age_all <- ggplot(df_long,
                                aes(x = Pred_ppm, y = SH20_age,
                                    colour = Model, group = Model)) +
              geom_path(linewidth = 0.4) +
              labs(
                x = bquote(.(el)~"(mg"~kg^{-1}*")"),
                y = "Age (cal a BP)",
                title = paste("Site", s, "-", el, "(age, all models)")
              ) +
              scale_color_manual(values = cols_all, breaks = names(cols_all)) +
              theme_bw() + theme_small
            
            p_age_best <- ggplot(df_long_best,
                                 aes(x = Pred_ppm, y = SH20_age,
                                     colour = Model, group = Model)) +
              geom_path(linewidth = 0.4) +
              labs(
                x = bquote(.(el)~"(mg"~kg^{-1}*")"),
                y = "Age (cal a BP)",
                title = paste("Site", s, "-", el, "(age, multivariate)")
              ) +
              scale_color_manual(
                values = cols_all,
                breaks = c(multivar_labels, "ICP-MS (Observed)")
              ) +
              theme_bw() + theme_small
            
            if (!is.null(ace_overlay_age)) {
              # ICP error bars (horizontal in ppm)
              if (has_sd_el) {
                p_age_all <- p_age_all +
                  geom_segment(
                    data = ace_overlay_age,
                    aes(y = SH20_age,
                        x = ICP_ppm - ICP_sd_ppm,
                        xend = ICP_ppm + ICP_sd_ppm),
                    inherit.aes = FALSE,
                    colour = "pink"
                  )
                
                p_age_best <- p_age_best +
                  geom_segment(
                    data = ace_overlay_age,
                    aes(y = SH20_age,
                        x = ICP_ppm - ICP_sd_ppm,
                        xend = ICP_ppm + ICP_sd_ppm),
                    inherit.aes = FALSE,
                    colour = "pink"
                  )
              }
              
              # ICP path + points
              p_age_all <- p_age_all +
                geom_path(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age,
                      colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age,
                      colour = "ICP-MS (Observed)"),
                  shape = 21,
                  fill  = "white",
                  stroke = 0.4,
                  size   = 1.4,
                  inherit.aes = FALSE
                )
              
              p_age_best <- p_age_best +
                geom_path(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age,
                      colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age,
                      colour = "ICP-MS (Observed)"),
                  shape = 21,
                  fill  = "white",
                  stroke = 0.4,
                  size   = 1.4,
                  inherit.aes = FALSE
                )
            }
            
            age_plots_all[[el]]  <- p_age_all
            age_plots_best[[el]] <- p_age_best
          }
          
          if (length(age_plots_all) > 0) {
            panels <- lapply(elements, function(e) {
              if (is.null(age_plots_all[[e]])) ggplot() + theme_void() else age_plots_all[[e]]
            })
            used <- sum(!vapply(panels, function(p) inherits(p$theme, "theme_void"), logical(1)))
            if (used > 0) {
              lay <- dynamic_layout(used)
              combined <- wrap_plots(
                panels,
                ncol   = lay$ncol,
                nrow   = lay$nrow,
                guides = "collect"
              ) & theme(legend.position = "bottom")
              
              ggsave(
                file.path(sites_dir_all,
                          paste0("Site_", s, "_Profiles_age_all.pdf")),
                combined, width = 21, height = 29.7, units = "cm"
              )
            }
          }
          
          if (length(age_plots_best) > 0) {
            panels <- lapply(elements, function(e) {
              if (is.null(age_plots_best[[e]])) ggplot() + theme_void() else age_plots_best[[e]]
            })
            used <- sum(!vapply(panels, function(p) inherits(p$theme, "theme_void"), logical(1)))
            if (used > 0) {
              lay <- dynamic_layout(used)
              combined <- wrap_plots(
                panels,
                ncol   = lay$ncol,
                nrow   = lay$nrow,
                guides = "collect"
              ) & theme(legend.position = "bottom")
              
              ggsave(
                file.path(sites_dir_best,
                          paste0("Site_", s, "_Profiles_age_multivariate.pdf")),
                combined, width = 21, height = 29.7, units = "cm"
              )
            }
          }
        } # end age
        
      } # end site loop
    }
  }
  
  # ===============================================================
  # 10. ELEMENT-LEVEL MULTI-SITE AGE–CONC COMPARISON (WITH ICP_ppm)
  # ===============================================================
  
  if (!is.null(XRF_new) &&
      length(preds_store) > 0 &&
      "Site"     %in% names(XRF_new) &&
      "SH20_age" %in% names(XRF_new) &&
      "Site"     %in% names(ICP_ppm) &&
      "SH20_age" %in% names(ICP_ppm)) {
    
    elements_dir_all  <- file.path(all_base,  "elements")
    elements_dir_best <- file.path(best_base, "elements")
    dir.create(elements_dir_all,  recursive = TRUE, showWarnings = FALSE)
    dir.create(elements_dir_best, recursive = TRUE, showWarnings = FALSE)
    
    # Global age range
    age_all <- c(XRF_new$SH20_age, ICP_ppm$SH20_age)
    age_all <- age_all[!is.na(age_all)]
    age_min <- if (length(age_all)) min(age_all) else NA_real_
    age_max <- if (length(age_all)) max(age_all) else NA_real_
    
    site_levels_all <- sort(unique(c(as.character(XRF_new$Site),
                                     as.character(ICP_ppm$Site))))
    
    multivar_labels <- c("Bayes","RF","PLS(LOO)","PLS(k)")
    
    for (el in elements) {
      
      df_el <- preds_store[[el]]
      if (is.null(df_el)) next
      
      pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                        names(df_el), value = TRUE)
      if (length(pred_cols) == 0) next
      
      df_long_all <- df_el %>%
        dplyr::filter(!is.na(SH20_age)) %>%
        dplyr::select(Site, SH20_age, dplyr::all_of(pred_cols)) %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(pred_cols),
          names_to  = "Model",
          values_to = "Pred_ppm"
        )
      
      df_long_all$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long_all$Model)
      df_long_all$Model <- factor(df_long_all$Model,
                                  levels = model_levels,
                                  labels = model_labels)
      
      df_long_best <- df_long_all %>% dplyr::filter(Model %in% multivar_labels)
      
      obs_col <- paste0(el, "_ICP")
      sd_col  <- paste0(el, "_ICP_sd")
      has_sd_el <- sd_col %in% names(ICP_ppm)
      
      ace_all <- ICP_ppm %>%
        dplyr::filter(!is.na(SH20_age),
                      !is.na(.data[[obs_col]])) %>%
        dplyr::mutate(
          ICP_ppm    = .data[[obs_col]],
          ICP_sd_ppm = if (has_sd_el) .data[[sd_col]] else NA_real_
        )
      
      # Global y-range for log-equal plots
      y_pred_vals <- df_long_all$Pred_ppm
      y_obs_vals  <- ace_all$ICP_ppm
      y_all_vals  <- c(y_pred_vals, y_obs_vals)
      y_all_vals  <- y_all_vals[is.finite(y_all_vals)]
      
      if (!length(y_all_vals) || min(y_all_vals) <= 0) {
        y_min_el <- NA_real_
        y_max_el <- NA_real_
      } else {
        y_min_el <- min(y_all_vals)
        y_max_el <- max(y_all_vals)
      }
      
      if (!is.na(y_min_el) && !is.na(y_max_el)) {
        y_log_breaks <- 10^(floor(log10(y_min_el)) : ceiling(log10(y_max_el)))
        y_log_breaks <- y_log_breaks[y_log_breaks >= y_min_el & y_log_breaks <= y_max_el]
        if (!length(y_log_breaks)) {
          y_log_breaks <- scales::log_breaks()(c(y_min_el, y_max_el))
        }
        y_log_labels <- scales::scientific_format()(y_log_breaks)
      } else {
        y_log_breaks <- NULL
        y_log_labels <- NULL
      }
      
      pls_pred_col <- paste0("PLS_kfold_", el, "_Pred_ppm")
      pls_l95_col  <- paste0("PLS_kfold_", el, "_L95_ppm")
      pls_u95_col  <- paste0("PLS_kfold_", el, "_U95_ppm")
      has_pls_cols <- all(c(pls_pred_col, pls_l95_col, pls_u95_col) %in% names(df_el))
      
      plot_list_all           <- list()
      plot_list_best          <- list()
      plot_list_all_y_equal   <- list()
      plot_list_best_y_equal  <- list()
      
      for (s in site_levels_all) {
        
        df_pred_site_all  <- df_long_all  %>% dplyr::filter(Site == s)
        df_pred_site_best <- df_long_best %>% dplyr::filter(Site == s)
        
        ace_site <- ace_all %>% dplyr::filter(Site == s)
        
        df_raw_xrf_site <- df_el   %>% dplyr::filter(Site == s)
        icp_site_raw    <- ICP_ppm %>% dplyr::filter(Site == s)
        n_xrf <- sum(!is.na(df_raw_xrf_site$SH20_age))
        n_icp <- sum(!is.na(icp_site_raw[[obs_col]]) & !is.na(icp_site_raw$SH20_age))
        n_label <- paste0("n = ", n_icp, " (ICP), ", n_xrf, " (XRF)")
        
        if (nrow(ace_site) > 0) {
          icp_lower_site <- ace_site$ICP_ppm - ace_site$ICP_sd_ppm
          icp_upper_site <- ace_site$ICP_ppm + ace_site$ICP_sd_ppm
          icp_lower_min  <- min(icp_lower_site, na.rm = TRUE)
          icp_upper_max  <- max(icp_upper_site, na.rm = TRUE)
          site_median    <- median(ace_site$ICP_ppm, na.rm = TRUE)
        } else {
          icp_lower_min <- icp_upper_max <- site_median <- NA_real_
        }
        
        if (has_pls_cols) {
          df_pls_site <- df_el %>%
            dplyr::filter(Site == s) %>%
            dplyr::select(
              Site,
              SH20_age,
              Pred_ppm = !!sym(pls_pred_col),
              L95_ppm  = !!sym(pls_l95_col),
              U95_ppm  = !!sym(pls_u95_col)
            )
        } else {
          df_pls_site <- NULL
        }
        
        # Exceedance test: upper CI > ICP + 2*sd
        ast_df <- NULL
        if (!is.null(df_pls_site) &&
            nrow(df_pls_site) > 0 &&
            !is.na(icp_upper_max)) {
          
          thresh <- icp_upper_max + 2 * max(ace_site$ICP_sd_ppm, na.rm = TRUE)
          
          exceed_idx <- with(
            df_pls_site,
            !is.na(U95_ppm) & U95_ppm > thresh
          )
          
          if (any(exceed_idx)) ast_df <- df_pls_site[exceed_idx, , drop = FALSE]
        }
        
        y_top_free <- suppressWarnings(
          max(c(df_pred_site_all$Pred_ppm, ace_site$ICP_ppm), na.rm = TRUE)
        )
        if (!is.finite(y_top_free)) y_top_free <- NA_real_
        
        # ===========================================================
        # ALL MODELS – NO RIBBON, NO ERROR BARS
        # ===========================================================
        if (nrow(df_pred_site_all) == 0 && nrow(ace_site) == 0) {
          p_all <- ggplot() + theme_void()
        } else {
          p_all <- ggplot(df_pred_site_all,
                          aes(x = SH20_age, y = Pred_ppm, colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(all models)")
            ) +
            scale_color_manual(values = cols_all, breaks = names(cols_all)) +
            theme_bw() + theme_small +
            coord_cartesian(xlim = c(age_min, age_max))
          
          if (nrow(ace_site) > 0) {
            p_all <- p_all +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white", stroke = 0.4, size = 1.4,
                inherit.aes = FALSE
              )
          }
          
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_all <- p_all +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red", size = 3, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          p_all <- p_all +
            annotate("text", x = age_max, y = y_top_free,
                     label = n_label, hjust = 1, vjust = 1, size = 2.7)
        }
        
        # ===========================================================
        # ALL MODELS – Y-EQUAL (LOG10), WITH GREY RIBBON, NO ERROR BARS
        # ===========================================================
        if (inherits(p_all$theme, "theme_void") ||
            is.na(y_min_el) || is.na(y_max_el)) {
          p_all_y <- p_all
        } else {
          p_all_y <- ggplot(df_pred_site_all,
                            aes(x = SH20_age, y = Pred_ppm, colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(all models, y-equal)")
            ) +
            scale_color_manual(values = cols_all, breaks = names(cols_all)) +
            theme_bw() + theme_small +
            coord_cartesian(xlim = c(age_min, age_max))
          
          if (!is.null(df_pls_site)) {
            p_all_y <- p_all_y +
              geom_ribbon(
                data = df_pls_site,
                aes(x = SH20_age, ymin = L95_ppm, ymax = U95_ppm),
                inherit.aes = FALSE,
                fill = "grey80", alpha = 0.5
              )
          }
          
          if (nrow(ace_site) > 0) {
            p_all_y <- p_all_y +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white", stroke = 0.4, size = 1.4,
                inherit.aes = FALSE
              )
          }
          
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_all_y <- p_all_y +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red", size = 3, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          p_all_y <- p_all_y +
            annotate("text", x = age_max, y = y_top_free,
                     label = n_label, hjust = 1, vjust = 1, size = 2.7) +
            scale_y_log10(
              limits = c(y_min_el, y_max_el),
              breaks = y_log_breaks,
              labels = y_log_labels,
              oob = scales::oob_squish
            )
        }
        
        # ===========================================================
        # MULTIVARIATE — WITH GREY RIBBON, NO ERROR BARS
        # ===========================================================
        if (nrow(df_pred_site_best) == 0 && nrow(ace_site) == 0) {
          p_best <- ggplot() + theme_void()
        } else {
          p_best <- ggplot(df_pred_site_best,
                           aes(x = SH20_age, y = Pred_ppm, colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(multivariate)")
            ) +
            scale_color_manual(values = cols_all,
                               breaks = c(multivar_labels, "ICP-MS (Observed)")) +
            theme_bw() + theme_small +
            coord_cartesian(xlim = c(age_min, age_max))
          
          if (!is.null(df_pls_site)) {
            p_best <- p_best +
              geom_ribbon(
                data = df_pls_site,
                aes(x = SH20_age, ymin = L95_ppm, ymax = U95_ppm),
                inherit.aes = FALSE,
                fill = "grey80", alpha = 0.5
              )
          }
          
          if (nrow(ace_site) > 0) {
            p_best <- p_best +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white", stroke = 0.4, size = 1.4,
                inherit.aes = FALSE
              )
          }
          
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_best <- p_best +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red", size = 3, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          p_best <- p_best +
            annotate("text", x = age_max, y = y_top_free,
                     label = n_label, hjust = 1, vjust = 1, size = 2.7)
        }
        
        # ===========================================================
        # MULTIVARIATE Y-EQUAL — NO ERROR BARS, WITH RIBBON
        # ===========================================================
        if (inherits(p_best$theme, "theme_void") ||
            is.na(y_min_el) || is.na(y_max_el)) {
          p_best_y <- p_best
        } else {
          p_best_y <- ggplot(df_pred_site_best,
                             aes(x = SH20_age, y = Pred_ppm, colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(multivariate, y-equal)")
            ) +
            scale_color_manual(values = cols_all,
                               breaks = c(multivar_labels, "ICP-MS (Observed)")) +
            theme_bw() + theme_small +
            coord_cartesian(xlim = c(age_min, age_max))
          
          if (!is.null(df_pls_site)) {
            p_best_y <- p_best_y +
              geom_ribbon(
                data = df_pls_site,
                aes(x = SH20_age, ymin = L95_ppm, ymax = U95_ppm),
                inherit.aes = FALSE,
                fill = "grey80", alpha = 0.5
              )
          }
          
          if (nrow(ace_site) > 0) {
            p_best_y <- p_best_y +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm, colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white", stroke = 0.4, size = 1.4,
                inherit.aes = FALSE
              )
          }
          
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_best_y <- p_best_y +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red", size = 3, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          p_best_y <- p_best_y +
            annotate("text", x = age_max, y = y_top_free,
                     label = n_label, hjust = 1, vjust = 1, size = 2.7) +
            scale_y_log10(
              limits = c(y_min_el, y_max_el),
              breaks = y_log_breaks,
              labels = y_log_labels,
              oob = scales::oob_squish
            )
        }
        
        # store
        plot_list_all[[s]]          <- p_all
        plot_list_all_y_equal[[s]]  <- p_all_y
        plot_list_best[[s]]         <- p_best
        plot_list_best_y_equal[[s]] <- p_best_y
      }
      
      keep_nonempty <- function(lst) {
        lst[vapply(lst, function(p)
          !inherits(p$theme, "theme_void"), logical(1))]
      }
      
      combine_1col <- function(plot_list) {
        keep <- keep_nonempty(plot_list)
        if (length(keep) == 0) return(NULL)
        n <- length(keep)
        keep <- Map(
          function(p, i) if (i < n) p + labs(x = NULL) else p,
          keep, seq_along(keep)
        )
        wrap_plots(keep, ncol = 1, nrow = n, guides = "collect") &
          theme(legend.position = "bottom")
      }
      
      combined_all <- combine_1col(plot_list_all)
      if (!is.null(combined_all)) {
        ggsave(
          file.path(elements_dir_all, paste0(el, "_model_age_comparison_all.pdf")),
          combined_all, width = 21, height = 29.7, units = "cm"
        )
      }
      
      combined_all_y <- combine_1col(plot_list_all_y_equal)
      if (!is.null(combined_all_y)) {
        ggsave(
          file.path(elements_dir_all, paste0(el, "_model_age_comparison_all_y_equal.pdf")),
          combined_all_y, width = 21, height = 29.7, units = "cm"
        )
      }
      
      combined_best <- combine_1col(plot_list_best)
      if (!is.null(combined_best)) {
        ggsave(
          file.path(elements_dir_best, paste0(el, "_model_age_comparison_multivariate.pdf")),
          combined_best, width = 21, height = 29.7, units = "cm"
        )
      }
      
      combined_best_y <- combine_1col(plot_list_best_y_equal)
      if (!is.null(combined_best_y)) {
        ggsave(
          file.path(elements_dir_best, paste0(el, "_model_age_comparison_multivariate_y_equal.pdf")),
          combined_best_y, width = 21, height = 29.7, units = "cm"
        )
      }
    }
  }
  
  # ===============================================================
  # 11. GLOBAL SUMMARY ACROSS ELEMENTS
  # ===============================================================
  
  if (nrow(global_summary) > 0) {
    
    global_ranked <- global_summary %>%
      dplyr::arrange(
        Element,
        dplyr::desc(R2),
        RMSEP,
        RMSE
      ) %>%
      dplyr::group_by(Element) %>%
      dplyr::mutate(Rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    write.csv(
      global_ranked,
      file.path(out_dir, summary_all_csv),
      row.names = FALSE
    )
    
    best_models <- global_ranked %>%
      dplyr::filter(Rank == 1)
    
    write.csv(
      best_models,
      file.path(out_dir, best_models_csv),
      row.names = FALSE
    )
  }
  
  # ===============================================================
  # 11.1 PPM-SCALE SUMMARY
  # ===============================================================
  
  if (exists("global_ranked") && nrow(global_ranked) > 0) {
    
    elem_means_list <- lapply(elements, function(el) {
      y_col <- paste0(el, "_ICP")
      if (!y_col %in% names(data)) {
        return(data.frame(
          Element      = el,
          Mean_log_ppm = NA_real_,
          Mean_ppm     = NA_real_
        ))
      }
      vals <- data[[y_col]]
      vals <- vals[!is.na(vals)]
      if (!length(vals)) {
        return(data.frame(
          Element      = el,
          Mean_log_ppm = NA_real_,
          Mean_ppm     = NA_real_
        ))
      }
      mean_log <- mean(vals)
      mean_ppm <- exp(mean_log)
      data.frame(
        Element      = el,
        Mean_log_ppm = mean_log,
        Mean_ppm     = mean_ppm
      )
    })
    
    elem_means_df <- do.call(rbind, elem_means_list)
    
    global_ranked_ppm <- global_ranked %>%
      dplyr::left_join(elem_means_df, by = "Element") %>%
      dplyr::mutate(
        RMSE_log      = RMSE,
        RMSEP_log     = RMSEP,
        RMSE_factor   = exp(RMSE_log),
        RMSEP_factor  = exp(RMSEP_log),
        RMSE_abs_ppm  = RMSE_log  * Mean_ppm,
        RMSEP_abs_ppm = RMSEP_log * Mean_ppm
      ) %>%
      dplyr::select(
        Element, Model, Rank, R2,
        RMSE_log, RMSE_factor, RMSE_abs_ppm,
        RMSEP_log, RMSEP_factor, RMSEP_abs_ppm,
        Mean_log_ppm, Mean_ppm,
        AIC, BIC
      )
    
    write.csv(
      global_ranked_ppm,
      file.path(out_dir, "AllElements_ModelSummary_ppm.csv"),
      row.names = FALSE
    )
  }
  
  # ===============================================================
  # 12. FINISH
  # ===============================================================
  
  message("=== COMPLETED SUCCESSFULLY ===")
  message("All results saved in: ", out_dir)
  
  # ===============================================================
  # USER OPTION FOR ACE PREDICTIONS (TRUE/FALSE)
  # ===============================================================
  
  # Default behaviour: run ACE predictions if ACE_dataset is provided
  do_ACE_predictions <- TRUE
  
  # If user supplies ACE_predict argument, override default
  if (!missing(ACE_predict)) {
    if (is.logical(ACE_predict)) {
      do_ACE_predictions <- ACE_predict
    } else {
      warning("ACE_predict must be TRUE/FALSE – using default TRUE.")
    }
  }
  
  if (do_ACE_predictions) {
    message("ACE predictions ENABLED.")
  } else {
    message("ACE predictions DISABLED.")
  }

  # ===============================================================
  # 13. ACE_dataset PREDICTION (OPTIONAL)
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== Running ACE prediction across all elements ===")
    
    # Ensure ACE has required predictors
    if (!"Site" %in% names(ACE_dataset))
      stop("ACE_dataset must contain a 'Site' column.")
    
    if (!"SH20_age" %in% names(ACE_dataset))
      stop("ACE_dataset must contain 'SH20_age'.")
    
    # Create output root
    ace_root <- file.path(all_base, "Calibration_pred")
    dir.create(ace_root, recursive = TRUE, showWarnings = FALSE)
    
    # Loop over elements
    for (el in elements) {
      
      el_dir <- file.path(ace_root, el)
      dir.create(el_dir, recursive = TRUE, showWarnings = FALSE)
      
      message("ACE predictions for element: ", el)
      
      # Prepare output df
      ace_df <- ACE_dataset %>%
        dplyr::select(Site, depth, SH20_age, dplyr::all_of(el)) %>%
        dplyr::rename(XRF_value = !!sym(el))
      
      # Access fitted models + RMSE
      models_el <- model_store[[el]]
      rmse_el   <- calibration_rmse[[el]]
      
      # Compute predictions, store per model
      for (mname in names(models_el)) {
        
        m <- models_el[[mname]]
        
        # Predictions in log space
        pred_log <- tryCatch({
          if (inherits(m, "mvr"))
            as.numeric(predict(m, ncomp = attr(m, "best_ncomp"), newdata = ACE_dataset))
          else
            as.numeric(predict(m, newdata = ACE_dataset))
        }, error = function(e) rep(NA_real_, nrow(ACE_dataset)))
        
        # Compute CI
        rmse_val <- rmse_el[[mname]]
        if (is.null(rmse_val) || is.na(rmse_val)) rmse_val <- NA_real_
        
        L95_log <- pred_log - 1.96 * rmse_val
        U95_log <- pred_log + 1.96 * rmse_val
        
        # Save ppm-scale values
        ace_df[[paste0(mname, "_", el, "_Pred_ppm")]] <- exp(pred_log)
        ace_df[[paste0(mname, "_", el, "_L95_ppm")]]  <- exp(L95_log)
        ace_df[[paste0(mname, "_", el, "_U95_ppm")]]  <- exp(U95_log)
      }
      
      # Save predictions
      write.csv(
        ace_df,
        file.path(el_dir, paste0(el, "_ACE_predictions.csv")),
        row.names = FALSE
      )
    }
    
    message("=== Finished ACE predictions ===")
  }
  
  # ===============================================================
  # 14. ACE DIAGNOSTIC EVALUATION vs ICP_ppm
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== Running ACE diagnostic evaluation ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    
    for (el in elements) {
      
      el_dir <- file.path(ace_root, el)
      pred_file <- file.path(el_dir, paste0(el, "_ACE_predictions.csv"))
      if (!file.exists(pred_file)) next
      
      ace_df <- read.csv(pred_file)
      
      obs_col <- paste0(el, "_ICP")
      sd_col  <- paste0(el, "_ICP_sd")
      
      if (!obs_col %in% names(ICP_ppm)) next
      
      # join ACE predictions with ICP measured data
      joined <- ace_df %>%
        left_join(ICP_ppm %>% dplyr::select(Site, depth, SH20_age, !!sym(obs_col), !!sym(sd_col)),
                  by = c("Site", "depth", "SH20_age"))
      
      # Create diagnostics and model comparison just like Section 8
      diag_out <- data.frame()
      
      for (mname in names(model_store[[el]])) {
        
        pred_nm <- paste0(mname, "_", el, "_Pred_ppm")
        if (!pred_nm %in% names(joined)) next
        
        pred <- joined[[pred_nm]]
        obs  <- joined[[obs_col]]
        
        ok <- !is.na(pred) & !is.na(obs)
        if (!any(ok)) next
        
        SS_res <- sum((obs[ok] - pred[ok])^2)
        SS_tot <- sum((obs[ok] - mean(obs[ok]))^2)
        
        R2   <- 1 - SS_res / SS_tot
        RMSE <- sqrt(mean((obs[ok] - pred[ok])^2))
        
        diag_out <- rbind(
          diag_out,
          data.frame(
            Model = mname,
            R2    = R2,
            RMSE  = RMSE
          )
        )
      }
      
      write.csv(
        diag_out,
        file.path(el_dir, paste0(el, "_ACE_diagnostics.csv")),
        row.names = FALSE
      )
    }
  }
  
  # ===============================================================
  # 15. ACE PREDICTION DIAGNOSTIC PLOTS
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== Creating ACE diagnostic plots ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    
    for (el in elements) {
      
      el_dir <- file.path(ace_root, el)
      pred_file <- file.path(el_dir, paste0(el, "_ACE_predictions.csv"))
      if (!file.exists(pred_file)) next
      
      ace_df <- read.csv(pred_file)
      
      obs_col <- paste0(el, "_ICP")
      sd_col  <- paste0(el, "_ICP_sd")
      
      if (!obs_col %in% names(ICP_ppm)) next
      
      # Merge predictions with observed ICP
      joined <- ace_df %>%
        left_join(ICP_ppm %>% dplyr::select(Site, depth, SH20_age, !!sym(obs_col), !!sym(sd_col)),
                  by = c("Site", "depth", "SH20_age"))
      
      for (mname in names(model_store[[el]])) {
        
        pred_nm <- paste0(mname, "_", el, "_Pred_ppm")
        if (!pred_nm %in% names(joined)) next
        
        df <- joined %>%
          dplyr::filter(!is.na(.data[[pred_nm]]),
                        !is.na(.data[[obs_col]]))
        
        if (nrow(df) == 0) next
        
        # Pred vs Obs
        p1 <- ggplot(df, aes_string(x = obs_col, y = pred_nm, colour = "Site")) +
          geom_point(size = 1) +
          geom_abline(slope = 1, intercept = 0) +
          theme_bw() + theme_small +
          labs(
            x = paste(el, "(ICP ppm)"),
            y = paste(mname, "(ACE ppm)"),
            title = paste(el, mname, "- ACE Pred vs Obs")
          )
        
        ggsave(
          file.path(el_dir, paste0(el, "_", mname, "_ACE_PredVsObs.pdf")),
          p1, width = 12, height = 9, units = "cm"
        )
        
        # Residuals
        df$resid <- df[[pred_nm]] - df[[obs_col]]
        p2 <- ggplot(df, aes(resid)) +
          geom_histogram(bins = 30, fill = "grey80") +
          geom_density(colour = "red") +
          theme_bw() + theme_small +
          labs(title = paste(el, mname, "- ACE Residuals"))
        
        ggsave(
          file.path(el_dir, paste0(el, "_", mname, "_ACE_Residuals.pdf")),
          p2, width = 12, height = 9, units = "cm"
        )
      }
    }
  }
  
  # ===============================================================
  # 17. ACE GLOBAL SUMMARY (like Section 11)
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== ACE Global Summary ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    summary_file <- file.path(ace_root, "AllElements_ModelSummary_pred.csv")
    best_file    <- file.path(ace_root, "BestModels_PerElement_pred.csv")
    
    ace_global <- data.frame()
    
    for (el in elements) {
      el_dir <- file.path(ace_root, el)
      diag_file <- file.path(el_dir, paste0(el, "_ACE_diagnostics.csv"))
      if (!file.exists(diag_file)) next
      
      d <- read.csv(diag_file)
      if (nrow(d) == 0) next
      
      d$Element <- el
      temp <- d %>%
        dplyr::arrange(desc(R2), RMSE) %>%
        dplyr::mutate(Rank = dplyr::row_number())
      
      ace_global <- rbind(ace_global, temp)
    }
    
    write.csv(ace_global, summary_file, row.names = FALSE)
    
    best_models <- ace_global %>% dplyr::filter(Rank == 1)
    write.csv(best_models, best_file, row.names = FALSE)
  }
  
  # ===============================================================
  # 18. ACE ppm-scale summary
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== ACE ppm-scale summary ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    summary_file <- file.path(ace_root, "AllElements_ModelSummary_ppm_pred.csv")
    
    if (!exists("ace_global") || nrow(ace_global) == 0) next
    
    # compute means of ICP ppm for scaling
    elem_means_list <- lapply(elements, function(el) {
      obs_col <- paste0(el, "_ICP")
      if (!(obs_col %in% names(ICP_ppm))) return(
        data.frame(Element = el, Mean_ppm = NA_real_)
      )
      
      vals <- ICP_ppm[[obs_col]]
      vals <- vals[!is.na(vals)]
      data.frame(Element = el, Mean_ppm = mean(vals))
    })
    
    elem_means_df <- do.call(rbind, elem_means_list)
    
    ace_ppm_summary <- ace_global %>%
      dplyr::left_join(elem_means_df, by = "Element") %>%
      dplyr::mutate(
        RMSE_factor  = RMSE / Mean_ppm,
        RMSE_abs_ppm = RMSE
      )
    
    write.csv(ace_ppm_summary, summary_file, row.names = FALSE)
    
    message("=== ACE diagnostics + summaries completed successfully ===")
  }

} # end run_full_regressions()
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Run it all ... --------------------------------------------------------

run_full_regressions(
  elements = c("Ca", "Ti","Fe","Mn","Sr","Zr"),
  data     = ACE_dataset,   # calibration (log-space)
  XRF_new  = XRF_pred,      # predictors (log-space)
  ICP_ppm  = ICP_obs,       # observed ICP data in ppm: Ti_ICP, Ti_ICP_sd, Site, SH20_age, etc.
  save_dir = "/Users/sjro/Desktop"  # function will create /Desktop/Regression_Multivariate
)

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Check it works ... ------------------------------------------------------
list.files("/Users/sjro/Desktop/Regression_multivariate_6/All_models_ppm/Ti", full.names = TRUE)
read.csv("/Users/sjro/Desktop/Regression_multivariate_6//All_models_ppm/Ti/Ti_Model_Ranking.csv")
# -------------------------------------------------------------------------


