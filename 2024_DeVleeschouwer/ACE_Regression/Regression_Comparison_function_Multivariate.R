# -------------------------------------------------------------------------
# Regression Comparison Multivariate function for XRF-CS log_inc dataset
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
# - OLS, WLS (var), weighted OLS (sd), 2-step WLS, Bayes, RF, PLS LOO & k-fold
# - R2 (true) as Sum of Squares = 1-(SSres / SS_tot)
# - writes per-element folders into your chosen directory
# Saves:
# - Residual distribution plots (PDF)
# - Predicted vs Observed plots (PDF)
# - PLS RMSEP validation plots (PDF)
# - Influence plots (Cook’s distance + residuals vs leverage) as PDF
# Writes a per-element diagnostics to a .csv file with:
# - AIC & BIC tests
# - Normality test (Shapiro–Wilk)
# - Heteroscedasticity (Breusch–Pagan) for LM/Bayes
# - Unified 10-fold CV RMSE for every model (where possible)
# Native CV metrics for:
# - RF OOB MSE & R² 
# - PLS RMSEP at best ncomp
# -	Explanations for each test in plain English
# Makes a Best Model performance summary
# - ranking based on R2, RMSE, RMSEP where AIC & BIC dont exist (PLS, RF)
# - where R2 higher, RMSE/RMSEP lower, AIC, BIC lower is better
# — sequential R that works on a Mac.

# Key features for each element:
# OLS, WLS, weighted OLS, 2-step weighted WLS
# Bayesian regression (bayesglm)
# Random Forest
# PLS (LOO and 10-fold CV)
# Residual plots, normality tests, heteroscedasticity tests, influence plots
# Predicted vs observed plots colored by Site using a JCO palette
# RF variable importance, PLS RMSEP validation plots
# Per-model/per-element diagnostics (PDF + CSV)
# Performance summary CSV (R², RMSE, RMSEP, AIC, BIC)
# Best-model table
# README.txt file expaling criteria and rankings

# Multivariate run set up
# 
# Single-predictor models (unchanged):
# •	OLS
# •	WLS (variance weighted OLS)
# •	OLS_weighted
# •	WLS_weighted
# 
# Multi-predictor models (updated from univariate function):
# •	Bayes → uses ALL 6 elements as predictors
# •	PLS_LOO → uses ALL 6 elements as predictors
# •	PLS_kfold → uses ALL 6 elements as predictors
# •	Random Forest → uses ALL 6 elements as predictors
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
# Based on original manual code by Steve Roberts
# Refined and made into a function with assistance from ChatGPT
# We acknowledge the use of ChatGPT (OpenAI, version 5.1) for assistance in 
# generating preliminary code and with constructing functions based on manual 
# code written by the authors. The final implementations of the code and 
# functions were produced and verified independently by the authors.

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
# Import ACE matched XRF-CS log (element/inc) file & checks names/data  ----------------------------
ACE_dataset <- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/ACE_Regression/Input/ACE_subsample_xrf_icp_matched_log_inc.csv",
  name_repair = "minimal") %>%   # equivalent to check.names = FALSE 
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>% # Replace infinite values with NA
  #mutate(Site = as.factor(Site)) %>%  # Convert Site to factor
  print()

# Check dataset is loaded OK & working before running function
names(ACE_dataset)
c("Ti", "Ti_ICP") %in% names(ACE_dataset)
c("Ca", "Ca_ICP") %in% names(ACE_dataset)
c("Fe", "Fe_ICP") %in% names(ACE_dataset)
summary(ACE_dataset$Ti)
var(ACE_dataset$Ti)
summary(ACE_dataset$Ti_ICP)
var(ACE_dataset$Ti_ICP)
cor(ACE_dataset$Ti, ACE_dataset$Ti_ICP)



# Import XRF_new for new ppm predictions from calibration ----------------------------------------

main_elements_xrf <- c("Ca", "Fe", "Mn", "Sr", "Ti", "Zr")

XRF_new <- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_log_inc.csv"
) %>%
  dplyr::select(Site, depth, SH20_age, dplyr::all_of(main_elements_xrf)) %>%
  dplyr::filter(Site %in% c("BI10","HER42PB","KER1","KER3","PB1")) %>%
  dplyr::mutate(across(all_of(main_elements_xrf), ~ ifelse(. <= -10, NA, .)))
XRF_new
str(XRF_new)

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Regression comparison function - multivariate 6 element --------------------

run_full_regressions <- function(
    elements,
    data,          # calibration dataset (log-transformed ICP)
    XRF_new = NULL, # new XRF dataset (log-scale predictors)
    var_suffix    = "_var",
    sd_suffix     = "_ICP_sd",
    save_dir      = "/Users/sjro/Desktop/Regression_multivariate_6",
    summary_all_csv = "AllElements_ModelSummary.csv",
    best_models_csv  = "BestModels_PerElement.csv"
) {
  
  ############################################################
  # 1. LOAD PACKAGES
  ############################################################
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
  })
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Base folders for all-model vs best-model outputs
  all_base  <- file.path(save_dir, "All_models_ppm")
  best_base <- file.path(save_dir, "Best_models_ppm")
  dir.create(all_base,  recursive = TRUE, showWarnings = FALSE)
  dir.create(best_base, recursive = TRUE, showWarnings = FALSE)
  
  ############################################################
  # 2. SITE FACTOR IN CALIBRATION DATA
  ############################################################
  if (!"Site" %in% names(data))
    stop("Calibration dataset must contain a 'Site' column.")
  
  data$Site <- factor(
    data$Site,
    levels = c("BI10","HER42PB","KER1","KER3","PB1")
  )
  
  ############################################################
  # 3. UNIVERSAL 8-POINT THEME
  ############################################################
  theme_small <- theme(
    text       = element_text(size = 8),
    axis.text  = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.title   = element_text(size = 8, face = "bold")
  )
  
  ############################################################
  # 4. MODEL LEVELS / LABELS AND COLOURS (NO VarWLS)
  ############################################################
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
  
  # npg colours for the 8 models
  cols_models <- ggsci::pal_npg("nrc")(length(model_labels))
  names(cols_models) <- model_labels
  
  # Force PLS(k) to be salmon3 for contrast with ICP-MS red
  cols_models["PLS(k)"] <- "salmon3"
  
  # Add ICP-MS (Observed) as ninth legend entry
  cols_all <- c(cols_models, "ICP-MS (Observed)" = "blue")
  
  ############################################################
  # 5. DATA INTEGRITY CHECK
  ############################################################
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
  
  ############################################################
  # 6. 10-FOLD CV RMSE HELPER FUNCTION (log-space)
  ############################################################
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
          init_sd <- lm(formula, data = train, weights = 1/(train[[sd_col]]^2))
          rm_sd   <- lm(abs(init_sd$residuals) ~ init_sd$fitted.values)
          w2_sd   <- 1/(rm_sd$fitted.values^2)
          lm(formula, data = train, weights = w2_sd)
          
        } else if (mname == "WLS") {
          ols_fit <- lm(formula, data = train)
          r_ols   <- residuals(ols_fit)
          rm_ols  <- lm(abs(r_ols) ~ fitted(ols_fit))
          w_ols   <- 1/(rm_ols$fitted.values^2)
          lm(formula, data = train, weights = w_ols)
          
        } else if (mname == "Bayes") {
          bayesglm(formula, data = train)
          
        } else if (mname == "RF") {
          randomForest(formula, data = train, ntree = 300)
          
        } else if (mname %in% c("PLS_LOO","PLS_kfold")) {
          plsr(formula, data = train, scale = TRUE, validation = "none")
          
        } else NULL
        
      }, error = function(e) NULL)
      
      if (is.null(fit)) return(NA_real_)
      
      pred <- tryCatch({
        if (inherits(fit, "mvr"))
          as.numeric(predict(fit, ncomp = best_ncomp, newdata = test))
        else
          predict(fit, newdata = test)
      }, error = function(e) NULL)
      
      if (is.null(pred)) return(NA_real_)
      
      y_test <- test[[y_col]]
      ok <- !is.na(y_test) & !is.na(pred)
      if (!any(ok)) next
      
      cv_sse <- cv_sse + sum((y_test[ok] - pred[ok])^2)
      cv_n   <- cv_n   + sum(ok)
    }
    
    if (cv_n == 0) return(NA_real_)
    sqrt(cv_sse / cv_n)
  }
  
  ############################################################
  # 7. EXPLANATION TEXT
  ############################################################
  expl_rank <- paste(
    "Ranking: R² highest → lowest; RMSEP lowest → highest; RMSE lowest → highest.",
    "AIC/BIC reported but not used for ranking."
  )
  
  expl_shapiro <- "Shapiro–Wilk tests residual normality; p > 0.05 desirable."
  expl_bp      <- "Breusch–Pagan tests heteroscedasticity; p < 0.05 indicates non-constant variance."
  expl_cv10    <- "10-fold CV RMSE (log-space); lower indicates better predictive performance."
  expl_r2      <- "R² = 1 − SS_res/SS_tot (log-space)."
  expl_rmse    <- "RMSE = calibration error in log-space."
  expl_rmsep   <- "RMSEP = 10-fold cross-validated error in log-space."
  
  ############################################################
  # 8. GLOBAL SUMMARY TABLE (log-space)
  ############################################################
  global_summary <- tibble(
    Element = character(),
    Model   = character(),
    R2      = numeric(),
    RMSE    = numeric(),   # log-space
    RMSEP   = numeric(),   # log-space
    AIC     = numeric(),
    BIC     = numeric()
  )
  
  ############################################################
  # 9. STORE PER-ELEMENT PREDICTION DATA
  ############################################################
  preds_store <- list()
  
  ############################################################
  # 10. LOOP OVER ELEMENTS
  ############################################################
  
  pb_elements <- progress_bar$new(
    format = "Processing Elements [:bar] :current/:total (:percent) eta::eta",
    total  = length(elements),
    clear  = FALSE,
    width  = 80
  )
  
  for (el in elements) {
    
    pb_elements$tick()
    
    y_col  <- paste0(el, "_ICP")
    x_col  <- el
    var_col <- paste0(el, var_suffix)
    sd_col  <- paste0(el, sd_suffix)
    
    # All-model element folder
    element_dir <- file.path(all_base, el)
    dir.create(element_dir, recursive = TRUE, showWarnings = FALSE)
    
    errors <- check_dataset_integrity(data, el, y_col, x_col, var_col, sd_col)
    if (length(errors) > 0) {
      warning(paste("Skipping", el, ":", paste(errors, collapse = "; ")))
      next
    }
    
    # SINGLE-PREDICTOR FORMULA (OLS / WLS)
    fmla_single <- as.formula(paste0(y_col, " ~ ", x_col))
    
    # MULTI-PREDICTOR FORMULA (Bayes, RF, PLS)
    fmla_multi  <- as.formula(
      paste0(y_col, " ~ ", paste(elements, collapse = " + "))
    )
    
    y <- data[[y_col]]
    
    ############################################################
    # 10.1 FIT MODELS
    ############################################################
    
    model_list <- list()
    
    # OLS
    model_list$OLS <- lm(fmla_single, data = data)
    
    # WLS (residual-based)
    {
      ols_fit <- model_list$OLS
      r_ols   <- residuals(ols_fit)
      rm_ols  <- lm(abs(r_ols) ~ fitted(ols_fit))
      w_ols   <- 1/(rm_ols$fitted.values^2)
      model_list$WLS <- lm(fmla_single, data = data, weights = w_ols)
    }
    
    # OLS_weighted & WLS_weighted (sd-based)
    if (sd_col %in% names(data)) {
      w_sd <- 1/(data[[sd_col]]^2)
      model_list$OLS_weighted <- lm(fmla_single, data = data, weights = w_sd)
      
      init_sd <- lm(fmla_single, data = data, weights = w_sd)
      rm_sd   <- lm(abs(init_sd$residuals) ~ init_sd$fitted.values)
      w2_sd   <- 1/(rm_sd$fitted.values^2)
      model_list$WLS_weighted <- lm(fmla_single, data = data, weights = w2_sd)
    }
    
    # Bayesian regression (multivariate)
    model_list$Bayes <- bayesglm(fmla_multi, data = data)
    
    # Random Forest (multivariate)
    set.seed(42)
    model_list$RF <- randomForest(
      fmla_multi,
      data = data,
      ntree = 500,
      importance = TRUE
    )
    
    # PLS LOOCV (multivariate)
    pls_loo <- plsr(
      fmla_multi,
      data = data,
      scale = TRUE,
      validation = "LOO"
    )
    best_loo <- selectNcomp(pls_loo, method = "onesigma", plot = FALSE)
    attr(pls_loo, "best_ncomp") <- best_loo
    model_list$PLS_LOO <- pls_loo
    
    # PLS k-fold (multivariate)
    pls_k <- plsr(
      fmla_multi,
      data = data,
      scale = TRUE,
      validation = "CV",
      segments = 10
    )
    best_k <- selectNcomp(pls_k, method = "onesigma", plot = FALSE)
    attr(pls_k, "best_ncomp") <- best_k
    model_list$PLS_kfold <- pls_k
    
    ############################################################
    # 10.2 RESULT TABLES FOR THIS ELEMENT
    ############################################################
    
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
    
    ############################################################
    # 10.3 PREDICTION DATAFRAME STARTS AS XRF_new
    ############################################################
    
    do_predictions <- !is.null(XRF_new)
    preds_df <- NULL
    
    if (do_predictions) {
      preds_df <- XRF_new   # retain original columns
    }
    
    ############################################################
    # 10.4 MODEL ANALYSIS LOOP
    ############################################################
    
    pb_models <- progress_bar$new(
      format = paste0("  Models for ", el, " [:bar] :current/:total (:percent)"),
      total  = length(model_list),
      clear  = FALSE,
      width  = 70
    )
    
    for (mname in names(model_list)) {
      
      pb_models$tick()
      m <- model_list[[mname]]
      
      # Predictions on calibration data (log-space)
      pred <- tryCatch({
        if (inherits(m, "mvr"))
          as.numeric(predict(m, ncomp = attr(m, "best_ncomp"), newdata = data))
        else
          predict(m)
      }, error = function(e) NA_real_)
      
      if (any(is.na(pred))) next
      
      resid_vals <- y - pred
      
      # R²
      SS_res <- sum((y - pred)^2)
      SS_tot <- sum((y - mean(y))^2)
      R2     <- 1 - SS_res / SS_tot
      
      # RMSE (log-space)
      RMSE <- sqrt(mean((y - pred)^2))
      
      # RMSEP (log-space)
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
      
      # AIC / BIC
      AIC_val <- if (inherits(m, c("lm","glm","bayesglm")))
        tryCatch(AIC(m), error=function(e) NA_real_) else NA_real_
      BIC_val <- if (inherits(m, c("lm","glm","bayesglm")))
        tryCatch(BIC(m), error=function(e) NA_real_) else NA_real_
      
      # Store in per-element ranking table
      ranking_df <- add_row(
        ranking_df,
        Model = mname,
        R2    = R2,
        RMSE  = RMSE,
        RMSEP = RMSEP_val,
        AIC   = AIC_val,
        BIC   = BIC_val
      )
      
      # Store in global summary
      global_summary <- add_row(
        global_summary,
        Element = el,
        Model   = mname,
        R2      = R2,
        RMSE    = RMSE,
        RMSEP   = RMSEP_val,
        AIC     = AIC_val,
        BIC     = BIC_val
      )
      
      ############################################################
      # Diagnostics
      ############################################################
      
      sh_W <- sh_p <- NA_real_
      if (length(resid_vals) >= 3 && length(resid_vals) <= 5000) {
        sh <- tryCatch(shapiro.test(resid_vals), error=function(e) NULL)
        if (!is.null(sh)) {
          sh_W <- sh$statistic
          sh_p <- sh$p.value
        }
      }
      
      bp_stat <- bp_p <- NA_real_
      if (inherits(m, c("lm","bayesglm"))) {
        bp <- tryCatch(bptest(m), error=function(e) NULL)
        if (!is.null(bp)) {
          bp_stat <- bp$statistic
          bp_p    <- bp$p.value
        }
      }
      
      diag_df <- add_row(
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
      
      ############################################################
      # PREDICTIONS ON XRF_new (log-space → ppm, stored as *_ppm cols)
      ############################################################
      
      if (do_predictions) {
        
        pred_new_log <- tryCatch({
          if (inherits(m, "mvr"))
            as.numeric(predict(m, ncomp = attr(m, "best_ncomp"), newdata = XRF_new))
          else
            as.numeric(predict(m, newdata = XRF_new))
        }, error = function(e) rep(NA_real_, nrow(XRF_new)))
        
        lower_log <- upper_log <- pred_new_log
        
        if (mname %in% c("OLS","WLS","OLS_weighted","WLS_weighted")) {
          
          pi_obj <- tryCatch(
            predict(m, newdata = XRF_new, interval = "prediction", level = 0.95),
            error = function(e) NULL
          )
          if (!is.null(pi_obj)) {
            lower_log <- as.numeric(pi_obj[, "lwr"])
            upper_log <- as.numeric(pi_obj[, "upr"])
          } else {
            lower_log <- pred_new_log - 1.96 * RMSE
            upper_log <- pred_new_log + 1.96 * RMSE
          }
          
        } else if (mname == "Bayes") {
          
          p_bayes <- tryCatch(
            predict(m, newdata = XRF_new, se.fit = TRUE, type = "link"),
            error = function(e) NULL
          )
          if (!is.null(p_bayes)) {
            fit <- as.numeric(p_bayes$fit)
            se  <- as.numeric(p_bayes$se.fit)
            pred_new_log <- fit
            lower_log <- fit - 1.96 * se
            upper_log <- fit + 1.96 * se
          } else {
            lower_log <- pred_new_log - 1.96 * RMSE
            upper_log <- pred_new_log + 1.96 * RMSE
          }
          
        } else if (mname %in% c("RF","PLS_LOO","PLS_kfold")) {
          lower_log <- pred_new_log - 1.96 * RMSE
          upper_log <- pred_new_log + 1.96 * RMSE
        }
        
        pred_ppm  <- exp(pred_new_log)
        lower_ppm <- exp(lower_log)
        upper_ppm <- exp(upper_log)
        
        base_name <- paste(mname, el, sep = "_")
        
        pred_ppm_col <- paste0(base_name, "_Pred_ppm")
        lwr_ppm_col  <- paste0(base_name, "_L95_ppm")
        upr_ppm_col  <- paste0(base_name, "_U95_ppm")
        
        preds_df[[pred_ppm_col]] <- pred_ppm
        preds_df[[lwr_ppm_col]]  <- lower_ppm
        preds_df[[upr_ppm_col]]  <- upper_ppm
      }
      
      ############################################################
      # PLOT LABEL WITH R2 / RMSEP / RMSE
      ############################################################
      stats_label <- paste0(
        "R² = ", signif(R2, 3),
        "\nRMSEP = ", signif(RMSEP_val, 3),
        "\nRMSE = ", signif(RMSE, 3)
      )
      
      ############################################################
      # 8-POINT PLOTS — 13 × 9 cm (Pred vs Obs, Residuals)
      ############################################################
      
      # 1) Predicted vs Observed (log space, unitless)
      df_pred <- data.frame(
        Observed  = y,
        Predicted = pred,
        Site      = data$Site
      )
      
      x_min_po <- min(df_pred$Observed,  na.rm = TRUE)
      y_max_po <- max(df_pred$Predicted, na.rm = TRUE)
      
      p_po <- ggplot(df_pred, aes(Observed, Predicted, color = Site)) +
        geom_point(size = 1.25, alpha = 0.85) +
        ggsci::scale_color_npg(name = "Site") +
        geom_abline(slope = 1, intercept = 0) +
        annotate(
          "text",
          x = x_min_po,
          y = y_max_po,
          hjust = 0,
          vjust = 1,
          label = stats_label,
          size = 2.7
        ) +
        theme_bw() + theme_small +
        ggtitle(paste(el, mname, "- Predicted vs Observed (log-space)"))
      
      ggsave(
        file.path(element_dir, paste0(el, "_", mname, "_PredVsObs.pdf")),
        p_po, width = 13, height = 9, units = "cm"
      )
      
      # 2) Residual distribution (log-space)
      res_df <- data.frame(resid = resid_vals)
      hist_obj <- hist(resid_vals, plot = FALSE)
      x_min_res <- min(resid_vals, na.rm = TRUE)
      y_max_res <- max(hist_obj$counts, na.rm = TRUE)
      
      p_res <- ggplot(res_df, aes(resid)) +
        geom_histogram(bins = 30, fill = "grey80") +
        geom_density(color = "red") +
        annotate(
          "text",
          x = x_min_res,
          y = y_max_res,
          hjust = 0,
          vjust = 1,
          label = stats_label,
          size = 2.7
        ) +
        theme_bw() + theme_small +
        ggtitle(paste(el, mname, "- Residuals (log-space)"))
      
      ggsave(
        file.path(element_dir, paste0(el, "_", mname, "_Residuals.pdf")),
        p_res, width = 13, height = 9, units = "cm"
      )
      
      # 3) Influence plots (lm only)
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
            ggtitle(paste(el, mname, "- Leverage Plot"))
        )
        
        dev.off()
      }
      
      # 4) PLS RMSEP vs Components
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
    
    ############################################################
    # 10.5 RANK MODELS — R2 ↓ → RMSEP ↑ → RMSE ↑
    ############################################################
    
    if (nrow(ranking_df) > 0) {
      
      ranking_df <- ranking_df %>%
        arrange(
          dplyr::desc(R2),
          RMSEP,
          RMSE
        ) %>%
        mutate(Rank = dplyr::row_number())
      
      ranking_expl_df <- ranking_df %>%
        mutate(
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
    
    ############################################################
    # 10.6 DIAGNOSTICS SUMMARY
    ############################################################
    
    write.csv(
      diag_df,
      file.path(element_dir, paste0(el, "_Diagnostics_Summary.csv")),
      row.names = FALSE
    )
    
    ############################################################
    # 10.7 SAVE PREDICTION CSV (PPM) — <Element>_preds.csv
    ############################################################
    
    if (!is.null(XRF_new) && !is.null(preds_df)) {
      preds_path <- file.path(element_dir, paste0(el, "_preds.csv"))
      write.csv(preds_df, preds_path, row.names = FALSE)
      
      # Also keep in memory for site-level & element-level plots
      preds_store[[el]] <- preds_df
    }
    
    ############################################################
    # 10.8 README PER ELEMENT
    ############################################################
    
    readme_txt <- c(
      paste("Element:", el),
      "",
      "DATA NOTE:",
      "  - Calibration responses (e.g., ", el, "_ICP) are natural-log transformed (ln mg kg^-1).",
      "  - All model fitting and errors (RMSE, RMSEP) are computed in log-space.",
      "  - Predictions on XRF_new are back-transformed with exp() to mg kg^-1.",
      "",
      "PREDICTION COLUMNS IN *_preds.csv:",
      "  - For each model, columns follow the pattern:",
      "      Model_Element_Pred_ppm : predicted concentration in mg kg^-1",
      "      Model_Element_L95_ppm  : lower 95% prediction bound in mg kg^-1",
      "      Model_Element_U95_ppm  : upper 95% prediction bound in mg kg^-1",
      "    e.g. PLS_LOO_Ti_Pred_ppm, PLS_LOO_Ti_L95_ppm, PLS_LOO_Ti_U95_ppm.",
      "",
      "MODELS FITTED:",
      "  - OLS, WLS, OLS(wt), WLS(wt), Bayes, RF, PLS(LOO), PLS(k).",
      "",
      "RANKING RULES:",
      "  1) Highest R² ranked best.",
      "  2) For similar R², lower RMSEP ranked higher.",
      "  3) For similar R² and RMSEP, lower RMSE ranked higher."
    )
    
    writeLines(readme_txt, file.path(element_dir, "README.txt"))
    
  } # end element loop
  
  
  ############################################################
  # 11. SITE-LEVEL A4 PLOTS (ALL & BEST MULTIVARIATE)
  ############################################################
  
  if (!is.null(XRF_new) && length(preds_store) > 0) {
    
    has_site     <- "Site" %in% names(XRF_new)
    has_depth    <- "depth" %in% names(XRF_new)
    has_age_new  <- "SH20_age" %in% names(XRF_new)
    has_age_cal  <- "SH20_mean_age" %in% names(data)
    
    sites_dir_all  <- file.path(all_base,  "Sites")
    sites_dir_best <- file.path(best_base, "Sites")
    dir.create(sites_dir_all,  recursive = TRUE, showWarnings = FALSE)
    dir.create(sites_dir_best, recursive = TRUE, showWarnings = FALSE)
    
    site_levels <- if (has_site) sort(unique(XRF_new$Site)) else character(0)
    
    multivar_labels <- c("Bayes","RF","PLS(LOO)","PLS(k)")
    
    for (s in site_levels) {
      
      ##########################################################
      # 11.1 DEPTH-BASED A4 PLOTS PER SITE (2×3 GRID)
      ##########################################################
      
      depth_plots_all  <- list()
      depth_plots_best <- list()
      
      if (has_depth) {
        for (el in elements) {
          
          if (is.null(preds_store[[el]])) next
          
          df_el <- preds_store[[el]]
          if (!all(c("Site", "depth") %in% names(df_el))) next
          
          df_site <- df_el %>%
            dplyr::filter(Site == s)
          
          if (nrow(df_site) == 0) next
          
          pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                            names(df_site), value = TRUE)
          if (length(pred_cols) == 0) next
          
          df_long <- df_site %>%
            dplyr::select(depth, Site, all_of(pred_cols)) %>%
            tidyr::pivot_longer(
              cols = all_of(pred_cols),
              names_to = "Model",
              values_to = "Pred_ppm"
            )
          
          df_long$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long$Model)
          df_long$Model <- factor(df_long$Model,
                                  levels = model_levels,
                                  labels = model_labels)
          
          df_long_best <- df_long %>%
            dplyr::filter(Model %in% multivar_labels)
          
          y_col <- paste0(el, "_ICP")
          ace_overlay <- NULL
          if (all(c("Site", "depth", y_col) %in% names(data))) {
            ace_overlay <- data %>%
              dplyr::filter(Site == s) %>%
              dplyr::filter(!is.na(.data[[y_col]]), !is.na(depth)) %>%
              dplyr::mutate(ICP_ppm = exp(.data[[y_col]]))
            if (nrow(ace_overlay) == 0) ace_overlay <- NULL
          }
          
          # All models depth profile
          p_depth_all <- ggplot(df_long, aes(x = Pred_ppm, y = depth,
                                             colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_y_reverse() +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Depth (cm)",
              title = paste0("Site ", s, " — ", el, " (depth, all models)")
            ) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = names(cols_all),
              guide  = guide_legend(nrow = 2)
            ) +
            theme_bw() + theme_small
          
          # Best multivariate depth
          p_depth_best <- ggplot(df_long_best, aes(x = Pred_ppm, y = depth,
                                                   colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_y_reverse() +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Depth (cm)",
              title = paste0("Site ", s, " — ", el, " (depth, multivariate)")
            ) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = c(multivar_labels, "ICP-MS (Observed)"),
              guide  = guide_legend(nrow = 2)
            ) +
            theme_bw() + theme_small
          
          if (!is.null(ace_overlay)) {
            p_depth_all <- p_depth_all +
              geom_path(
                data = ace_overlay,
                aes(x = ICP_ppm, y = depth,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_overlay,
                aes(x = ICP_ppm, y = depth,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                fill = "white",
                shape = 21,
                stroke = 0.4,
                size = 1.4
              )
            
            p_depth_best <- p_depth_best +
              geom_path(
                data = ace_overlay,
                aes(x = ICP_ppm, y = depth,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_overlay,
                aes(x = ICP_ppm, y = depth,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                fill = "white",
                shape = 21,
                stroke = 0.4,
                size = 1.4
              )
          }
          
          depth_plots_all[[el]]  <- p_depth_all
          depth_plots_best[[el]] <- p_depth_best
        }
        
        if (length(depth_plots_all) > 0) {
          
          ord_all <- lapply(elements, function(e) depth_plots_all[[e]])
          ord_all <- lapply(ord_all, function(p) {
            if (is.null(p)) ggplot() + theme_void() else p
          })
          
          combined_depth_all <- (ord_all[[1]] | ord_all[[2]]) /
            (ord_all[[3]] | ord_all[[4]]) /
            (ord_all[[5]] | ord_all[[6]])
          
          combined_depth_all <- combined_depth_all +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          out_depth_all <- file.path(sites_dir_all,
                                     paste0("Site_", s, "_Profiles_A4_depth_all.pdf"))
          
          ggsave(
            out_depth_all,
            combined_depth_all,
            width = 21, height = 29.7, units = "cm"
          )
        }
        
        if (length(depth_plots_best) > 0) {
          
          ord_best <- lapply(elements, function(e) depth_plots_best[[e]])
          ord_best <- lapply(ord_best, function(p) {
            if (is.null(p)) ggplot() + theme_void() else p
          })
          
          combined_depth_best <- (ord_best[[1]] | ord_best[[2]]) /
            (ord_best[[3]] | ord_best[[4]]) /
            (ord_best[[5]] | ord_best[[6]])
          
          combined_depth_best <- combined_depth_best +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          out_depth_best <- file.path(sites_dir_best,
                                      paste0("Site_", s, "_Profiles_A4_depth_multivariate.pdf"))
          
          ggsave(
            out_depth_best,
            combined_depth_best,
            width = 21, height = 29.7, units = "cm"
          )
        }
      }
      
      ##########################################################
      # 11.2 AGE-BASED A4 PLOTS PER SITE (Y = Age, X = mg kg^-1)
      ##########################################################
      
      age_plots_all  <- list()
      age_plots_best <- list()
      
      if (has_age_new && has_age_cal) {
        
        for (el in elements) {
          
          if (is.null(preds_store[[el]])) next
          
          df_el <- preds_store[[el]]
          if (!all(c("Site", "SH20_age") %in% names(df_el))) next
          
          df_site <- df_el %>%
            dplyr::filter(Site == s)
          
          if (nrow(df_site) == 0) next
          
          pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                            names(df_site), value = TRUE)
          if (length(pred_cols) == 0) next
          
          df_long <- df_site %>%
            dplyr::select(SH20_age, Site, all_of(pred_cols)) %>%
            tidyr::pivot_longer(
              cols = all_of(pred_cols),
              names_to = "Model",
              values_to = "Pred_ppm"
            )
          
          df_long$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long$Model)
          df_long$Model <- factor(df_long$Model,
                                  levels = model_levels,
                                  labels = model_labels)
          
          df_long_best <- df_long %>%
            dplyr::filter(Model %in% multivar_labels)
          
          y_col <- paste0(el, "_ICP")
          ace_overlay_age <- NULL
          if (all(c("Site", "SH20_mean_age", y_col) %in% names(data))) {
            ace_overlay_age <- data %>%
              dplyr::filter(Site == s) %>%
              dplyr::filter(!is.na(.data[[y_col]]), !is.na(SH20_mean_age)) %>%
              dplyr::mutate(ICP_ppm = exp(.data[[y_col]]))
            if (nrow(ace_overlay_age) == 0) ace_overlay_age <- NULL
          }
          
          # All models age vs ppm
          p_age_all <- ggplot(df_long, aes(x = Pred_ppm, y = SH20_age,
                                           colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Age (cal a BP)",
              title = paste0("Site ", s, " — ", el, " (age vs mg kg^-1, all)")
            ) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = names(cols_all),
              guide  = guide_legend(nrow = 2)
            ) +
            theme_bw() + theme_small
          
          # Best multivariate age vs ppm
          p_age_best <- ggplot(df_long_best, aes(x = Pred_ppm, y = SH20_age,
                                                 colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Age (cal a BP)",
              title = paste0("Site ", s, " — ", el, " (age vs mg kg^-1, multivariate)")
            ) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = c(multivar_labels, "ICP-MS (Observed)"),
              guide  = guide_legend(nrow = 2)
            ) +
            theme_bw() + theme_small
          
          if (!is.null(ace_overlay_age)) {
            p_age_all <- p_age_all +
              geom_path(
                data = ace_overlay_age,
                aes(x = ICP_ppm, y = SH20_mean_age,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_overlay_age,
                aes(x = ICP_ppm, y = SH20_mean_age,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                fill = "white",
                shape = 21,
                stroke = 0.4,
                size = 1.4
              )
            
            p_age_best <- p_age_best +
              geom_path(
                data = ace_overlay_age,
                aes(x = ICP_ppm, y = SH20_mean_age,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_overlay_age,
                aes(x = ICP_ppm, y = SH20_mean_age,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                fill = "white",
                shape = 21,
                stroke = 0.4,
                size = 1.4
              )
          }
          
          age_plots_all[[el]]  <- p_age_all
          age_plots_best[[el]] <- p_age_best
        }
        
        if (length(age_plots_all) > 0) {
          
          ord_a <- lapply(elements, function(e) age_plots_all[[e]])
          ord_a <- lapply(ord_a, function(p) {
            if (is.null(p)) ggplot() + theme_void() else p
          })
          
          combined_age_all <- (ord_a[[1]] | ord_a[[2]]) /
            (ord_a[[3]] | ord_a[[4]]) /
            (ord_a[[5]] | ord_a[[6]])
          
          combined_age_all <- combined_age_all +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          out_age_all <- file.path(sites_dir_all,
                                   paste0("Site_", s, "_Profiles_A4_age_all.pdf"))
          
          ggsave(
            out_age_all,
            combined_age_all,
            width = 21, height = 29.7, units = "cm"
          )
        }
        
        if (length(age_plots_best) > 0) {
          
          ord_ab <- lapply(elements, function(e) age_plots_best[[e]])
          ord_ab <- lapply(ord_ab, function(p) {
            if (is.null(p)) ggplot() + theme_void() else p
          })
          
          combined_age_best <- (ord_ab[[1]] | ord_ab[[2]]) /
            (ord_ab[[3]] | ord_ab[[4]]) /
            (ord_ab[[5]] | ord_ab[[6]])
          
          combined_age_best <- combined_age_best +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          out_age_best <- file.path(sites_dir_best,
                                    paste0("Site_", s, "_Profiles_A4_age_multivariate.pdf"))
          
          ggsave(
            out_age_best,
            combined_age_best,
            width = 21, height = 29.7, units = "cm"
          )
        }
      }
    } # end site loop
  }
  
  ############################################################
  # 12. ELEMENT-LEVEL MULTI-SITE AGE–CONCENTRATION COMPARISON
  #     (ALL MODELS & BEST MULTIVARIATE)
  ############################################################
  
  if (!is.null(XRF_new) &&
      length(preds_store) > 0 &&
      "Site" %in% names(XRF_new) &&
      "SH20_age" %in% names(XRF_new) &&
      "SH20_mean_age" %in% names(data)) {
    
    elements_dir_all  <- file.path(all_base,  "elements")
    elements_dir_best <- file.path(best_base, "elements")
    dir.create(elements_dir_all,  recursive = TRUE, showWarnings = FALSE)
    dir.create(elements_dir_best, recursive = TRUE, showWarnings = FALSE)
    
    age_all <- c(XRF_new$SH20_age, data$SH20_mean_age)
    age_all <- age_all[!is.na(age_all)]
    if (length(age_all) > 0) {
      age_min <- min(age_all)
      age_max <- max(age_all)
    } else {
      age_min <- NA_real_
      age_max <- NA_real_
    }
    
    site_levels_all <- sort(unique(c(as.character(XRF_new$Site),
                                     as.character(data$Site))))
    
    multivar_labels <- c("Bayes","RF","PLS(LOO)","PLS(k)")
    
    for (el in elements) {
      
      if (is.null(preds_store[[el]])) next
      
      df_el <- preds_store[[el]]
      
      pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                        names(df_el), value = TRUE)
      if (length(pred_cols) == 0) next
      
      df_long_all <- df_el %>%
        dplyr::filter(!is.na(SH20_age)) %>%
        dplyr::select(Site, SH20_age, all_of(pred_cols)) %>%
        tidyr::pivot_longer(
          cols = all_of(pred_cols),
          names_to = "Model",
          values_to = "Pred_ppm"
        )
      
      df_long_all$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long_all$Model)
      df_long_all$Model <- factor(df_long_all$Model,
                                  levels = model_levels,
                                  labels = model_labels)
      
      df_long_best <- df_long_all %>%
        dplyr::filter(Model %in% multivar_labels)
      
      y_col <- paste0(el, "_ICP")
      ace_all <- NULL
      if (all(c("Site", "SH20_mean_age", y_col) %in% names(data))) {
        ace_all <- data %>%
          dplyr::filter(!is.na(SH20_mean_age), !is.na(.data[[y_col]])) %>%
          dplyr::mutate(ICP_ppm = exp(.data[[y_col]]))
      }
      
      if (is.null(ace_all) || nrow(df_long_all) == 0) next
      
      # ALL models page
      plot_list_all  <- list()
      # BEST multivariate-only page
      plot_list_best <- list()
      
      for (s in site_levels_all) {
        
        df_pred_site_all  <- df_long_all  %>% dplyr::filter(Site == s)
        df_pred_site_best <- df_long_best %>% dplyr::filter(Site == s)
        
        ace_site <- ace_all %>%
          dplyr::filter(Site == s)
        
        # ALL models
        if (nrow(df_pred_site_all) == 0 && nrow(ace_site) == 0) {
          plot_list_all[[s]] <- ggplot() + theme_void()
        } else {
          p_site_all <- ggplot(df_pred_site_all,
                               aes(x = SH20_age, y = Pred_ppm,
                                   colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste0("Site ", s, " (all models)")
            ) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = names(cols_all),
              guide  = guide_legend(nrow = 2)
            ) +
            theme_bw() + theme_small
          
          if (!is.na(age_min) && !is.na(age_max)) {
            p_site_all <- p_site_all + coord_cartesian(xlim = c(age_min, age_max))
          }
          
          if (nrow(ace_site) > 0) {
            p_site_all <- p_site_all +
              geom_path(
                data = ace_site,
                aes(x = SH20_mean_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_mean_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                fill = "white",
                shape = 21,
                stroke = 0.4,
                size = 1.4
              )
          }
          
          plot_list_all[[s]] <- p_site_all
        }
        
        # BEST multivariate
        if (nrow(df_pred_site_best) == 0 && nrow(ace_site) == 0) {
          plot_list_best[[s]] <- ggplot() + theme_void()
        } else {
          p_site_best <- ggplot(df_pred_site_best,
                                aes(x = SH20_age, y = Pred_ppm,
                                    colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste0("Site ", s, " (multivariate)")
            ) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = c(multivar_labels, "ICP-MS (Observed)"),
              guide  = guide_legend(nrow = 2)
            ) +
            theme_bw() + theme_small
          
          if (!is.na(age_min) && !is.na(age_max)) {
            p_site_best <- p_site_best + coord_cartesian(xlim = c(age_min, age_max))
          }
          
          if (nrow(ace_site) > 0) {
            p_site_best <- p_site_best +
              geom_path(
                data = ace_site,
                aes(x = SH20_mean_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_mean_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                fill = "white",
                shape = 21,
                stroke = 0.4,
                size = 1.4
              )
          }
          
          plot_list_best[[s]] <- p_site_best
        }
      }
      
      if (length(plot_list_all) > 0) {
        combined_all <- Reduce(`/`, plot_list_all)
        combined_all <- combined_all +
          plot_layout(guides = "collect") &
          theme(legend.position = "bottom")
        
        out_el_all <- file.path(elements_dir_all,
                                paste0(el, "_model_age_comparison_all.pdf"))
        
        ggsave(
          out_el_all,
          combined_all,
          width = 21, height = 29.7, units = "cm"
        )
      }
      
      if (length(plot_list_best) > 0) {
        combined_best <- Reduce(`/`, plot_list_best)
        combined_best <- combined_best +
          plot_layout(guides = "collect") &
          theme(legend.position = "bottom")
        
        out_el_best <- file.path(elements_dir_best,
                                 paste0(el, "_model_age_comparison_multivariate.pdf"))
        
        ggsave(
          out_el_best,
          combined_best,
          width = 21, height = 29.7, units = "cm"
        )
      }
    }
  }
  
  ############################################################
  # 13. GLOBAL SUMMARY ACROSS ALL ELEMENTS (LOG-SPACE + PPM)
  ############################################################
  
  if (nrow(global_summary) > 0) {
    
    global_ranked <- global_summary %>%
      arrange(
        Element,
        dplyr::desc(R2),
        RMSEP,
        RMSE
      ) %>%
      group_by(Element) %>%
      mutate(Rank = dplyr::row_number()) %>%
      ungroup()
    
    write.csv(
      global_ranked,
      file.path(save_dir, summary_all_csv),
      row.names = FALSE
    )
    
    best_models <- global_ranked %>%
      dplyr::filter(Rank == 1)
    
    write.csv(
      best_models,
      file.path(save_dir, best_models_csv),
      row.names = FALSE
    )
    
    ##########################################################
    # 13.1 PPM-SCALE SUMMARY: AllElements_ModelSummary_ppm.csv
    ##########################################################
    
    elem_means_list <- lapply(elements, function(el) {
      y_col <- paste0(el, "_ICP")
      if (!y_col %in% names(data)) {
        return(data.frame(Element = el,
                          Mean_log_ppm = NA_real_,
                          Mean_ppm = NA_real_))
      }
      vals <- data[[y_col]]
      vals <- vals[!is.na(vals)]
      if (!length(vals)) {
        return(data.frame(Element = el,
                          Mean_log_ppm = NA_real_,
                          Mean_ppm = NA_real_))
      }
      mean_log <- mean(vals)
      mean_ppm <- exp(mean_log)
      data.frame(Element = el,
                 Mean_log_ppm = mean_log,
                 Mean_ppm     = mean_ppm)
    })
    
    elem_means_df <- do.call(rbind, elem_means_list)
    
    global_ranked_ppm <- global_ranked %>%
      left_join(elem_means_df, by = "Element") %>%
      mutate(
        RMSE_log      = RMSE,
        RMSEP_log     = RMSEP,
        RMSE_factor   = exp(RMSE_log),   # multiplicative error factor
        RMSEP_factor  = exp(RMSEP_log),  # multiplicative error factor
        RMSE_abs_ppm  = RMSE_log  * Mean_ppm,  # approximate abs error (mg kg^-1)
        RMSEP_abs_ppm = RMSEP_log * Mean_ppm
      ) %>%
      dplyr::select(
        Element, Model, Rank, R2,
        RMSE_log,  RMSE_factor,  RMSE_abs_ppm,
        RMSEP_log, RMSEP_factor, RMSEP_abs_ppm,
        Mean_log_ppm, Mean_ppm,
        AIC, BIC
      )
    
    write.csv(
      global_ranked_ppm,
      file.path(save_dir, "AllElements_ModelSummary_ppm.csv"),
      row.names = FALSE
    )
  }
  
  ############################################################
  # 14. FINISH
  ############################################################
  
  message("=== COMPLETED SUCCESSFULLY ===")
  message("All results saved in: ", save_dir)
}

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Run it ... --------------------------------------------------------

# Run code
elements <- c("Ti","Ca","Mn","Fe","Sr","Zr")
run_full_regressions(
  elements = c("Ca", "Ti","Mn","Fe","Sr","Zr"),
  data     = ACE_dataset, # matched calibration dataset 
  XRF_new  = XRF_new     # new XRF-CS prediction dataset
)

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Check it works ... ------------------------------------------------------
list.files("/Users/sjro/Desktop/Regression_multivariate_6/All_models_ppm/Ti", full.names = TRUE)
read.csv("/Users/sjro/Desktop/Regression_multivariate_6//All_models_ppm/Ti/Ti_Model_Ranking.csv")
# -------------------------------------------------------------------------


