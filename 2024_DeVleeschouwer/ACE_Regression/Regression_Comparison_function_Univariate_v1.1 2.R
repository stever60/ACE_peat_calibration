# Regression Comparison function for XRF-CS log_inc dataset

# -------------------------------------------------------------------------
# ReadMe description ----------------------------------------

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
# - RF OOB MSE & R² (Random Forest OOB (Out-of-Bag) is a method for estimating a random forest model's performance without using a separate test set. It works by having each tree in the forest evaluate the data points that were not included in its own training sample (the "out-of-bag" data). By averaging these predictions across all trees and their respective OOB samples, the out-of-bag error provides a robust and unbiased estimate of the model's generalization error)
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
# •	WLS_var
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
# Set up - Clear all ------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()
# Install and load libraries needed ---------------------------------------------------------------
install.packages(c(
  "broom", "performance", "see", "ggplot2", "dplyr", "purrr",
  "boot", "psych", "lmtest", "arm", "randomForest", "pls", "ggpubr", 
  "qqplotr", "future", "progressr"
))

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
# ACE_dataset <- read.csv("Papers_R/2024_DeVleeschouwer/ACE_Regression/Input/ACE_subsample_xrf_icp_matched_log_inc.csv", check.names=FALSE)  # safest
# is.na(ACE_subsample_icp_xrf_matched_log_inc)<-sapply(ACE_subsample_icp_xrf_matched_log_inc, is.infinite) # replace infinite values with NA
# #ACE_dataset$Site <- as.factor(ACE_dataset$Site) # Convert Site to grouping variable
# ACE_dataset

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

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Regression comparison - univariate - function ---------------------------------------------------------------------

run_full_regressions <- function(
    elements,
    data,
    var_suffix      = "_var",
    sd_suffix       = "_ICP_sd",
    save_dir        = "/Users/sjro/Desktop/Regression_univariate",
    summary_all_csv = "AllElements_ModelSummary.csv",
    best_models_csv = "BestModels_PerElement.csv"
) {
  
  ############################################################
  # 1. REQUIRED PACKAGES
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
  })
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  ############################################################
  # 2. ENSURE 'Site' COLUMN EXISTS AND IS ORDERED
  ############################################################
  if (!"Site" %in% names(data))
    stop("Input dataset must contain a column named 'Site'.")
  
  data$Site <- factor(
    data$Site,
    levels = c("BI10","HER42PB","KER1","KER3","PB1")
  )
  
  
  ############################################################
  # 3. DATA INTEGRITY CHECKER
  ############################################################
  check_dataset_integrity <- function(data, el, y_col, x_col, var_col, sd_col) {
    
    errors <- c()
    
    required <- c(y_col, x_col)
    for (v in required) {
      if (!v %in% names(data)) errors <- c(errors, paste("Missing column:", v))
    }
    
    if (!is.numeric(data[[y_col]])) errors <- c(errors, paste(y_col, "not numeric"))
    if (!is.numeric(data[[x_col]])) errors <- c(errors, paste(x_col, "not numeric"))
    
    if (any(is.na(data[[y_col]]))) errors <- c(errors, paste("NA values in", y_col))
    if (any(is.na(data[[x_col]]))) errors <- c(errors, paste("NA values in", x_col))
    
    if (var(data[[x_col]], na.rm = TRUE) == 0)
      errors <- c(errors, paste("Predictor", x_col, "has zero variance."))
    
    if (sd_col %in% names(data)) {
      if (!is.numeric(data[[sd_col]]) ||
          any(data[[sd_col]] <= 0, na.rm = TRUE))
        errors <- c(errors, paste(sd_col, "must contain positive numeric values."))
    }
    
    if (var_col %in% names(data)) {
      if (!is.numeric(data[[var_col]]) ||
          any(data[[var_col]] <= 0, na.rm = TRUE))
        errors <- c(errors, paste(var_col, "must contain positive numeric values."))
    }
    
    errors
  }
  
  
  ############################################################
  # 4. CROSS-VALIDATION RMSE FUNCTION
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
      
      # Fit appropriate model for CV
      fit <- tryCatch({
        
        if (mname == "OLS") {
          lm(formula, data = train)
          
        } else if (mname == "WLS_var" && var_col %in% names(train)) {
          lm(formula, data = train, weights = 1/train[[var_col]])
          
        } else if (mname == "OLS_weighted" && sd_col %in% names(train)) {
          lm(formula, data = train, weights = 1/(train[[sd_col]]^2))
          
        } else if (mname == "WLS_weighted" && sd_col %in% names(train)) {
          init_sd <- lm(formula, data = train, weights = 1/(train[[sd_col]]^2))
          rm_sd   <- lm(abs(init_sd$residuals) ~ init_sd$fitted.values)
          w2      <- 1/(rm_sd$fitted.values^2)
          lm(formula, data = train, weights = w2)
          
        } else if (mname == "Bayes") {
          bayesglm(formula, data = train)
          
        } else if (mname == "RF") {
          randomForest(formula, data = train, ntree = 300)
          
        } else if (mname %in% c("PLS_LOO","PLS_kfold")) {
          plsr(formula, data = train, scale = TRUE, validation = "none")
          
        } else NULL
        
      }, error = function(e) NULL)
      
      # Prediction
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
  # 5. EXPLANATION STRINGS
  ############################################################
  expl_rank <- paste(
    "Models are ranked by: (1) R² highest → lowest; (2) RMSEP lowest → highest;",
    "(3) RMSE lowest → highest. RMSEP is emphasised over RMSE because it reflects",
    "predictive performance under cross-validation. AIC/BIC are reported but not",
    "used for ranking."
  )
  
  expl_shapiro <- "Shapiro–Wilk tests residual normality; p > 0.05 desirable."
  expl_bp      <- "Breusch–Pagan tests heteroscedasticity; p < 0.05 suggests non-constant variance."
  expl_cv10    <- "10-fold cross-validation RMSE; lower values indicate better predictive performance."
  expl_r2      <- "R² = 1 − SS_res/SS_tot; fraction of variance explained (higher is better)."
  expl_rmse    <- "RMSE: calibration error (lower is better)."
  expl_rmsep   <- "RMSEP: cross-validated predictive error (lowest values preferred)."
  
  
  ############################################################
  # 6. INITIALISE GLOBAL SUMMARY DATAFRAME
  ############################################################
  global_summary <- tibble(
    Element = character(),
    Model   = character(),
    R2      = numeric(),
    RMSE    = numeric(),
    RMSEP   = numeric(),
    AIC     = numeric(),
    BIC     = numeric()
  )
  
  ############################################################
  # 7. MAIN LOOP OVER ELEMENTS (with progress bar)
  ############################################################
  
  pb_elements <- progress_bar$new(
    format = "Elements [:bar] :current/:total (:percent) eta::eta",
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
    
    element_dir <- file.path(save_dir, el)
    dir.create(element_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Integrity check
    errors <- check_dataset_integrity(data, el, y_col, x_col, var_col, sd_col)
    if (length(errors) > 0) {
      warning(paste("Skipping", el, ":", paste(errors, collapse = "; ")))
      next
    }
    
    fmla <- as.formula(paste0(y_col, " ~ ", x_col))
    y <- data[[y_col]]
    
    ############################################################
    # 7.1 FIT ALL MODELS FOR THIS ELEMENT
    ############################################################
    
    model_list <- list()
    
    # OLS
    model_list$OLS <- lm(fmla, data = data)
    
    # WLS using variance column if present
    if (var_col %in% names(data)) {
      model_list$WLS_var <- lm(fmla, data = data, weights = 1 / data[[var_col]])
    }
    
    # Weighted OLS and 2-step WLS using SD if present
    if (sd_col %in% names(data)) {
      w_sd <- 1 / (data[[sd_col]]^2)
      model_list$OLS_weighted <- lm(fmla, data = data, weights = w_sd)
      
      init_sd <- lm(fmla, data = data, weights = w_sd)
      rm_sd   <- lm(abs(init_sd$residuals) ~ init_sd$fitted.values)
      w2_sd   <- 1 / (rm_sd$fitted.values^2)
      model_list$WLS_weighted <- lm(fmla, data = data, weights = w2_sd)
    }
    
    # Bayesian regression (Gaussian family)
    model_list$Bayes <- bayesglm(fmla, data = data)
    
    # Random Forest
    set.seed(42)
    model_list$RF <- randomForest(fmla, data = data, ntree = 500)
    
    # PLS with LOOCV
    pls_loo <- plsr(fmla, data = data, scale = TRUE, validation = "LOO")
    best_loo <- selectNcomp(pls_loo, method = "onesigma", plot = FALSE)
    attr(pls_loo, "best_ncomp") <- best_loo
    model_list$PLS_LOO <- pls_loo
    
    # PLS with k-fold CV (10 segments)
    pls_k <- plsr(fmla, data = data, scale = TRUE,
                  validation = "CV", segments = 10)
    best_k <- selectNcomp(pls_k, method = "onesigma", plot = FALSE)
    attr(pls_k, "best_ncomp") <- best_k
    model_list$PLS_kfold <- pls_k
    
    
    ############################################################
    # 7.2 DATAFRAMES FOR THIS ELEMENT
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
    # 7.3 PER-MODEL ANALYSIS LOOP
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
      
      # Predictions for full dataset
      pred <- tryCatch({
        if (inherits(m, "mvr")) {
          as.numeric(predict(m, ncomp = attr(m, "best_ncomp"), newdata = data))
        } else {
          predict(m)
        }
      }, error = function(e) NA_real_)
      
      if (any(is.na(pred))) next
      
      resid_vals <- y - pred
      
      # Correct R²
      SS_res <- sum((y - pred)^2)
      SS_tot <- sum((y - mean(y))^2)
      R2     <- 1 - SS_res / SS_tot
      
      # RMSE
      RMSE   <- sqrt(mean((y - pred)^2))
      
      # RMSEP from 10-fold CV
      cv10 <- tryCatch(
        cv_rmse_for_model(
          mname, fmla, data, y_col, var_col, sd_col,
          best_ncomp = if (inherits(m, "mvr")) attr(m, "best_ncomp") else NA
        ),
        error = function(e) NA_real_
      )
      RMSEP_val <- cv10
      
      # AIC / BIC only for LM/GLM/Bayes
      AIC_val <- if (inherits(m, c("lm", "glm", "bayesglm")))
        tryCatch(AIC(m), error = function(e) NA_real_) else NA_real_
      
      BIC_val <- if (inherits(m, c("lm", "glm", "bayesglm")))
        tryCatch(BIC(m), error = function(e) NA_real_) else NA_real_
      
      # Save performance into per-element table
      ranking_df <- add_row(
        ranking_df,
        Model = mname,
        R2    = R2,
        RMSE  = RMSE,
        RMSEP = RMSEP_val,
        AIC   = AIC_val,
        BIC   = BIC_val
      )
      
      # Also append to global summary
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
      
      ##########################################################
      # Diagnostics: Shapiro, Breusch–Pagan, CV info
      ##########################################################
      
      sh_W <- sh_p <- NA_real_
      if (length(resid_vals) >= 3 && length(resid_vals) <= 5000) {
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
          bp_stat <- unname(bp$statistic)
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
        CV10_RMSE         = cv10,
        CV10_expl         = expl_cv10,
        AIC               = AIC_val,
        AIC_expl          = "AIC is reported but not used for ranking.",
        BIC               = BIC_val,
        BIC_expl          = "BIC is reported but not used for ranking."
      )
      
      ##########################################################
      # PLOTS
      ##########################################################
      
      # 1) Predicted vs Observed
      df_pred <- data.frame(
        Observed  = y,
        Predicted = pred,
        Site      = data$Site
      )
      
      p_po <- ggplot(df_pred, aes(Observed, Predicted, color = Site)) +
        geom_point(size = 3, alpha = 0.85) +
        ggsci::scale_color_jco() +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() +
        ggtitle(paste(el, mname, "- Predicted vs Observed"))
      
      ggsave(
        filename = file.path(element_dir, paste0(el, "_", mname, "_PredVsObs.pdf")),
        plot     = p_po,
        width    = 8,
        height   = 8
      )
      
      # 2) Residual distribution
      p_res <- ggplot(data.frame(resid = resid_vals), aes(resid)) +
        geom_histogram(bins = 30, fill = "grey80") +
        geom_density(color = "red") +
        theme_bw() +
        ggtitle(paste(el, mname, "- Residuals"))
      
      ggsave(
        filename = file.path(element_dir, paste0(el, "_", mname, "_Residuals.pdf")),
        plot     = p_res,
        width    = 8,
        height   = 8
      )
      
      # 3) Influence plots for classical LM models (not Bayesian)
      if (inherits(m, "lm") && !inherits(m, "bayesglm")) {
        
        cooks <- cooks.distance(m)
        h     <- hatvalues(m)
        rstd  <- rstandard(m)
        
        inf_df <- data.frame(
          obs      = seq_along(cooks),
          Cook     = cooks,
          Leverage = h,
          StdResid = rstd
        )
        
        p_cook <- ggplot(inf_df, aes(obs, Cook)) +
          geom_col() +
          theme_bw() +
          ggtitle(paste(el, mname, "- Cook's Distance"))
        
        p_lev <- ggplot(inf_df, aes(Leverage, StdResid)) +
          geom_point() +
          theme_bw() +
          ggtitle(paste(el, mname, "- Leverage Plot"))
        
        pdf(file.path(element_dir, paste0(el, "_", mname, "_Influence.pdf")),
            width = 8, height = 8)
        print(p_cook)
        print(p_lev)
        dev.off()
      }
      
      # 4) PLS RMSEP vs components
      if (inherits(m, "mvr")) {
        pdf(file.path(element_dir, paste0(el, "_", mname, "_PLS_RMSEP.pdf")),
            width = 8, height = 8)
        validationplot(
          m,
          val.type = "RMSEP",
          main = paste(el, mname, "- RMSEP vs Components")
        )
        dev.off()
      }
      
    } # end model loop for element
    
    ############################################################
    # 7.4 RANK MODELS PER ELEMENT — R² → RMSEP → RMSE
    ############################################################
    
    if (nrow(ranking_df) > 0) {
      
      # Rank: 1) highest R², 2) lowest RMSEP, 3) lowest RMSE
      ranking_df <- ranking_df %>%
        arrange(
          dplyr::desc(R2),   # R² high → low
          RMSEP,             # RMSEP low → high
          RMSE               # RMSE low → high
        ) %>%
        mutate(Rank = dplyr::row_number())
      
      # Detailed ranking + explanations
      ranking_expl_df <- ranking_df %>%
        mutate(
          Rank_explanation  = expl_rank,
          R2_explanation    = expl_r2,
          RMSE_explanation  = expl_rmse,
          RMSEP_explanation = expl_rmsep
        )
      
      # Save per-element ranking CSV
      write.csv(
        ranking_expl_df,
        file.path(element_dir, paste0(el, "_Model_Ranking.csv")),
        row.names = FALSE
      )
      
      # Keep slimmer version in memory if needed
      ranking_df <- ranking_df %>%
        dplyr::select(Model, R2, RMSE, RMSEP, AIC, BIC, Rank)
    }
    
    
    ############################################################
    # 7.5 DIAGNOSTICS CSV
    ############################################################
    
    if (nrow(diag_df) > 0) {
      write.csv(
        diag_df,
        file.path(element_dir, paste0(el, "_Diagnostics_Summary.csv")),
        row.names = FALSE
      )
    }
    
    
    ############################################################
    # 7.6 README.txt PER ELEMENT
    ############################################################
    
    readme_path <- file.path(element_dir, "README.txt")
    
    readme_txt <- c(
      paste("Element:", el),
      "",
      "FILES PRODUCED:",
      paste0("- ", el, "_Model_Ranking.csv : Summary of all regression models for this element,"),
      "  including R², RMSEP, RMSE, AIC and BIC, ranked as:",
      "    1) R² (highest → lowest),",
      "    2) RMSEP (lowest → highest),",
      "    3) RMSE (lowest → highest).",
      paste0("- ", el, "_Diagnostics_Summary.csv :"),
      "  Shapiro–Wilk normality test for residuals,",
      "  Breusch–Pagan heteroscedasticity test,",
      "  and 10-fold cross-validation RMSE values for each model.",
      "- *_PredVsObs.pdf : Predicted vs observed plots for each model.",
      "- *_Residuals.pdf : Residual distribution plots for each model.",
      "- *_Influence.pdf : Cook's distance & leverage plots (for LM models only).",
      "- *_PLS_RMSEP.pdf : RMSEP vs number of components for PLS models.",
      "",
      "RANKING LOGIC (PER ELEMENT):",
      "1. Rank 1 always corresponds to the model with the highest R².",
      "2. When models share similar R², those with lower RMSEP are ranked higher.",
      "3. When R² and RMSEP are similar, models with lower RMSE are ranked higher.",
      "4. AIC and BIC are reported for context only and are not used in ranking.",
      "",
      "DIAGNOSTICS:",
      "- Shapiro–Wilk: tests residual normality; p > 0.05 is usually desirable.",
      "- Breusch–Pagan: tests for heteroscedasticity; p < 0.05 suggests",
      "  non-constant error variance.",
      "- CV10_RMSE: 10-fold cross-validation RMSE; lower values indicate",
      "  better predictive performance.",
      "",
      "NOTES:",
      "- AIC/BIC may be NA where they are not defined (e.g., RF, PLS).",
      "- Influence plots are only computed for classical lm() models.",
      "- Site colours in PredVsObs plots use the 'jco' palette with Site",
      "  ordered as: BI10, HER42PB, KER1, KER3, PB1.",
      "",
      "Generated automatically by run_full_regressions()."
    )
    
    writeLines(readme_txt, readme_path)
    
  }  # end element loop
  
  
  ############################################################
  # 8. GLOBAL SUMMARY ACROSS ALL ELEMENTS — R² → RMSEP → RMSE
  ############################################################
  
  if (nrow(global_summary) > 0) {
    
    global_ranked <- global_summary %>%
      arrange(
        Element,
        dplyr::desc(R2),   # R² high → low within each element
        RMSEP,             # RMSEP low → high
        RMSE               # RMSE low → high
      ) %>%
      group_by(Element) %>%
      mutate(Rank = dplyr::row_number()) %>%
      ungroup()
    
    # Save all-model summary across elements
    write.csv(
      global_ranked,
      file.path(save_dir, summary_all_csv),
      row.names = FALSE
    )
    
    # Best model per element (Rank 1)
    best_models <- global_ranked %>%
      dplyr::filter(Rank == 1) %>%
      dplyr::select(Element, Model, R2, RMSE, RMSEP, AIC, BIC)
    
    write.csv(
      best_models,
      file.path(save_dir, best_models_csv),
      row.names = FALSE
    )
  }
  
  
  ############################################################
  # 9. DONE
  ############################################################
  
  message("=== COMPLETED SUCCESSFULLY ===")
  message("Results saved to directory: ", save_dir)
}

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Run it & check it works! --------------------------------------------------------
elements <- c("Ti","Ca","Mn","Fe","Sr","Zr")
run_full_regressions(elements, ACE_dataset)
read.csv("/Users/sjro/Desktop/Regression/Ti/Ti_Model_Ranking.csv") # chank Rank order 
# -------------------------------------------------------------------------


