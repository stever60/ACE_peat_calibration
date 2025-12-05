# Regression Comparison function

# Log_inc

# Regression Comparison function description ----------------------------------------

# This function performs a complete multi-model calibration and diagnostic workflow 
# for each chemical element in a dataset, fitting Ordinary Least Squares (OLS), 
# variance-weighted Weighted Least Squares (WLS), 
# SD-weighted OLS, two-step weighted WLS, Bayesian regression (MAP), 
# Random Forest (RF) and Partial Least Squares (PLS) regression. 

# It runs simple and weighted models and for every model and each element (e.g., 
# inputtted as Ti for Ti_XRF and Ti_ICP for TI_ICP-MS - it then generates full 
# diagnostic plots, and undertakes residual and influence analyses, normality and 
# heteroscedasticity  tests, bootstrap confidence intervals, 10-fold cross-validation 
# performance, and produces multi-model comparison plots. 
# 
# This function performs eight calibration models (OLS, WLS, weighted OLS, 
# weighted WLS, Bayesian regression (Bays), Random Forest (RF), PLS-LOO, PLS-10-fold CV) for 
# each element and generates a complete diagnostic suite: residual analysis, 
# normality testing, heteroscedasticity tests, influence diagnostics, bootstrap 
# CIs, predicted vs observed plots coloured by Site using ggsci::color_palette("jco"). 
# Variable importance (RF), PLS validation curves, individual diagnostic PDFs, 
# a combined summary for all models/elements, a .csv file per element of all model 
# metrics (R², RMSE, RMSEP, AIC, BIC) are also produced.

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



# Set up ------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()

# -------------------------------------------------------------------------
# Libraries ---------------------------------------------------------------
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

# Define elements to use -------------------------------------------------------

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

# Import ACE matched log (element/inc) file  ----------------------------
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

names(ACE_dataset)

c("Ti", "Ti_ICP") %in% names(ACE_dataset)
c("Ca", "Ca_ICP") %in% names(ACE_dataset)
c("Fe", "Fe_ICP") %in% names(ACE_dataset)


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Check dataset is loaded & working before runnnig function -------------
summary(ACE_dataset$Ti)
var(ACE_dataset$Ti)
summary(ACE_dataset$Ti_ICP)
var(ACE_dataset$Ti_ICP)
cor(ACE_dataset$Ti, ACE_dataset$Ti_ICP)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Regression comparison function ---------------------------------------------------------------------

run_full_regressions <- function(
    elements,
    data,
    var_suffix      = "_var",
    sd_suffix       = "_ICP_sd",
    save_dir        = "/Users/sjro/Desktop/Regression_outputs",
    summary_all_csv = "AllElements_ModelSummary.csv",
    best_models_csv = "BestModels_PerElement.csv"
) {
  
  ############################################################
  # 1. SETUP & PACKAGE LOADING
  ############################################################
  
  library(ggplot2)
  library(dplyr)
  library(randomForest)
  library(pls)
  library(arm)
  library(ggsci)
  library(lmtest)
  library(progress)
  
  save_dir <- normalizePath(save_dir, mustWork = FALSE)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  writeLines("TEST", file.path(save_dir, "TEST_WRITE.txt"))
  
  log_file <- file.path(save_dir, "run_log.txt")
  writeLines("=== Regression Run Log ===", log_file)
  
  if (!"Site" %in% names(data))
    stop("Dataset must contain a 'Site' column.")
  
  data$Site <- factor(
    data$Site,
    levels = c("BI10", "HER42PB", "KER1", "KER3", "PB1")
  )
  
  ############################################################
  # 2. HELPER FUNCTIONS & EXPLANATIONS
  ############################################################
  
  # 2.1 Dataset integrity check
  check_dataset_integrity <- function(data, el, y_col, x_col, var_col, sd_col) {
    errors <- c()
    
    required <- c(y_col, x_col)
    for (col in required)
      if (!col %in% names(data))
        errors <- c(errors, paste("Missing:", col))
    
    if (!is.numeric(data[[y_col]])) errors <- c(errors, paste(y_col, "not numeric"))
    if (!is.numeric(data[[x_col]])) errors <- c(errors, paste(x_col, "not numeric"))
    
    if (any(is.na(data[[y_col]]))) errors <- c(errors, paste("NA in", y_col))
    if (any(is.na(data[[x_col]]))) errors <- c(errors, paste("NA in", x_col))
    
    if (var(data[[x_col]], na.rm = TRUE) == 0)
      errors <- c(errors, paste("Predictor", x_col, "has zero variance"))
    
    if (sd_col %in% names(data))
      if (any(data[[sd_col]] <= 0, na.rm = TRUE))
        errors <- c(errors, paste(sd_col, "contains non-positive values"))
    
    if (var_col %in% names(data))
      if (any(data[[var_col]] <= 0, na.rm = TRUE))
        errors <- c(errors, paste(var_col, "contains non-positive values"))
    
    errors
  }
  
  # 2.2 Cross-validation RMSE (10-fold)
  cv_rmse_for_model <- function(
    mname, formula, data, y_col, var_col, sd_col,
    best_ncomp = NA, k = 10
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
          
        } else if (mname == "WLS_var" && var_col %in% names(train)) {
          lm(formula, data = train, weights = 1/train[[var_col]])
          
        } else if (mname == "OLS_weighted" && sd_col %in% names(train)) {
          lm(formula, data = train, weights = 1/(train[[sd_col]]^2))
          
        } else if (mname == "WLS_weighted" && sd_col %in% names(train)) {
          init <- lm(formula, data = train, weights = 1/(train[[sd_col]]^2))
          rm   <- lm(abs(init$residuals) ~ init$fitted.values)
          w2   <- 1/(rm$fitted.values^2)
          lm(formula, data = train, weights = w2)
          
        } else if (mname == "Bayes") {
          bayesglm(formula, data = train)
          
        } else if (mname == "RF") {
          randomForest(formula, data = train, ntree = 300)
          
        } else if (mname %in% c("PLS_LOO", "PLS_kfold")) {
          plsr(formula, data = train, scale = TRUE, validation = "none")
          
        } else {
          NULL
        }
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
    sqrt(cv_sse/cv_n)
  }
  
  # 2.3 Explanation text (for CSVs)
  
  expl_shapiro <- "Shapiro–Wilk test of residual normality: p > 0.05 suggests residuals are approximately normal."
  expl_bp      <- "Breusch–Pagan test for heteroscedasticity: p < 0.05 suggests non-constant residual variance."
  expl_cv10    <- "10-fold cross-validation RMSE: lower values indicate better out-of-sample predictive performance."
  expl_aic     <- "AIC (Akaike Information Criterion): lower values indicate better trade-off between fit and complexity."
  expl_bic     <- "BIC (Bayesian Information Criterion): like AIC but with stronger penalty for complexity; lower is better."
  
  expl_rank    <- paste(
    "Ranking logic: (1) All PLS models are ranked above non-PLS models;",
    "(2) among non-PLS models, if R2, RMSE and RMSEP are equal,",
    "WLS models are preferred over OLS_weighted, which are preferred over OLS;",
    "(3) otherwise, within each group, models are ordered by higher R2,",
    "lower RMSE, lower RMSEP, then lower AIC/BIC where available.",
    sep = " "
  )
  
  expl_r2 <- "R2 is the fraction of variance in the ICP data explained by the model (1 − SS_res / SS_tot); higher is better."
  expl_rmse <- "RMSE is the root mean squared error of calibration; lower values indicate better fit on the observed data."
  expl_rmsep <- "RMSEP is RMSE estimated by cross-validation (10-fold here when available); lower values indicate better predictive performance."
  expl_aicbic <- "AIC and BIC penalise model complexity; lower values indicate a better balance between fit and complexity. NA means these criteria are not available for that model type."
  
  ############################################################
  # 3. GLOBAL SUMMARY INITIALISATION
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
  # 4. MAIN LOOP OVER ELEMENTS (with progress bar)
  ############################################################
  
  pb_elements <- progress_bar$new(
    format = "Elements [:bar] :current/:total (:percent) eta::eta",
    total  = length(elements),
    clear  = FALSE,
    width  = 80
  )
  
  
  ############################################################
  # 4. PRIORITY HELPER FOR RANKING
  ############################################################
  
  # Smaller number = higher priority
  ranking_priority <- function(model_name) {
    if (grepl("^PLS", model_name))        return(1)  # All PLS best
    if (grepl("^WLS", model_name))        return(2)  # WLS above OLS
    if (grepl("^OLS_weighted", model_name)) return(3)
    if (model_name == "OLS")              return(4)
    return(5)  # everything else last
  }
  
  ############################################################
  # 5. MAIN LOOP OVER ELEMENTS (with progress bar)
  ############################################################
  
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
      write(paste("Skipping", el, ":", paste(errors, collapse = "; ")),
            file = log_file, append = TRUE)
      next
    }
    
    fmla <- as.formula(paste0(y_col, " ~ ", x_col))
    y <- data[[y_col]]
    
    ############################################################
    # 5.1 FIT MODELS (NO WEIGHTED PLS)
    ############################################################
    
    model_list <- list()
    
    # OLS
    model_list$OLS <- lm(fmla, data = data)
    
    # WLS using variance column if present
    if (var_col %in% names(data)) {
      model_list$WLS_var <- lm(fmla, data = data, weights = 1/data[[var_col]])
    }
    
    # Weighted OLS and 2-step WLS using SD if present
    if (sd_col %in% names(data)) {
      w_sd <- 1/(data[[sd_col]]^2)
      model_list$OLS_weighted <- lm(fmla, data = data, weights = w_sd)
      
      init_sd <- lm(fmla, data = data, weights = w_sd)
      rm_sd   <- lm(abs(init_sd$residuals) ~ init_sd$fitted.values)
      w2_sd   <- 1/(rm_sd$fitted.values^2)
      model_list$WLS_weighted <- lm(fmla, data = data, weights = w2_sd)
    }
    
    # Bayesian regression
    model_list$Bayes <- bayesglm(fmla, data = data)
    
    # Random Forest
    set.seed(42)
    model_list$RF <- randomForest(fmla, data = data, ntree = 500)
    
    # Unweighted PLS (LOO and k-fold)
    pls_loo <- plsr(fmla, data = data, scale = TRUE, validation = "LOO")
    best_loo <- selectNcomp(pls_loo, method = "onesigma", plot = FALSE)
    attr(pls_loo, "best_ncomp") <- best_loo
    model_list$PLS_LOO <- pls_loo
    
    pls_k <- plsr(fmla, data = data, scale = TRUE,
                  validation = "CV", segments = 10)
    best_k <- selectNcomp(pls_k, method = "onesigma", plot = FALSE)
    attr(pls_k, "best_ncomp") <- best_k
    model_list$PLS_kfold <- pls_k
    
    ############################################################
    # 5.2 PER-ELEMENT DATAFRAMES (add_row-safe)
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
    # 5.3 MODEL ANALYSIS LOOP (with progress bar)
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
      
      # Predictions
      pred <- tryCatch({
        if (inherits(m, "mvr"))
          as.numeric(predict(m, ncomp = attr(m, "best_ncomp"), newdata = data))
        else
          predict(m)
      }, error = function(e) NA_real_)
      
      if (any(is.na(pred))) next
      
      resid_vals <- y - pred
      
      # Correct R²
      SS_res <- sum((y - pred)^2)
      SS_tot <- sum((y - mean(y))^2)
      R2     <- 1 - SS_res/SS_tot
      
      RMSE <- sqrt(mean((y - pred)^2))
      
      # CV RMSE
      cv10 <- tryCatch(
        cv_rmse_for_model(
          mname, fmla, data, y_col, var_col, sd_col,
          best_ncomp = if (inherits(m, "mvr")) attr(m, "best_ncomp") else NA
        ),
        error = function(e) NA_real_
      )
      
      RMSEP_val <- if (!is.na(cv10)) cv10 else RMSE
      
      # AIC / BIC
      AIC_val <- if (inherits(m, c("lm", "glm", "bayesglm")))
        tryCatch(AIC(m), error = function(e) NA) else NA
      
      BIC_val <- if (inherits(m, c("lm", "glm", "bayesglm")))
        tryCatch(BIC(m), error = function(e) NA) else NA
      
      # Per-element ranking data
      ranking_df <- add_row(
        ranking_df,
        Model = mname, R2 = R2, RMSE = RMSE,
        RMSEP = RMSEP_val, AIC = AIC_val, BIC = BIC_val
      )
      
      # Global summary data
      global_summary <- add_row(
        global_summary,
        Element = el, Model = mname, R2 = R2, RMSE = RMSE,
        RMSEP = RMSEP_val, AIC = AIC_val, BIC = BIC_val
      )
      
      # Diagnostics
      sh_W <- sh_p <- NA
      if (length(resid_vals) >= 3 && length(resid_vals) <= 5000) {
        sh <- tryCatch(shapiro.test(resid_vals), error = function(e) NULL)
        if (!is.null(sh)) {
          sh_W <- sh$statistic
          sh_p <- sh$p.value
        }
      }
      
      bp_stat <- bp_p <- NA
      if (inherits(m, c("lm", "bayesglm"))) {
        bp <- tryCatch(bptest(m), error = function(e) NULL)
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
        CV10_RMSE         = cv10,
        CV10_expl         = expl_cv10,
        AIC               = AIC_val,
        AIC_expl          = expl_aic,
        BIC               = BIC_val,
        BIC_expl          = expl_bic
      )
      
      # Plots: Pred vs Obs
      df_pred <- data.frame(
        Observed = y,
        Predicted = pred,
        Site = data$Site
      )
      
      p_po <- ggplot(df_pred, aes(Observed, Predicted, color = Site)) +
        geom_point(size = 3, alpha = 0.85) +
        ggsci::scale_color_jco() +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() +
        ggtitle(paste(el, mname, "- Predicted vs Observed"))
      
      ggsave(file.path(element_dir, paste0(el, "_", mname, "_PredVsObs.pdf")),
             p_po, width = 8, height = 8)
      
      # Residuals
      p_res <- ggplot(data.frame(resid = resid_vals), aes(resid)) +
        geom_histogram(bins = 30, fill = "grey80") +
        geom_density(color = "red") +
        theme_bw() +
        ggtitle(paste(el, mname, "- Residuals"))
      
      ggsave(file.path(element_dir, paste0(el, "_", mname, "_Residuals.pdf")),
             p_res, width = 8, height = 8)
      
      # Influence plots for classical LM models (not Bayes)
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
      
      # PLS RMSEP plot
      if (inherits(m, "mvr")) {
        pdf(file.path(element_dir, paste0(el, "_", mname, "_PLS_RMSEP.pdf")),
            width = 8, height = 8)
        validationplot(m, val.type = "RMSEP",
                       main = paste(el, mname, "- RMSEP"))
        dev.off()
      }
      
    } # end model loop
    
    
    ############################################################
    # 5.4 RANK MODELS PER ELEMENT (priority + metrics)
    ############################################################
    
    if (nrow(ranking_df) > 0) {
      
      ranking_df <- ranking_df %>%
        mutate(Priority = sapply(Model, ranking_priority)) %>%
        arrange(
          Priority,      # 1 = PLS, 2 = WLS, 3 = OLS_weighted, 4 = OLS, 5 = others
          dplyr::desc(R2),
          RMSE,
          RMSEP,
          is.na(AIC), AIC,
          is.na(BIC), BIC
        ) %>%
        mutate(Rank = dplyr::row_number())
      
      # Add explanatory text columns
      ranking_df <- ranking_df %>%
        mutate(
          Rank_explanation      = expl_rank,
          R2_explanation        = expl_r2,
          RMSE_explanation      = expl_rmse,
          RMSEP_explanation     = expl_rmsep,
          AIC_BIC_explanation   = expl_aicbic
        )
      
      ranking_df <- dplyr::select(
        ranking_df,
        Rank, Model, R2, RMSE, RMSEP, AIC, BIC,
        Rank_explanation, R2_explanation, RMSE_explanation,
        RMSEP_explanation, AIC_BIC_explanation
      )
      
      write.csv(
        ranking_df,
        file.path(element_dir, paste0(el, "_Model_Ranking.csv")),
        row.names = FALSE
      )
    }
    
    ############################################################
    # 5.5 DIAGNOSTICS CSV (already includes explanations)
    ############################################################
    
    if (nrow(diag_df) > 0) {
      write.csv(
        diag_df,
        file.path(element_dir, paste0(el, "_Diagnostics_Summary.csv")),
        row.names = FALSE
      )
    }
    
    ############################################################
    # 5.6 README.txt FOR ELEMENT
    ############################################################
    
    readme_path <- file.path(element_dir, "README.txt")
    
    readme_txt <- c(
      paste("Element:", el),
      "",
      "FILES PRODUCED:",
      paste0("- ", el, "_Model_Ranking.csv : Ranked list of all regression models for this element."),
      paste0("- ", el, "_Diagnostics_Summary.csv : Normality, heteroscedasticity, cross-validation diagnostics per model."),
      "- *_PredVsObs.pdf : Predicted vs observed plots by model.",
      "- *_Residuals.pdf : Residual distribution plots by model.",
      "- *_Influence.pdf : Cook's distance & leverage plots (LM models only).",
      "- *_PLS_RMSEP.pdf : PLS RMSEP vs components (PLS models only).",
      "",
      "RANKING LOGIC:",
      "1. All PLS models (PLS_LOO, PLS_kfold) are ranked above OLS and WLS models,",
      "   because PLS is better suited to multicollinearity and can cope better",
      "   with heteroscedastic behaviour in geochemical data.",
      "2. Among non-PLS models, if R2, RMSE and RMSEP are equal,",
      "   WLS models are preferred over OLS_weighted, which are preferred over OLS,",
      "   because they explicitly account for measurement error/variance.",
      "3. Within each priority group, models are ranked by:",
      "   - Highest R2 (better fit),",
      "   - Lowest RMSE (better calibration error),",
      "   - Lowest RMSEP (better cross-validated prediction error),",
      "   - Lowest AIC (if available),",
      "   - Lowest BIC (if available).",
      "",
      "DIAGNOSTICS EXPLANATIONS:",
      "- Shapiro–Wilk: Tests residual normality; p > 0.05 suggests residuals are approximately normal.",
      "- Breusch–Pagan: Tests heteroscedasticity; p < 0.05 suggests non-constant residual variance.",
      "- CV10 RMSE: 10-fold cross-validation error; lower values imply better predictive power.",
      "- AIC/BIC: Penalised likelihood criteria; lower values indicate better trade-off between fit and complexity.",
      "",
      "NOTES:",
      "- Missing AIC/BIC (NA) means those criteria are not defined for that model type (e.g. random forest, PLS).",
      "- Residual-based diagnostics and influence plots are only defined for LM-class models.",
      "- All Site colours in PredVsObs plots follow ggpubr/ggsci 'jco' palette with Site-level ordering BI10, HER42PB, KER1, KER3, PB1.",
      "",
      "Generated automatically by run_full_regressions()."
    )
    
    writeLines(readme_txt, readme_path)
    
  }  # end element loop
  
  ############################################################
  # 6. GLOBAL SUMMARY ACROSS ALL ELEMENTS
  ############################################################
  
  if (nrow(global_summary) > 0) {
    
    # Apply same priority ranking logic across elements
    global_ranked <- global_summary %>%
      mutate(Priority = sapply(Model, ranking_priority)) %>%
      group_by(Element) %>%
      arrange(
        Element,
        Priority,
        dplyr::desc(R2),
        RMSE,
        RMSEP,
        is.na(AIC), AIC,
        is.na(BIC), BIC
      ) %>%
      mutate(Rank = dplyr::row_number()) %>%
      ungroup() %>%
      arrange(Element, Rank)
    
    write.csv(
      global_ranked,
      file.path(save_dir, summary_all_csv),
      row.names = FALSE
    )
    
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
  # 7. DONE
  ############################################################
  
  message("=== COMPLETED SUCCESSFULLY ===")
  message("Results saved to directory: ", save_dir)
}


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Run it  --------------------------------------------------------
elements <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
run_full_regressions(elements, ACE_dataset)
# -------------------------------------------------------------------------

# ***************************** END *************************************** 












# Old versions
# Regression comparison function - eaerlier version - without R2 (true) ---------------------------------------------------------------------

run_full_regressions <- function(
    elements,
    data,
    var_suffix = "_var",
    sd_suffix  = "_ICP_sd",
    save_dir   = "/Users/sjro/Desktop/Regression_old",
    summary_csv = "ModelSummary_AllElements.csv"  # reserved if you want later
) {
  
  # ==== PACKAGES (NO PARALLEL, NO PROGRESSR, NO SEE/PERFORMANCE) ====
  library(ggplot2)
  library(dplyr)
  library(randomForest)
  library(pls)
  library(arm)
  library(ggsci)
  
  # =====================================================
  # CREATE SAVE DIRECTORY + TEST WRITE
  # =====================================================
  message("Requested save_dir: ", save_dir)
  save_dir <- normalizePath(save_dir, mustWork = FALSE)
  message("Normalized save_dir: ", save_dir)
  
  if (!dir.exists(save_dir)) {
    ok <- dir.create(save_dir, recursive = TRUE)
    if (!ok) stop("Could not create directory: ", save_dir)
    message("Created directory: ", save_dir)
  } else {
    message("Directory already exists: ", save_dir)
  }
  
  # Try writing a simple test file
  test_file <- file.path(save_dir, "TEST_WRITE.txt")
  writeLines("If you can see this file, the path is correct.", test_file)
  message("Wrote test file: ", test_file)
  
  # Log file
  log_file <- file.path(save_dir, "run_log.txt")
  writeLines("=== Regression Run Log ===", log_file)
  message("Initialized log file: ", log_file)
  
  # =====================================================
  # ENSURE Site FACTOR ORDER
  # =====================================================
  if ("Site" %in% names(data)) {
    data$Site <- factor(
      data$Site,
      levels = c("BI10", "HER42PB", "KER1", "KER3", "PB1")
    )
  } else {
    stop("Dataset must contain a 'Site' column.")
  }
  
  # =====================================================
  # DATASET INTEGRITY CHECK
  # =====================================================
  check_dataset_integrity <- function(data, el, y_col, x_col, var_col, sd_col) {
    errors <- c()
    
    needed <- c(y_col, x_col)
    for (col in needed) {
      if (!col %in% names(data))
        errors <- c(errors, paste("Missing column:", col))
    }
    if (length(errors)) return(errors)
    
    if (!is.numeric(data[[y_col]]))
      errors <- c(errors, paste(y_col, "is not numeric"))
    if (!is.numeric(data[[x_col]]))
      errors <- c(errors, paste(x_col, "is not numeric"))
    
    if (any(is.na(data[[y_col]])))
      errors <- c(errors, paste("NA values in", y_col))
    if (any(is.na(data[[x_col]])))
      errors <- c(errors, paste("NA values in", x_col))
    
    if (var(data[[x_col]], na.rm = TRUE) == 0)
      errors <- c(errors, paste(x_col, "has zero variance"))
    
    if (sd_col %in% names(data)) {
      if (any(data[[sd_col]] <= 0, na.rm = TRUE))
        errors <- c(errors, paste(sd_col, "contains non-positive values"))
    }
    if (var_col %in% names(data)) {
      if (any(data[[var_col]] <= 0, na.rm = TRUE))
        errors <- c(errors, paste(var_col, "contains non-positive values"))
    }
    
    errors
  }
  
  # =====================================================
  # MAIN SEQUENTIAL LOOP
  # =====================================================
  for (el in elements) {
    
    message("\n==========================")
    message("Processing element: ", el)
    message("==========================")
    
    element_log <- paste0("\n--- Element: ", el, " ---\n")
    
    y_col  <- paste0(el, "_ICP")
    x_col  <- el
    var_col <- paste0(el, var_suffix)
    sd_col  <- paste0(el, sd_suffix)
    
    # Element-specific directory
    element_dir <- file.path(save_dir, el)
    dir.create(element_dir, recursive = TRUE, showWarnings = FALSE)
    message("Element directory: ", element_dir)
    
    # ---- Integrity check ----
    errors <- check_dataset_integrity(data, el, y_col, x_col, var_col, sd_col)
    
    if (length(errors) > 0) {
      element_log <- paste0(
        element_log,
        "Dataset integrity FAILED:\n",
        paste(" -", errors, collapse = "\n"),
        "\nSkipping element.\n"
      )
      write(element_log, file = log_file, append = TRUE)
      next
    }
    
    element_log <- paste0(element_log, "Dataset integrity OK.\n")
    
    # =====================================================
    # FIT MODELS
    # =====================================================
    fmla <- as.formula(paste0(y_col, " ~ ", x_col))
    y <- data[[y_col]]
    
    model_list <- list()
    
    # OLS
    model_list$OLS <- lm(fmla, data = data)
    
    # WLS on variance
    if (var_col %in% names(data)) {
      model_list$WLS_var <- lm(fmla, data = data, weights = 1 / data[[var_col]])
    }
    
    # Weighted OLS & 2-step WLS on SD
    if (sd_col %in% names(data)) {
      w_sd <- 1 / (data[[sd_col]]^2)
      model_list$OLS_weighted <- lm(fmla, data = data, weights = w_sd)
      
      init <- lm(fmla, data = data, weights = w_sd)
      rm   <- lm(abs(init$residuals) ~ init$fitted.values)
      w2   <- 1 / (rm$fitted.values^2)
      model_list$WLS_weighted <- lm(fmla, data = data, weights = w2)
    }
    
    # Bayesian
    model_list$Bayes <- bayesglm(fmla, data = data)
    
    # Random forest
    set.seed(42)
    model_list$RF <- randomForest(fmla, data = data, ntree = 500)
    
    # PLS LOO
    pls_loo <- plsr(fmla, data = data, scale = TRUE, validation = "LOO")
    best_loo <- selectNcomp(pls_loo, method = "onesigma", plot = FALSE)
    attr(pls_loo, "best_ncomp") <- best_loo
    model_list$PLS_LOO <- pls_loo
    
    # PLS k-fold
    pls_k <- plsr(fmla, data = data, scale = TRUE,
                  validation = "CV", segments = 10)
    best_k <- selectNcomp(pls_k, method = "onesigma", plot = FALSE)
    attr(pls_k, "best_ncomp") <- best_k
    model_list$PLS_kfold <- pls_k
    
    # =====================================================
    # DIAGNOSTICS + RANKING
    # =====================================================
    ranking_df <- dplyr::tibble()
    
    for (mname in names(model_list)) {
      
      m <- model_list[[mname]]
      message("  - Model: ", mname)
      
      # Safe predictions
      pred <- tryCatch({
        if (inherits(m, "mvr")) {
          ncomp <- attr(m, "best_ncomp")
          as.numeric(predict(m, ncomp = ncomp, newdata = data))
        } else {
          predict(m)
        }
      }, error = function(e) {
        message("    Prediction failed for ", mname, ": ", e$message)
        NA
      })
      
      if (!is.numeric(pred) || length(pred) != nrow(data) || any(is.na(pred))) {
        element_log <- paste0(element_log,
                              "Prediction invalid for ", mname, " (skipped).\n")
        next
      }
      
      resid_vals <- y - pred
      
      # -------- Residual distribution plot --------
      if (length(resid_vals) > 1 && !all(is.na(resid_vals))) {
        res_df <- data.frame(resid = resid_vals)
        p_res <- ggplot(res_df, aes(x = resid)) +
          geom_histogram(fill = "grey80", bins = 30) +
          geom_density(color = "red") +
          theme_bw() +
          ggtitle(paste(el, mname, "Residual Distribution"))
        
        ggsave(file.path(element_dir, paste0(el, "_", mname, "_Residuals.pdf")),
               p_res, width = 8, height = 8)
      } else {
        element_log <- paste0(element_log,
                              "Residuals invalid for ", mname, " (plot skipped).\n")
      }
      
      # -------- Predicted vs Observed --------
      df_pred <- data.frame(Observed = y, Predicted = pred, Site = data$Site)
      
      p_po <- ggplot(df_pred, aes(x = Observed, y = Predicted, color = Site)) +
        geom_point(size = 3, alpha = 0.8) +
        ggsci::scale_color_jco() +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw() +
        ggtitle(paste(el, mname, "Predicted vs Observed"))
      
      ggsave(file.path(element_dir, paste0(el, "_", mname, "_PredVsObs.pdf")),
             p_po, width = 8, height = 8)
      
      # -------- PLS RMSEP validation --------
      if (inherits(m, "mvr")) {
        pdf(file.path(element_dir, paste0(el, "_", mname, "_PLS_RMSEP.pdf")),
            width = 8, height = 8)
        validationplot(m, val.type = "RMSEP",
                       main = paste(el, mname, "RMSEP vs Components"))
        dev.off()
      }
      
      # -------- metrics for ranking --------
      R2    <- cor(y, pred)^2
      RMSE  <- sqrt(mean((y - pred)^2))
      RMSEP <- RMSE  # placeholder; could be replaced by CV-based metric
      
      ranking_df <- dplyr::bind_rows(ranking_df, dplyr::tibble(
        Model = mname,
        R2    = R2,
        RMSE  = RMSE,
        RMSEP = RMSEP
      ))
    }
    
    # =====================================================
    # RANK MODELS
    # =====================================================
    if (nrow(ranking_df) > 0) {
      ranking_df <- ranking_df %>%
        mutate(score = RMSE + RMSEP + (1 - R2)) %>%
        arrange(score) %>%
        mutate(Rank = dplyr::row_number()) %>%
        dplyr::select(Rank, Model, R2, RMSE, RMSEP)
      
      ranking_file <- file.path(element_dir, paste0(el, "_Model_Ranking.csv"))
      write.csv(ranking_df, ranking_file, row.names = FALSE)
      
      element_log <- paste0(
        element_log,
        "Model ranking written to: ", ranking_file, "\n"
      )
    } else {
      element_log <- paste0(element_log, "No valid models to rank.\n")
    }
    
    element_log <- paste0(
      element_log,
      "Element outputs saved in: ", element_dir, "\n"
    )
    write(element_log, file = log_file, append = TRUE)
    
  } # end for(elements)
  
  write("\n=== COMPLETED SUCCESSFULLY ===\n", file = log_file, append = TRUE)
  message("\nAll done. Outputs are under: ", save_dir)
}


# Input to run it! --------------------------------------------------------
elements <- c("Ti","Ca","Mn","Fe","Sr","Zr")
run_full_regressions(elements, ACE_dataset)