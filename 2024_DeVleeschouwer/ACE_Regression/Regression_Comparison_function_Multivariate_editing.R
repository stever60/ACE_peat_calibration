# **************************************************************************
# * Regression Comparison Multivariate function for XRF-CS log_inc dataset *
# *                                                                        *
# **************************************************************************
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
# - RF OOB MSE & R² 
# - PLS RMSEP at best ncomp
# -	Explanations for each test in plain English
# Makes a Best Model performance summary
# - ranking based on R2, RMSE, RMSEP where AIC & BIC dont exist (PLS, RF)
# - where R2 higher, RMSE/RMSEP lower, AIC, BIC lower is better
# — sequential R that works on a Mac.

# Random Forest OOB (Out-of-Bag) is a method for estimating a random forest 
# model's performance without using a separate test set. It works by having each 
# tree in the forest evaluate the data points that were not included in its own 
# training sample (the "out-of-bag" data). By averaging these predictions across 
# all trees and their respective OOB samples, the out-of-bag error provides a 
# robust and unbiased estimate of the model's generalization error)

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
# Install libraries needed ---------------------------------------------------------------
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
    data,          # calibration dataset, ln-transformed ICP (ln mg kg^-1)
    XRF_new = NULL,
    var_suffix = "_var",
    sd_suffix  = "_ICP_sd",
    save_dir   = "/Users/sjro/Desktop/Regression_multivariate_6",
    summary_all_csv = "AllElements_ModelSummary.csv",
    best_models_csv  = "BestModels_PerElement.csv"
) {
  
  ## 1. Packages ----
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(pls)
    library(randomForest)
    library(arm)
    library(ggsci)
    library(patchwork)
    library(lmtest)
    library(progress)
  })
  
  ## 2. Output dirs ----
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_base  <- file.path(save_dir, "All_models_ppm")
  best_base <- file.path(save_dir, "Best_models_ppm")
  
  dir.create(all_base,  recursive = TRUE, showWarnings = FALSE)
  dir.create(best_base, recursive = TRUE, showWarnings = FALSE)
  
  ## 3. Site factor ----
  if (!"Site" %in% names(data))
    stop("Calibration dataset must contain a 'Site' column.")
  
  data$Site <- factor(
    data$Site,
    levels = c("BI10","HER42PB","KER1","KER3","PB1")
  )
  
  ## 4. Theme ----
  theme_small <- theme(
    text       = element_text(size = 8),
    axis.text  = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.title   = element_text(size = 8, face = "bold")
  )
  
  ## 5. Models & colours ----
  model_levels <- c(
    "OLS","WLS","OLS_weighted","WLS_weighted",
    "Bayes","RF","PLS_LOO","PLS_kfold"
  )
  
  model_labels <- c(
    "OLS","WLS","OLS(wt)","WLS(wt)",
    "Bayes","RF","PLS(LOO)","PLS(k)"
  )
  
  cols_models <- ggsci::pal_npg("nrc")(length(model_labels))
  names(cols_models) <- model_labels
  
  # PLS(k) = salmon3
  cols_models["PLS(k)"] <- "salmon3"
  
  # ICP-MS observed in red
  cols_all <- c(cols_models, "ICP-MS (Observed)" = "red")
  
  # Multivariate models for “best” plots
  multivar_labels <- c("Bayes","RF","PLS(LOO)","PLS(k)")
  
  ## 6. Integrity check ----
  check_dataset_integrity <- function(data, el, y_col, x_col, var_col, sd_col) {
    errs <- c()
    
    req <- c(y_col, x_col)
    for (v in req) {
      if (!v %in% names(data))
        errs <- c(errs, paste("Missing required column:", v))
    }
    
    if (!is.numeric(data[[y_col]]))
      errs <- c(errs, paste(y_col, "is not numeric"))
    if (!is.numeric(data[[x_col]]))
      errs <- c(errs, paste(x_col, "is not numeric"))
    
    if (any(is.na(data[[y_col]])))
      errs <- c(errs, paste("NA values in", y_col))
    if (any(is.na(data[[x_col]])))
      errs <- c(errs, paste("NA values in", x_col))
    
    if (sd_col %in% names(data)) {
      if (any(data[[sd_col]] <= 0, na.rm = TRUE))
        errs <- c(errs, paste(sd_col, "contains non-positive values"))
    }
    
    errs
  }
  
  ## 7. 10-fold CV RMSE helper ----
  cv_rmse_for_model <- function(
    mname, formula, data, y_col,
    var_col, sd_col, best_ncomp = NA, k = 10
  ) {
    
    set.seed(123)
    n <- nrow(data)
    if (n < k) return(NA_real_)
    
    folds <- sample(rep(1:k, length.out = n))
    total_sse <- 0
    total_n   <- 0
    
    for (fold in 1:k) {
      
      train <- data[folds != fold, , drop = FALSE]
      test  <- data[folds == fold, , drop = FALSE]
      
      fit <- tryCatch({
        if (mname == "OLS") {
          lm(formula, data = train)
          
        } else if (mname == "WLS") {
          base <- lm(formula, data = train)
          r    <- residuals(base)
          rm   <- lm(abs(r) ~ fitted(base))
          w    <- 1/(rm$fitted.values^2)
          lm(formula, data = train, weights = w)
          
        } else if (mname == "OLS_weighted" && sd_col %in% names(train)) {
          w_sd <- 1/(train[[sd_col]]^2)
          lm(formula, data = train, weights = w_sd)
          
        } else if (mname == "WLS_weighted" && sd_col %in% names(train)) {
          w_sd <- 1/(train[[sd_col]]^2)
          base_sd <- lm(formula, data = train, weights = w_sd)
          rm_sd   <- lm(abs(base_sd$residuals) ~ base_sd$fitted.values)
          w2_sd   <- 1/(rm_sd$fitted.values^2)
          lm(formula, data = train, weights = w2_sd)
          
        } else if (mname == "Bayes") {
          bayesglm(formula, data = train)
          
        } else if (mname == "RF") {
          randomForest(formula, data = train, ntree = 300)
          
        } else if (mname %in% c("PLS_LOO","PLS_kfold")) {
          plsr(formula, data = train, scale = TRUE)
          
        } else NULL
      }, error = function(e) NULL)
      
      if (is.null(fit)) return(NA_real_)
      
      pred <- tryCatch({
        if (inherits(fit,"mvr"))
          as.numeric(predict(fit, ncomp = best_ncomp, newdata = test))
        else
          as.numeric(predict(fit, newdata = test))
      }, error = function(e) NULL)
      
      if (is.null(pred)) return(NA_real_)
      
      y_test <- test[[y_col]]
      ok <- !is.na(y_test) & !is.na(pred)
      if (!any(ok)) next
      
      total_sse <- total_sse + sum((y_test[ok] - pred[ok])^2)
      total_n   <- total_n   + sum(ok)
    }
    
    if (total_n == 0) return(NA_real_)
    sqrt(total_sse / total_n)
  }
  
  ## 8. Global summary + preds storage ----
  global_summary <- data.frame(
    Element = character(),
    Model   = character(),
    R2      = numeric(),
    RMSE    = numeric(),
    RMSEP   = numeric(),
    AIC     = numeric(),
    BIC     = numeric(),
    stringsAsFactors = FALSE
  )
  
  preds_store <- list()
  
  
  ## 9. Loop over elements ----
  pb_elements <- progress_bar$new(
    format = "Processing elements [:bar] :current/:total (:percent) eta::eta",
    total  = length(elements),
    clear  = FALSE,
    width  = 80
  )
  
  for (el in elements) {
    
    pb_elements$tick()
    
    y_col  <- paste0(el, "_ICP")       # ln(mg kg^-1)
    x_col  <- el                       # ln(XRF)
    var_col <- paste0(el, var_suffix)
    sd_col  <- paste0(el, sd_suffix)
    
    element_dir <- file.path(all_base, el)
    dir.create(element_dir, recursive = TRUE, showWarnings = FALSE)
    
    errs <- check_dataset_integrity(data, el, y_col, x_col, var_col, sd_col)
    if (length(errs) > 0) {
      warning(paste("Skipping", el, ":", paste(errs, collapse = "; ")))
      next
    }
    
    fmla_single <- as.formula(paste0(y_col, " ~ ", x_col))
    fmla_multi  <- as.formula(
      paste0(y_col, " ~ ", paste(elements, collapse = " + "))
    )
    
    y <- data[[y_col]]
    
    ## 9.1 Fit models ----
    model_list <- list()
    
    # OLS
    model_list$OLS <- lm(fmla_single, data = data)
    
    # WLS (residual-based)
    {
      base <- model_list$OLS
      r    <- residuals(base)
      rm   <- lm(abs(r) ~ fitted(base))
      w    <- 1/(rm$fitted.values^2)
      model_list$WLS <- lm(fmla_single, data = data, weights = w)
    }
    
    # Weighted models if sd_col exists
    if (sd_col %in% names(data)) {
      w_sd <- 1/(data[[sd_col]]^2)
      
      model_list$OLS_weighted <- lm(fmla_single, data = data, weights = w_sd)
      
      base_sd <- lm(fmla_single, data = data, weights = w_sd)
      rm_sd   <- lm(abs(base_sd$residuals) ~ base_sd$fitted.values)
      w2_sd   <- 1/(rm_sd$fitted.values^2)
      model_list$WLS_weighted <- lm(fmla_single, data = data, weights = w2_sd)
    } else {
      model_list$OLS_weighted <- NULL
      model_list$WLS_weighted <- NULL
    }
    
    # Bayes (multivariate)
    model_list$Bayes <- bayesglm(fmla_multi, data = data)
    
    # Random Forest (multivariate)
    set.seed(42)
    model_list$RF <- randomForest(
      fmla_multi,
      data = data,
      ntree = 500,
      importance = TRUE
    )
    
    # PLS_LOO (multivariate)
    pls_loo <- plsr(
      fmla_multi,
      data = data,
      scale = TRUE,
      validation = "LOO"
    )
    best_loo <- selectNcomp(pls_loo, method = "onesigma", plot = FALSE)
    attr(pls_loo, "best_ncomp") <- best_loo
    model_list$PLS_LOO <- pls_loo
    
    # PLS_kfold (multivariate)
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
    
    ## 9.2 Per-element tables ----
    ranking_df <- data.frame(
      Model = character(),
      R2    = numeric(),
      RMSE  = numeric(),
      RMSEP = numeric(),
      AIC   = numeric(),
      BIC   = numeric(),
      stringsAsFactors = FALSE
    )
    
    diag_df <- data.frame(
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
      BIC_expl          = character(),
      stringsAsFactors  = FALSE
    )
    
    do_predictions <- !is.null(XRF_new)
    preds_df <- if (do_predictions) XRF_new else NULL
    
    ## 9.3 Model loop ----
    for (mname in names(model_list)) {
      
      fit <- model_list[[mname]]
      if (is.null(fit)) next
      
      pred_log <- tryCatch({
        if (inherits(fit,"mvr")) {
          ncomp <- attr(fit,"best_ncomp")
          as.numeric(predict(fit, ncomp = ncomp, newdata = data))
        } else {
          as.numeric(predict(fit, newdata = data))
        }
      }, error = function(e) rep(NA_real_, nrow(data)))
      
      ok <- !is.na(pred_log) & !is.na(y)
      if (!any(ok)) next
      
      SS_res <- sum((y[ok] - pred_log[ok])^2)
      SS_tot <- sum((y[ok] - mean(y[ok]))^2)
      R2     <- 1 - SS_res / SS_tot
      RMSE   <- sqrt(mean((y[ok] - pred_log[ok])^2))
      
      best_ncomp <- if (inherits(fit,"mvr")) attr(fit,"best_ncomp") else NA
      
      RMSEP10 <- cv_rmse_for_model(
        mname      = mname,
        formula    = if (mname %in% c("OLS","WLS","OLS_weighted","WLS_weighted"))
          fmla_single else fmla_multi,
        data       = data,
        y_col      = y_col,
        var_col    = var_col,
        sd_col     = sd_col,
        best_ncomp = best_ncomp,
        k          = 10
      )
      
      if (mname %in% c("RF","PLS_LOO","PLS_kfold")) {
        AICv <- NA_real_
        BICv <- NA_real_
      } else {
        AICv <- tryCatch(AIC(fit), error = function(e) NA_real_)
        BICv <- tryCatch(BIC(fit), error = function(e) NA_real_)
      }
      
      sh <- tryCatch(shapiro.test(residuals(fit)), error = function(e) NULL)
      if (is.null(sh)) {
        ShW <- NA_real_; Shp <- NA_real_
        Sh_expl <- "Shapiro–Wilk not applicable."
      } else {
        ShW <- as.numeric(sh$statistic)
        Shp <- sh$p.value
        Sh_expl <- if (Shp < 0.05)
          "Residuals deviate from normality."
        else
          "Residuals approx. normal."
      }
      
      bp <- tryCatch(bptest(fit), error = function(e) NULL)
      if (is.null(bp)) {
        BPstat <- NA_real_; BPp <- NA_real_
        BP_expl <- "Breusch–Pagan not applicable."
      } else {
        BPstat <- as.numeric(bp$statistic)
        BPp    <- bp$p.value
        BP_expl <- if (BPp < 0.05)
          "Heteroscedastic residuals."
        else
          "No evidence of heteroscedasticity."
      }
      
      CV_expl  <- if (is.na(RMSEP10)) "10-fold CV not available." else "10-fold CV computed."
      AIC_expl <- if (is.na(AICv))    "AIC not available."       else "AIC computed."
      BIC_expl <- if (is.na(BICv))    "BIC not available."       else "BIC computed."
      
      diag_df <- rbind(
        diag_df,
        data.frame(
          Model             = mname,
          Shapiro_W         = ShW,
          Shapiro_p         = Shp,
          Shapiro_expl      = Sh_expl,
          BreuschPagan_stat = BPstat,
          BreuschPagan_p    = BPp,
          BreuschPagan_expl = BP_expl,
          CV10_RMSE         = RMSEP10,
          CV10_expl         = CV_expl,
          AIC               = AICv,
          AIC_expl          = AIC_expl,
          BIC               = BICv,
          BIC_expl          = BIC_expl,
          stringsAsFactors  = FALSE
        )
      )
      
      runlab <- dplyr::recode(
        mname,
        OLS         = "OLS",
        WLS         = "WLS",
        OLS_weighted= "OLS(wt)",
        WLS_weighted= "WLS(wt)",
        Bayes       = "Bayes",
        RF          = "RF",
        PLS_LOO     = "PLS(LOO)",
        PLS_kfold   = "PLS(k)"
      )
      
      ranking_df <- rbind(
        ranking_df,
        data.frame(
          Model = runlab,
          R2    = R2,
          RMSE  = RMSE,
          RMSEP = RMSEP10,
          AIC   = AICv,
          BIC   = BICv,
          stringsAsFactors = FALSE
        )
      )
      
      ## Predictions for XRF_new in ppm ----
      if (do_predictions) {
        pred_new_log <- tryCatch({
          if (inherits(fit,"mvr")) {
            ncomp <- attr(fit,"best_ncomp")
            as.numeric(predict(fit, ncomp = ncomp, newdata = XRF_new))
          } else {
            as.numeric(predict(fit, newdata = XRF_new))
          }
        }, error = function(e) rep(NA_real_, nrow(XRF_new)))
        
        lwr_log <- pred_new_log - 1.96 * RMSE
        upr_log <- pred_new_log + 1.96 * RMSE
        
        preds_df[[paste0(runlab,"_",el,"_Pred_ppm")]] <- exp(pred_new_log)
        preds_df[[paste0(runlab,"_",el,"_L95_ppm")]]  <- exp(lwr_log)
        preds_df[[paste0(runlab,"_",el,"_U95_ppm")]]  <- exp(upr_log)
      }
      
      ## Pred vs obs (log-space) ----
      stats_label <- paste0(
        "R² = ", signif(R2,3),
        "\nRMSEP = ", signif(RMSEP10,3),
        "\nRMSE = ", signif(RMSE,3)
      )
      
      df_po <- data.frame(Observed = y, Predicted = pred_log)
      x_min_po <- min(df_po$Observed,  na.rm = TRUE)
      y_max_po <- max(df_po$Predicted, na.rm = TRUE)
      
      p_po <- ggplot(df_po, aes(Observed, Predicted)) +
        geom_point(size = 1, alpha = 0.85) +
        geom_abline(slope = 1, intercept = 0) +
        annotate(
          "text",
          x = x_min_po, y = y_max_po,
          hjust = 0, vjust = 1,
          label = stats_label,
          size = 2.7
        ) +
        theme_bw() + theme_small +
        ggtitle(paste(el, "-", runlab, "Predicted vs Observed (log-space)"))
      
      ggsave(
        file.path(element_dir, paste0(el,"_",runlab,"_PredVsObs_log.pdf")),
        p_po, width = 13, height = 9, units = "cm"
      )
      
      ## Residuals (if numeric) ----
      res <- tryCatch(as.numeric(residuals(fit)), error = function(e) NULL)
      
      if (!is.null(res) &&
          is.numeric(res) &&
          length(res) > 1 &&
          !all(is.na(res))) {
        
        hist_obj  <- hist(res, plot = FALSE)
        x_min_res <- min(res, na.rm = TRUE)
        y_max_res <- max(hist_obj$counts, na.rm = TRUE)
        
        res_df <- data.frame(resid = res)
        
        p_res <- ggplot(res_df, aes(resid)) +
          geom_histogram(bins = 30, fill = "grey80") +
          geom_density(color = "red") +
          annotate(
            "text",
            x = x_min_res, y = y_max_res,
            hjust = 0, vjust = 1,
            label = stats_label,
            size = 2.7
          ) +
          theme_bw() + theme_small +
          ggtitle(paste(el, "-", runlab, "Residuals (log-space)"))
        
        ggsave(
          file.path(element_dir, paste0(el,"_",runlab,"_Residuals_log.pdf")),
          p_res, width = 13, height = 9, units = "cm"
        )
      }
      
      ## Influence (lm only) ----
      if (inherits(fit,"lm") && !inherits(fit,"bayesglm")) {
        infl_dir <- file.path(element_dir, "Influence")
        dir.create(infl_dir, showWarnings = FALSE)
        
        cooks <- cooks.distance(fit)
        lev   <- hatvalues(fit)
        rstd  <- rstandard(fit)
        
        infl_df <- data.frame(
          obs      = seq_along(cooks),
          Cook     = cooks,
          Leverage = lev,
          StdResid = rstd
        )
        
        p_cook <- ggplot(infl_df, aes(obs, Cook)) +
          geom_col() +
          theme_bw() + theme_small +
          ggtitle(paste(el, "-", runlab, "Cook's distance"))
        
        p_lev <- ggplot(infl_df, aes(Leverage, StdResid)) +
          geom_point(size = 1) +
          theme_bw() + theme_small +
          ggtitle(paste(el, "-", runlab, "Leverage vs Std. Residual"))
        
        pdf(
          file.path(infl_dir, paste0(el,"_",runlab,"_Influence.pdf")),
          width = 13/2.54, height = 9/2.54
        )
        print(p_cook)
        print(p_lev)
        dev.off()
      }
      
      ## PLS RMSEP ----
      if (inherits(fit,"mvr")) {
        pdf(
          file.path(element_dir, paste0(el,"_",runlab,"_PLS_RMSEP.pdf")),
          width = 13/2.54, height = 9/2.54
        )
        validationplot(
          fit,
          val.type = "RMSEP",
          main = paste(el,"-",runlab,"RMSEP vs components")
        )
        dev.off()
      }
      
    } # end model loop
    
    ## 9.4 Ranking, diagnostics, preds ----
    if (nrow(ranking_df) > 0) {
      
      ranking_df <- ranking_df[order(-ranking_df$R2,
                                     ranking_df$RMSE,
                                     ranking_df$RMSEP), ]
      ranking_df$Rank <- seq_len(nrow(ranking_df))
      
      write.csv(
        ranking_df,
        file.path(element_dir, paste0(el, "_Model_Ranking.csv")),
        row.names = FALSE
      )
      
      mean_log <- mean(y, na.rm = TRUE)
      mean_ppm <- exp(mean_log)
      
      ranking_ppm <- ranking_df
      ranking_ppm$RMSE_ppm_factor  <- exp(ranking_ppm$RMSE)
      ranking_ppm$RMSEP_ppm_factor <- exp(ranking_ppm$RMSEP)
      ranking_ppm$Mean_log_ppm     <- mean_log
      ranking_ppm$Mean_ppm         <- mean_ppm
      ranking_ppm$RMSE_ppm_abs     <- ranking_ppm$RMSE  * mean_ppm
      ranking_ppm$RMSEP_ppm_abs    <- ranking_ppm$RMSEP * mean_ppm
      
      write.csv(
        ranking_ppm,
        file.path(element_dir, paste0(el, "_Model_Ranking_ppm.csv")),
        row.names = FALSE
      )
      
      add_global <- ranking_df
      add_global$Element <- el
      add_global <- add_global[, c("Element","Model","R2","RMSE","RMSEP","AIC","BIC")]
      
      global_summary <- rbind(global_summary, add_global)
    }
    
    if (nrow(diag_df) > 0) {
      write.csv(
        diag_df,
        file.path(element_dir, paste0(el, "_Diagnostics_Summary.csv")),
        row.names = FALSE
      )
    }
    
    if (do_predictions) {
      write.csv(
        preds_df,
        file.path(element_dir, paste0(el, "_preds.csv")),
        row.names = FALSE
      )
      preds_store[[el]] <- preds_df
    }
    
  } # end element loop
  
  ############################################################
  # 10. SITE-LEVEL DEPTH & AGE PLOTS
  ############################################################
  
  if (!is.null(XRF_new) && length(preds_store) > 0 &&
      "Site" %in% names(XRF_new)) {
    
    has_depth   <- "depth" %in% names(XRF_new)
    has_age_new <- "SH20_age" %in% names(XRF_new)
    has_age_cal <- "SH20_mean_age" %in% names(data)
    
    sites_dir_all  <- file.path(all_base,  "Sites")
    sites_dir_best <- file.path(best_base, "Sites")
    dir.create(sites_dir_all,  recursive = TRUE, showWarnings = FALSE)
    dir.create(sites_dir_best, recursive = TRUE, showWarnings = FALSE)
    
    sites <- sort(unique(XRF_new$Site))
    
    for (s in sites) {
      
      ## 10.1 Depth plots (2×3) ----
      depth_panels_all  <- list()
      depth_panels_best <- list()
      
      if (has_depth) {
        
        for (el in elements) {
          
          df_pred <- preds_store[[el]]
          if (is.null(df_pred)) next
          if (!all(c("Site","depth") %in% names(df_pred))) next
          
          df_s <- df_pred[df_pred$Site == s, , drop = FALSE]
          if (nrow(df_s) == 0) next
          
          pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                            names(df_s), value = TRUE)
          if (length(pred_cols) == 0) next
          
          df_long <- df_s[, c("depth","Site", pred_cols), drop = FALSE] %>%
            tidyr::pivot_longer(
              cols = all_of(pred_cols),
              names_to = "Model",
              values_to = "Pred_ppm"
            )
          
          df_long$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long$Model)
          df_long$Model <- factor(df_long$Model,
                                  levels = model_levels,
                                  labels = model_labels)
          
          df_long_best <- df_long[df_long$Model %in% multivar_labels, , drop = FALSE]
          
          y_col <- paste0(el, "_ICP")
          ace_overlay <- data[data$Site == s &
                                !is.na(data$depth) &
                                !is.na(data[[y_col]]), ]
          if (nrow(ace_overlay) > 0) {
            ace_overlay$ICP_ppm <- exp(ace_overlay[[y_col]])
          } else {
            ace_overlay <- NULL
          }
          
          p_all <- ggplot(df_long,
                          aes(x = Pred_ppm, y = depth,
                              colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_y_reverse() +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = names(cols_all),
              guide  = guide_legend(nrow = 2)
            ) +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Depth (cm)",
              title = paste0("Site ", s, " — ", el, " (depth, all models)")
            ) +
            theme_bw() + theme_small
          
          p_best <- ggplot(df_long_best,
                           aes(x = Pred_ppm, y = depth,
                               colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_y_reverse() +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = c(multivar_labels, "ICP-MS (Observed)"),
              guide  = guide_legend(nrow = 2)
            ) +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Depth (cm)",
              title = paste0("Site ", s, " — ", el, " (depth, multivariate)")
            ) +
            theme_bw() + theme_small
          
          if (!is.null(ace_overlay)) {
            p_all <- p_all +
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
            
            p_best <- p_best +
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
          
          depth_panels_all[[el]]  <- p_all
          depth_panels_best[[el]] <- p_best
        }
        
        if (length(depth_panels_all) > 0) {
          ord_all <- lapply(elements, function(e) depth_panels_all[[e]])
          ord_all <- lapply(ord_all, function(p)
            if (is.null(p)) ggplot() + theme_void() else p
          )
          
          depth_all <- (ord_all[[1]] | ord_all[[2]]) /
            (ord_all[[3]] | ord_all[[4]]) /
            (ord_all[[5]] | ord_all[[6]])
          
          depth_all <- depth_all +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          ggsave(
            file.path(sites_dir_all,
                      paste0("Site_", s, "_Profiles_A4_depth_all.pdf")),
            depth_all, width = 21, height = 29.7, units = "cm"
          )
        }
        
        if (length(depth_panels_best) > 0) {
          ord_best <- lapply(elements, function(e) depth_panels_best[[e]])
          ord_best <- lapply(ord_best, function(p)
            if (is.null(p)) ggplot() + theme_void() else p
          )
          
          depth_best <- (ord_best[[1]] | ord_best[[2]]) /
            (ord_best[[3]] | ord_best[[4]]) /
            (ord_best[[5]] | ord_best[[6]])
          
          depth_best <- depth_best +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          ggsave(
            file.path(sites_dir_best,
                      paste0("Site_", s, "_Profiles_A4_depth_multivariate.pdf")),
            depth_best, width = 21, height = 29.7, units = "cm"
          )
        }
      }
      
      ## 10.2 Age plots (2×3) ----
      if (has_age_new && has_age_cal) {
        
        age_panels_all  <- list()
        age_panels_best <- list()
        
        for (el in elements) {
          
          df_pred <- preds_store[[el]]
          if (is.null(df_pred)) next
          if (!all(c("Site","SH20_age") %in% names(df_pred))) next
          
          df_s <- df_pred[df_pred$Site == s, , drop = FALSE]
          if (nrow(df_s) == 0) next
          
          pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                            names(df_s), value = TRUE)
          if (length(pred_cols) == 0) next
          
          df_long <- df_s[, c("SH20_age","Site", pred_cols), drop = FALSE] %>%
            tidyr::pivot_longer(
              cols = all_of(pred_cols),
              names_to = "Model",
              values_to = "Pred_ppm"
            )
          
          df_long$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long$Model)
          df_long$Model <- factor(df_long$Model,
                                  levels = model_levels,
                                  labels = model_labels)
          
          df_long_best <- df_long[df_long$Model %in% multivar_labels, , drop = FALSE]
          
          y_col <- paste0(el, "_ICP")
          ace_overlay_age <- data[data$Site == s &
                                    !is.na(data$SH20_mean_age) &
                                    !is.na(data[[y_col]]), ]
          if (nrow(ace_overlay_age) > 0) {
            ace_overlay_age$ICP_ppm <- exp(ace_overlay_age[[y_col]])
          } else {
            ace_overlay_age <- NULL
          }
          
          p_all <- ggplot(df_long,
                          aes(x = Pred_ppm, y = SH20_age,
                              colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = names(cols_all),
              guide  = guide_legend(nrow = 2)
            ) +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Age (cal a BP)",
              title = paste0("Site ", s, " — ", el, " (age vs mg kg^-1, all models)")
            ) +
            theme_bw() + theme_small
          
          p_best <- ggplot(df_long_best,
                           aes(x = Pred_ppm, y = SH20_age,
                               colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = c(multivar_labels, "ICP-MS (Observed)"),
              guide  = guide_legend(nrow = 2)
            ) +
            labs(
              x = bquote(.(el)~"predicted (mg"~kg^{-1}*")"),
              y = "Age (cal a BP)",
              title = paste0("Site ", s, " — ", el, " (age vs mg kg^-1, multivariate)")
            ) +
            theme_bw() + theme_small
          
          if (!is.null(ace_overlay_age)) {
            p_all <- p_all +
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
            
            p_best <- p_best +
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
          
          age_panels_all[[el]]  <- p_all
          age_panels_best[[el]] <- p_best
        }
        
        if (length(age_panels_all) > 0) {
          ord_all <- lapply(elements, function(e) age_panels_all[[e]])
          ord_all <- lapply(ord_all, function(p)
            if (is.null(p)) ggplot() + theme_void() else p
          )
          
          age_all <- (ord_all[[1]] | ord_all[[2]]) /
            (ord_all[[3]] | ord_all[[4]]) /
            (ord_all[[5]] | ord_all[[6]])
          
          age_all <- age_all +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          ggsave(
            file.path(sites_dir_all,
                      paste0("Site_", s, "_Profiles_A4_age_all.pdf")),
            age_all, width = 21, height = 29.7, units = "cm"
          )
        }
        
        if (length(age_panels_best) > 0) {
          ord_best <- lapply(elements, function(e) age_panels_best[[e]])
          ord_best <- lapply(ord_best, function(p)
            if (is.null(p)) ggplot() + theme_void() else p
          )
          
          age_best <- (ord_best[[1]] | ord_best[[2]]) /
            (ord_best[[3]] | ord_best[[4]]) /
            (ord_best[[5]] | ord_best[[6]])
          
          age_best <- age_best +
            plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
          
          ggsave(
            file.path(sites_dir_best,
                      paste0("Site_", s, "_Profiles_A4_age_multivariate.pdf")),
            age_best, width = 21, height = 29.7, units = "cm"
          )
        }
      }
      
    } # end site loop
  } # end if XRF_new & preds_store
  
  ############################################################
  # 11. ELEMENT-LEVEL MULTI-SITE AGE–PPM PLOTS
  ############################################################
  
  if (!is.null(XRF_new) &&
      length(preds_store) > 0 &&
      "Site" %in% names(XRF_new) &&
      "SH20_age" %in% names(XRF_new) &&
      "SH20_mean_age" %in% names(data)) {
    
    elem_dir_all  <- file.path(all_base,  "Elements")
    elem_dir_best <- file.path(best_base, "Elements")
    dir.create(elem_dir_all,  recursive = TRUE, showWarnings = FALSE)
    dir.create(elem_dir_best, recursive = TRUE, showWarnings = FALSE)
    
    age_vals <- c(XRF_new$SH20_age, data$SH20_mean_age)
    age_vals <- age_vals[!is.na(age_vals)]
    if (length(age_vals) > 0) {
      min_age <- min(age_vals)
      max_age <- max(age_vals)
    } else {
      min_age <- NA_real_
      max_age <- NA_real_
    }
    
    site_levels_all <- sort(unique(c(as.character(XRF_new$Site),
                                     as.character(data$Site))))
    
    for (el in elements) {
      
      df_pred <- preds_store[[el]]
      if (is.null(df_pred)) next
      
      pred_cols <- grep(paste0("_", el, "_Pred_ppm$"),
                        names(df_pred), value = TRUE)
      if (length(pred_cols) == 0) next
      
      df_long_all <- df_pred[!is.na(df_pred$SH20_age),
                             c("Site","SH20_age", pred_cols), drop = FALSE] %>%
        tidyr::pivot_longer(
          cols = all_of(pred_cols),
          names_to = "Model",
          values_to = "Pred_ppm"
        )
      
      df_long_all$Model <- sub(paste0("_", el, "_Pred_ppm$"), "", df_long_all$Model)
      df_long_all$Model <- factor(df_long_all$Model,
                                  levels = model_levels,
                                  labels = model_labels)
      
      df_long_best <- df_long_all[df_long_all$Model %in% multivar_labels, , drop = FALSE]
      
      y_col <- paste0(el, "_ICP")
      ace_all <- data[!is.na(data$SH20_mean_age) &
                        !is.na(data[[y_col]]), ]
      if (nrow(ace_all) > 0) {
        ace_all$ICP_ppm <- exp(ace_all[[y_col]])
      } else {
        ace_all <- NULL
      }
      
      plot_list_all  <- list()
      plot_list_best <- list()
      
      for (s in site_levels_all) {
        
        df_all_s  <- df_long_all[df_long_all$Site == s, , drop = FALSE]
        df_best_s <- df_long_best[df_long_best$Site == s, , drop = FALSE]
        
        ace_s <- if (!is.null(ace_all)) ace_all[ace_all$Site == s, , drop = FALSE] else NULL
        
        ## All models panel ----
        if (nrow(df_all_s) == 0 && (is.null(ace_s) || nrow(ace_s) == 0)) {
          plot_list_all[[s]] <- ggplot() + theme_void()
        } else {
          p_all <- ggplot(df_all_s,
                          aes(x = SH20_age, y = Pred_ppm,
                              colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = names(cols_all),
              guide  = guide_legend(nrow = 2)
            ) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste0("Site ", s, " (all models)")
            ) +
            theme_bw() + theme_small
          
          if (!is.na(min_age) && !is.na(max_age)) {
            p_all <- p_all + coord_cartesian(xlim = c(min_age, max_age))
          }
          
          if (!is.null(ace_s) && nrow(ace_s) > 0) {
            p_all <- p_all +
              geom_path(
                data = ace_s,
                aes(x = SH20_mean_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_s,
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
          
          plot_list_all[[s]] <- p_all
        }
        
        ## Multivariate-only panel ----
        if (nrow(df_best_s) == 0 && (is.null(ace_s) || nrow(ace_s) == 0)) {
          plot_list_best[[s]] <- ggplot() + theme_void()
        } else {
          p_best <- ggplot(df_best_s,
                           aes(x = SH20_age, y = Pred_ppm,
                               colour = Model, group = Model)) +
            geom_path(linewidth = 0.4) +
            scale_color_manual(
              name   = "Model",
              values = cols_all,
              breaks = c(multivar_labels, "ICP-MS (Observed)"),
              guide  = guide_legend(nrow = 2)
            ) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste0("Site ", s, " (multivariate)")
            ) +
            theme_bw() + theme_small
          
          if (!is.na(min_age) && !is.na(max_age)) {
            p_best <- p_best + coord_cartesian(xlim = c(min_age, max_age))
          }
          
          if (!is.null(ace_s) && nrow(ace_s) > 0) {
            p_best <- p_best +
              geom_path(
                data = ace_s,
                aes(x = SH20_mean_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)",
                    group  = "ICP-MS (Observed)"),
                inherit.aes = FALSE,
                linewidth = 0.6
              ) +
              geom_point(
                data = ace_s,
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
          
          plot_list_best[[s]] <- p_best
        }
        
      } # end per-site loop
      
      if (length(plot_list_all) > 0) {
        combined_all <- Reduce(`/`, plot_list_all)
        combined_all <- combined_all +
          plot_layout(guides = "collect") &
          theme(legend.position = "bottom")
        
        ggsave(
          file.path(elem_dir_all,
                    paste0(el, "_model_age_comparison_all.pdf")),
          combined_all, width = 21, height = 29.7, units = "cm"
        )
      }
      
      if (length(plot_list_best) > 0) {
        combined_best <- Reduce(`/`, plot_list_best)
        combined_best <- combined_best +
          plot_layout(guides = "collect") &
          theme(legend.position = "bottom")
        
        ggsave(
          file.path(elem_dir_best,
                    paste0(el, "_model_age_comparison_multivariate.pdf")),
          combined_best, width = 21, height = 29.7, units = "cm"
        )
      }
      
    } # end element loop
  } # end if multi-site age available
  
  ############################################################
  # 12. GLOBAL SUMMARY TABLES (LOG + PPM)
  ############################################################
  
  if (nrow(global_summary) > 0) {
    
    global_ranked <- global_summary[order(global_summary$Element,
                                          -global_summary$R2,
                                          global_summary$RMSE,
                                          global_summary$RMSEP), ]
    
    global_ranked$Rank <- ave(
      global_ranked$R2,
      global_ranked$Element,
      FUN = function(z) seq_along(z)
    )
    
    write.csv(
      global_ranked,
      file.path(save_dir, summary_all_csv),
      row.names = FALSE
    )
    
    best_models <- subset(global_ranked, Rank == 1)
    write.csv(
      best_models,
      file.path(save_dir, best_models_csv),
      row.names = FALSE
    )
    
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
      ml <- mean(vals)
      mp <- exp(ml)
      data.frame(
        Element      = el,
        Mean_log_ppm = ml,
        Mean_ppm     = mp
      )
    })
    
    elem_means_df <- do.call(rbind, elem_means_list)
    
    global_ext <- merge(global_ranked, elem_means_df,
                        by = "Element", all.x = TRUE)
    
    global_ext$RMSE_ppm_factor  <- exp(global_ext$RMSE)
    global_ext$RMSEP_ppm_factor <- exp(global_ext$RMSEP)
    global_ext$RMSE_ppm_abs     <- global_ext$RMSE  * global_ext$Mean_ppm
    global_ext$RMSEP_ppm_abs    <- global_ext$RMSEP * global_ext$Mean_ppm
    
    cols_keep <- c(
      "Element","Model","Rank","R2",
      "RMSE","RMSE_ppm_factor","RMSE_ppm_abs",
      "RMSEP","RMSEP_ppm_factor","RMSEP_ppm_abs",
      "Mean_log_ppm","Mean_ppm",
      "AIC","BIC"
    )
    
    global_summary_ppm <- global_ext[, cols_keep]
    
    write.csv(
      global_summary_ppm,
      file.path(save_dir, "AllElements_ModelSummary_ppm.csv"),
      row.names = FALSE
    )
  }
  
  ############################################################
  # 13. DONE
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
list.files("/Users/sjro/Desktop/Regression_multivariate_6/Ti", full.names = TRUE)
read.csv("/Users/sjro/Desktop/Regression_multivariate_6/Ti/Ti_Model_Ranking.csv")
# -------------------------------------------------------------------------



