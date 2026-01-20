# -------------------------------------------------------------------------
# Calibration Model Comparison - Multivariate function
# Using Matched XRF-CS log_inc & log ICP-MS matched dataset
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# ReadMe description and updates  ----------------------------------------

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
# • Predicted vs observed plots coloured by Site using a JCO palette
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


# To fix in v. 2.1
# Radar plots print to screen but not to pdf even though pdf created 
# - couldnt resolve issue but not a big deal .. just nice to have 


# Improvements in v.2.1-2.2 ----------------------------------------------------

# Added Diagnostics.csv file writing - accidentally deleted!

# CONFIDENCE FLAGGING * NEW in version 2.1
# Interpretation
# •	Uncertainty-dominated → prediction intervals dominate variance
# •	Noise-dominated → unstable downcore structure
# •	Signal-dominated → coherent downcore signal with constrained uncertainty
# 
# This language aligns with palaeoclimate proxy evaluation literature.
#
# # Strong model: high robustness, not unstable
# Robustness_Score >= 0.75 &
#   !Unstable_flag ~
#   "High confidence (signal-dominated)",
# 
# # Uncertainty-dominated: wide CI relative to variance
# SNR_CI_s < 0.4 &
#   Robustness_Score < 0.6 ~
#   "Low confidence (uncertainty-dominated)",
# 
# # Noise-dominated: rough downcore behaviour
# Unstable_flag == TRUE ~
#   "Low confidence (noise-dominated)",
# 
# # Intermediate but usable
# Robustness_Score >= 0.6 &
#   Robustness_Score < 0.75 ~
#   "Moderate confidence",

# Improvements in v.2.0 ----------------------------------------------------
# Model robustness was evaluated using a composite signal-to-noise framework 
# incorporating prediction uncertainty, downcore smoothness, and calibration 
# performance. Signal-to-noise ratios were calculated from prediction confidence 
# intervals and profile roughness, scaled within each element, and combined with
# R² into a weighted robustness score overall and per site. 
# Overall robustness across the whole calibration dataset was useful, but models 
# exhibiting unstable downcore behaviour were assessed on a per site basis and 
# excluded from selection if classified as unstable (robustness <0.75). 
# For each element, and at each site, the highest-ranked 
# stable model was designated as the production model and used for subsequent 
# interpretation.

# Metrics used
# 
# Quantity	Definition
# Signal 	  Variability of predicted concentrations
# Noise 	  Model uncertainty + downcore roughness
# Scale	    ppm (mg kg⁻¹), not log
# 
# Metrics computed per element × model:
# 1.	SNR_model = mean predicted ppm / RMSE_ppm
# 2.	SNR_CI (preferred) = sd(predicted ppm) / mean(U95 - L95)
# 3.	SNR_smooth (downcore stability) = sd(predicted ppm) / sd(diff(predicted ppm)
# 
# 
# For each element × model, robustness is derived from:
# 1.	Accuracy
# •	R² (calibration or validation performance)
# 2.	Noise / uncertainty
# •	SNR_CI
# Signal-to-noise ratio based on prediction uncertainty
# (signal variability ÷ mean CI width)
# 3.	Downcore stability
# •	SNR_smooth
# Penalises high point-to-point roughness (spiky profiles)
# 
# 
# Robustness score
# 
# Each metric is scaled 0–1 within each element, then combined:
#   
# Robustness Score = 0.4 × SNR_CI_scaled + 0.3 × SNR_smooth_scaled + 0.3 × R²_scaled
# 
# Models are ranked per element by:
# 1.	Highest SNR_CI
# 2.	Highest SNR_smooth
# 3.	Highest R²
# 
# This ensures:
# •	Models with narrow CIs,
# •	smooth downcore behaviour, and
# •	high explanatory power
# are favoured.
# 
# Class	Definition
# Preferred:	Highest robustness score and not flagged unstable
# Acceptable:	Rank 2–3 robustness, stable but slightly noisier
# Unstable:	High roughness or poor SNR despite reasonable R²
# 
# Unstable models are never selected as production models, even if R² is high.

# Generated a summary Heatmap matrix for visualising robustness quickly for all
# models vs elements and radar plots for each element

# Fixed performace ranking issue - rank now appears in global summary csv files 
# Small bug fixes and pacthes to make plotting include full error bar range
# Major change to ensure erro bars are cources and plpotted correctly in 
# Sections 8 and 15 

# Prediction uncertainty in Predicted vs Observed plots (Sections 8 & 15)
# Predicted vs Observed (Pred–Obs) plots include uncertainty information that
# reflects both measurement error and model prediction uncertainty, handled
# consistently but differently between calibration (Section 8) and external
# validation / ACE prediction (Section 15).
# 
# Observed values (horizontal error bars)
# Horizontal error bars always represent ±2 × ICP-MS analytical standard
# deviation (ICP_sd) in ppm - remains measured ±1 SD in log scale model
# This reflects uncertainty in the reference (observed) ICP-MS measurements
# plotted on the x-axis and relates to 95% CI 
# 
# Predicted values (vertical error bars):
# Vertical error bars represent 95% prediction uncertainty of the calibrated
# model, centred on each predicted value.
# •	For linear models (OLS, WLS, weighted variants), 95% prediction intervals
# are taken directly from the fitted model where available.
# •	For non-linear or multivariate models (Bayesian regression, Random Forest, PLS),
# prediction uncertainty is approximated using ±1.96 × RMSE, assuming normally
# distributed prediction errors.
# •	All prediction intervals are converted to the plotting scale
# (log space in Section 8; ppm in Section 15) before display.
# 
# Error bar centring:
# All vertical prediction error bars are explicitly constructed so that
# upper and lower bounds are symmetric around the predicted value, ensuring
# that error bars are visually centred on each data point.
# 
# Plot scaling and geometry:
# •	All Pred–Obs plots use equal x and y axes and square panels,
# enforcing a true 1:1 reference.
# •	Axis limits are expanded to fully include the extent of both horizontal
# and vertical error bars, with a small padding to avoid clipping.
# 
# Section-specific notes:
# - Section 8 (Calibration):
#   Pred–Obs plots are shown in log space, comparing observed vs
#   fitted ICP-MS values during model calibration.
# Section 15 (ACE / External validation):
#   Pred–Obs plots are shown in ppm, comparing observed ICP-MS values to model
#   predictions applied to an independent dataset.
# 
# This approach ensures that uncertainty shown in Pred–Obs plots is
# model-appropriate, scale-consistent, and directly comparable across elements
# and modelling methods.

# Improvements in v1 (1.2-1.9)  ------------------------------------------------

# v. 1.9
# Relabled axes and plots for predicted vs obs

# Caption example for Training Model
# Summary predicted (measured XRF-CS matched Ln[Ti/inc.]) vs observed 
# (measured ICP-MS matched Ln[Ti]) log-space plots for 8 training models 
# (response variable y = Ln ICP-MS; predictors x = Ln XRF-CS) 
# – 4 univariate (OLS, WLS, and weighted OLS, WLS), and 4 multivariate 
# Bayes GLM (Bayes), Random Forest (RF) and PLS (k-fold and LOO CV) based on 
# Ti individually for univariate models and 6-elements (Ca, Ti, Fe, Mn, Sr, Zr) 
# for multivariate models (n = 268).
# 
# Caption example for Prediction model 
# Predicted (XRF-CS in mg kg-1) vs observed (measured XRF-CS matched Ln[Ti/inc.]) 
# prediction model plots for Ti, applied to the ACE dataset to produce predicted 
# concentrations in mg kg-1 (see Methods for further details) (n = 268).

# v.1.8
# Global summary tables now include %based error summaries - not used for ranking
# ACE section element folder pred vs obs dimensions corrected 
# Corrected this calculation in Section 11 and 17 for RMSE values summay csv 
# RMSE_ppm_pc = ifelse(
# is.na(mean_obs_ppm) || mean_obs_ppm == 0, NA_real_, (RMSE_ppm * (RMSE_log_pc / 100))
# Added RMSE_ppm and RSMEP_ppm to AllModels & BestModels csv files 

# v. 1.7
# Created percentage based 95% upper and lower log space CI
# Rather than 1.96 * RMSE which creates a much larger CI 
# Applied % based CI calculation to XRF outputs in Section 8 and Section 15
# New Section 16 creates a RMSE pred error and log space %based pred error comparison
# Plotted these in age comparison plots rather than previous errors
# Element age comparison plots now only show red dashes where 
# % based prediction 95% CI > ICP-MS +2SD value - i/e. out of range
# light grey shading removed by assigning to _U95_pc_ppm because it was confusing 
# in section 8 - code for adding grey shading is still in there
# could be added back in in v1.8 perhaps but this is cleaner visually
# 95CI shading for RMSE based methods still exists in v1.5 & 1.6 outputs if needed 
# Added Sections 19, 20, 22 listed below
# Rmemoved check it work - now part of section 20 (TEster) writes check to console

# v. 1.6
# Various bug fixed 
# Including correcting reversed assignemnt of SD values for IPC and XRF in pred vs obs plots in Section 8 and 15

# v. 1.5
# Correlation diagnostic and element pairwise bi plots added
# pink errors bars removed from section 9
# all_y_equal plotted as log10 scale for easier visualiation / comparison
# red dots added where CI exceeds 2 x ICP_sd - i.e., noisy
# Added PLS(k) RMSEP 95% CI upper and lower added to Sites age comparison all_y_equal 
# This can be changed to any of the 8 models in section 10 
# Model 95% CI plotted can be changed by changing these three lines of code 
# pls_pred_col <- paste0("PLS_kfold_", el, "_Pred_ppm")
# pls_l95_col  <- paste0("PLS_kfold_", el, "_L95_ppm")
# pls_u95_col  <- paste0("PLS_kfold_", el, "_U95_ppm")
# has_pls_cols <- all(c(pls_pred_col, pls_l95_col, pls_u95_col) %in% names(df_el))
# Best to only plot 1 because too messy visually otherwise
# +/-SD error bars added to pred vs obs calibration models

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
# Section 1–7   setup, checks, helpers
# Section 8     calibration element loop
# Section 9     site-level depth plots
# Section 10    element-level age–conc comparisons
# Section 11    global calibration summary & ppm-scale calibration summary
# Section 12    finish message for calibration
# Section 13    ACE prediction
# Section 14    ACE diagnostic comparison vs ICP
# Section 15    ACE PREDICTION DIAGNOSTIC PLOTS
# Section 16    ACE SUMMARY PANELS (ORIGINAL + PERCENTAGE PI, SQUARE PANELS)
# Section 17    ACE GLOBAL SUMMARY (NA-safe)
# Section 18    ACE ppm-scale summary (NA-safe)
# Section 19    ACE DIAGNOSTICS DASHBOARD & VALIDATOR (New)
# Section 20    TESTER — CHECK ALL OUTPUTS EXIST
# Section 21    MULTI-ELEMENT CORRELATION DIAGNOSTICS
# Section 22    WRITES GLOBAL README FILE (placed in Regression_Multivariate)
# Workflow processing pipeline -------------------------------------

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
#   │   Pearson / Spearman CSVs
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
  "qqplotr", "future", "progressr", "fmsb"
))


# Element lists - not needed  -------------------------------------------------------

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
# Load libraries used  --------------------------------------------------
packages <-c(
  "broom", "performance", "see", "ggplot2", "dplyr", "purrr",
  "boot", "psych", "lmtest", "arm", "randomForest", "pls", "ggpubr", 
  "qqplotr", "ggsci", "car", "performance", "ggpmisc", "future", "progressr"
)
lapply(packages, library, character.only=TRUE)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Import datasets
# • ACE_dataset  --------
# - matched log ICP-MS & log_inc XRF-CS for calibration models
ACE_dataset <- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/ACE_Calibration/Input/ACE_subsample_xrf_icp_matched_log_inc.csv",
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

# • XRF_pred  ----------------------
# = XRF-CS log_inc dataset for new ppm predictions Section 1-12
main_elements_xrf <- c("Ca", "Fe", "Mn", "Sr", "Ti", "Zr")

XRF_pred <- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/ACE_Calibration/Input/ACE_ITRAX_qc_acf_log_inc.csv"
) %>%
  dplyr::select(Site, depth, SH20_age, dplyr::all_of(main_elements_xrf)) %>%
  dplyr::filter(Site %in% c("BI10","HER42PB","KER1","KER3","PB1")) %>%
  dplyr::mutate(across(all_of(main_elements_xrf), ~ ifelse(. <= -10, NA, .))) %>% 
  print()
str(XRF_pred)
names(ACE_dataset)

# • ICP_obs ---------------
# - matched ICP-MS measured - XRF-CS cps & errors for Section 13-21
ICP_obs<- readr::read_csv(
  "Papers_R/2024_DeVleeschouwer/ACE_Calibration/Input/ACE_subsample_icp_xrf_matched_cps.csv",
  name_repair = "minimal") %>%   # equivalent to check.names = FALSE 
  rename (SH20_age = SH20_mean_age) %>% 
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>% # Replace infinite values with NA
  #mutate(Site = as.factor(Site)) %>%  # Convert Site to factor
  print()
str(ICP_obs)
names(ICP_obs)

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Load Calibration model comparison function  --------------------

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
  
  all_base  <- file.path(out_dir, "All_models_ppm")
  best_base <- file.path(out_dir, "Best_models_ppm")
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
    plot.title   = element_text(size = 8, face = "bold"),
    panel.grid   = element_blank()      # removes all internal gridlines
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
  cols_models["OLS"]  <- "grey"
  cols_models["OLS(wt)"]  <- "firebrick1"
  cols_models["WLS"]  <- "darkslategrey"
  cols_models["WLS(wt)"]  <- "darkcyan"
  cols_models["Bayes"]  <- "gold2"
  cols_models["RF"]  <- "pink"
  cols_models["PLS(k)"] <- "salmon"
  cols_models["PLS(LOO)"] <- "salmon3"
  
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
    RMSE_ppm    = numeric(),
    RMSEP_ppm   = numeric(),
    AIC     = numeric(),
    BIC     = numeric()
  )
  
  preds_store <- list()
  
  # ===============================================================
  # 8. ELEMENT LOOP & Pred vs Obs PLOTS fro each element & model
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
    
    p_po_list <- list()
    p_po_list_noerr <- list()
    
    stopifnot(is.list(p_po_list), is.list(p_po_list_noerr))
    
    pb_elements$tick()
    
    y_col   <- paste0(el, "_ICP")
    x_col   <- el
    var_col <- paste0(el, var_suffix)
    sd_col  <- paste0(el, sd_suffix)
    
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
      RMSE_ppm  = numeric(),
      RMSEP_ppm = numeric(),
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
      
      # ------------------------------------------------------------------
      # Compute RMSE in PPM scale (for comparability with ACE)
      # ------------------------------------------------------------------
      
      # Convert observed + predicted back to ppm scale
      obs_ppm  <- exp(y)
      pred_ppm <- exp(pred)
      
      # Absolute RMSE in ppm
      RMSE_ppm <- tryCatch(
        sqrt(mean((obs_ppm - pred_ppm)^2, na.rm = TRUE)),
        error = function(e) NA_real_
      )
      
      # % RMSE relative to mean observed ppm
      mean_obs_ppm <- mean(obs_ppm, na.rm = TRUE)
      
      RMSE_ppm_pc <- ifelse(
        is.na(mean_obs_ppm) || mean_obs_ppm == 0,
        NA_real_,
        100 * RMSE_ppm / mean_obs_ppm
      )
      
      # ------------------------------------------------------------------
      # RMSEP in ppm scale (convert CV log predictions into ppm)
      # ------------------------------------------------------------------
      
      RMSEP_ppm <- tryCatch(
        {
          # Use same CV predictions used for RMSEP_val, but convert them:
          cv_log <- cv_predictions_for_model(      # You already have RMSEP_val fn
            mname,
            cv_formula,
            data,
            y_col,
            var_col,
            sd_col,
            best_ncomp = if (inherits(m, "mvr")) attr(m, "best_ncomp") else NA
          )
          
          if (is.null(cv_log)) {
            NA_real_
          } else {
            cv_ppm <- exp(cv_log$pred)
            obs_cv_ppm <- exp(cv_log$obs)
            sqrt(mean((obs_cv_ppm - cv_ppm)^2, na.rm = TRUE))
          }
        },
        error = function(e) NA_real_
      )
      
      AIC_val <- if (inherits(m, c("lm","glm","bayesglm")))
        tryCatch(AIC(m), error = function(e) NA_real_) else NA_real_
      
      BIC_val <- if (inherits(m, c("lm","glm","bayesglm")))
        tryCatch(BIC(m), error = function(e) NA_real_) else NA_real_
      
      ranking_df <- dplyr::add_row(
        ranking_df,
        Model       = mname,
        R2          = R2,
        RMSE        = RMSE,
        RMSEP       = RMSEP_val,
        RMSE_ppm    = RMSE_ppm,
        RMSEP_ppm   = RMSEP_ppm,
        AIC         = AIC_val,
        BIC         = BIC_val
      )
      
      global_summary <- dplyr::add_row(
        global_summary,
        Element = el,
        Model   = mname,
        R2      = R2,
        RMSE    = RMSE,
        RMSEP   = RMSEP_val,
        RMSE_ppm    = RMSE_ppm,
        RMSEP_ppm   = RMSEP_ppm,
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
      
      # -----------------------------------------------------------
      # Predictions for XRF_new  (log-space and ppm)
      # -----------------------------------------------------------
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
        
        base_nm <- paste0(mname, "_", el)
        
        # ppm-scale values (original PI)
        preds_df[[paste0(base_nm, "_Pred_ppm")]] <- exp(pred_new_log)
        preds_df[[paste0(base_nm, "_L95_ppm")]]  <- exp(lower_log)
        preds_df[[paste0(base_nm, "_U95_ppm")]]  <- exp(upper_log)
        
        # ---------------------------------------------------------
        # NEW: log-space % error AND ppm-scale % error intervals
        # ---------------------------------------------------------
        safe_pred_log <- ifelse(is.na(pred_new_log) | pred_new_log == 0, NA, pred_new_log)
        
        L95_log_pc <- (safe_pred_log - lower_log) / safe_pred_log * 100
        U95_log_pc <- (upper_log - safe_pred_log) / safe_pred_log * 100
        
        preds_df[[paste0(base_nm, "_L95_log_pc")]] <- L95_log_pc
        preds_df[[paste0(base_nm, "_U95_log_pc")]] <- U95_log_pc
        
        pred_ppm <- exp(pred_new_log)
        L95_pc_ppm <- pred_ppm - (pred_ppm * (L95_log_pc / 100))
        U95_pc_ppm <- pred_ppm + (pred_ppm * (U95_log_pc / 100))
        
        preds_df[[paste0(base_nm, "_L95_pc_ppm")]] <- L95_pc_ppm
        preds_df[[paste0(base_nm, "_U95_pc_ppm")]] <- U95_pc_ppm
      }
      
      # -----------------------------------------------------------
      # PLOTS: Pred vs Obs, Residuals, Influence, PLS RMSEP
      # -----------------------------------------------------------
      stats_label <- paste0(
        "R² = ", signif(R2, 3),
        "\nRMSE = ", signif(RMSE, 4),
        "\nRMSEP = ", signif(RMSEP_val, 4),
        "\nRMSE_ppm = ", signif(RMSE_ppm, 4)
      )
      
      df_pred <- data.frame(
        Observed  = y,
        Predicted = pred,
        Site      = data$Site
      )
      
      # -----------------------------------------------------------
      # Add ±SD error bars (Observed: ICP-MS SD)
      # -----------------------------------------------------------
      
      # ICP-MS SD → horizontal bars (Observed / x-axis)
      icp_sd_col <- paste0(el, "_ICP_sd")   # e.g. "Ti_ICP_sd"
      
      df_pred$Obs_low   <- df_pred$Observed
      df_pred$Obs_high  <- df_pred$Observed
      
      if (icp_sd_col %in% names(data)) {
        sd_icp <- data[[icp_sd_col]]
        df_pred$Obs_low  <- df_pred$Observed - sd_icp
        df_pred$Obs_high <- df_pred$Observed + sd_icp
      }
      
      # -----------------------------------------------------------
      # Add model CI (95%) for vertical prediction error bars
      # FULL CORRECTED CI DECISION LOGIC
      # -----------------------------------------------------------
      
      # Default: no CI
      df_pred$Pred_low  <- df_pred$Predicted
      df_pred$Pred_high <- df_pred$Predicted
      
      # Use prediction intervals ONLY for lm / glm (not bayesglm)
      if (inherits(m, c("lm", "glm")) && !inherits(m, "bayesglm")) {
        
        pi_obj <- tryCatch(
          predict(
            m,
            newdata  = data,
            interval = "prediction",
            level    = 0.95
          ),
          error = function(e) NULL
        )
        
        # Use only if prediction intervals are actually returned
        if (!is.null(pi_obj) &&
            is.matrix(pi_obj) &&
            ncol(pi_obj) >= 3) {
          
          # Columns are positionally: fit | lwr | upr
          df_pred$Pred_low  <- pi_obj[, 2]
          df_pred$Pred_high <- pi_obj[, 3]
        }
        
      } else {
        # bayesglm, RF, PLS → RMSE-based fallback
        df_pred$Pred_low  <- df_pred$Predicted - 1.96 * RMSE
        df_pred$Pred_high <- df_pred$Predicted + 1.96 * RMSE
      }
      
      # -----------------------------------------------------------
      # Convert CI bounds into centred asymmetric errors
      # -----------------------------------------------------------
      
      df_pred$Pred_err_low  <- df_pred$Predicted - df_pred$Pred_low
      df_pred$Pred_err_high <- df_pred$Pred_high - df_pred$Predicted
      
      # -----------------------------------------------------------
      # 1. Pred vs Obs plots WITH error bars (per element & model)
      # ----------------------------------------------------------- 
      
      # Axis limits must include full error bar extents
      x_min <- min(df_pred$Obs_low,  na.rm = TRUE)
      x_max <- max(df_pred$Obs_high, na.rm = TRUE)
      
      y_min <- min(df_pred$Predicted - df_pred$Pred_err_low,  na.rm = TRUE)
      y_max <- max(df_pred$Predicted + df_pred$Pred_err_high, na.rm = TRUE)
      
      # Force shared range so plots are square
      axis_min <- min(x_min, y_min)
      axis_max <- max(x_max, y_max)
      
      # Add small padding so points / error bars don't touch edges
      pad <- 0.03 * (axis_max - axis_min)
      axis_min_pad <- axis_min - pad
      axis_max_pad <- axis_max + pad
      
      # --- annotation coordinates: top-left INSIDE the padded limits ---
      ann_x <- axis_min_pad + 0.02 * (axis_max_pad - axis_min_pad)
      ann_y <- axis_max_pad - 0.02 * (axis_max_pad - axis_min_pad)
      
      p_po <- ggplot(df_pred, aes(Observed, Predicted, colour = Site)) +
        geom_segment(aes(
          x = Obs_low, xend = Obs_high,
          y = Predicted, yend = Predicted
        ),
        inherit.aes = FALSE,
        colour = "grey90", linewidth = 0.4) +
        
        geom_segment(aes(
          x = Observed, xend = Observed,
          y = Predicted - Pred_err_low,
          yend = Predicted + Pred_err_high
        ),
        inherit.aes = FALSE,
        colour = "grey90", linewidth = 0.4) +
        
        geom_point(size = 1.25) +
        ggsci::scale_color_jco(name = "Site") +
        geom_abline(slope = 1, intercept = 0) +
        scale_x_continuous(limits = c(axis_min_pad, axis_max_pad)) +
        scale_y_continuous(limits = c(axis_min_pad, axis_max_pad)) +
        annotate(
          "text",
          x = ann_x, y = ann_y,
          hjust = 0, vjust = 1,
          label = stats_label,
          size = 2.7
        ) +
        labs(
          x = expression("Observed ICP-MS [ln(mg kg"^{-1}*")]"),
          y = expression("Predicted ICP-MS [ln(mg kg"^{-1}*")]"),
          title = paste(el, mname, "- Calibration")
        ) +
        theme_bw() + theme_small +
        coord_fixed()
      
      ggsave(
        file.path(element_dir, paste0(el, "_", mname, "_PredVsObs_Calib_err.pdf")),
        p_po, width = 13, height = 9, units = "cm"
      )
      
      # STORE ERROR-BAR PLOT FOR SUMMARY PANELS  
      if (!exists("p_po_list")) p_po_list <- list()
      p_po_list[[mname]] <- p_po
      
      # -----------------------------------------------------------
      # 2. Pred vs Obs plots WITHOUT error bars (per element & model)
      # ----------------------------------------------------------- 
      
      # Axis limits must include full error bar extents
      x_min <- min(df_pred$Observed,  na.rm = TRUE)
      x_max <- max(df_pred$Observed, na.rm = TRUE)
      
      y_min <- min(df_pred$Predicted,  na.rm = TRUE)
      y_max <- max(df_pred$Predicted, na.rm = TRUE)
      
      # Force shared range so plots are square
      axis_min <- min(x_min, y_min)
      axis_max <- max(x_max, y_max)
      
      # Add small padding so points / error bars don't touch edges
      pad <- 0.03 * (axis_max - axis_min)
      axis_min_pad <- axis_min - pad
      axis_max_pad <- axis_max + pad
      
      # --- annotation coordinates: top-left INSIDE the padded limits ---
      ann_x <- axis_min_pad + 0.02 * (axis_max_pad - axis_min_pad)
      ann_y <- axis_max_pad - 0.02 * (axis_max_pad - axis_min_pad)
      
      # Square aspect
      p_noerr <- ggplot(df_pred, aes(Observed, Predicted, colour = Site)) +
        geom_point(size = 1.25) +
        ggsci::scale_color_jco(name = "Site") +
        geom_abline(slope = 1, intercept = 0) +
        scale_x_continuous(limits = c(axis_min_pad, axis_max_pad)) +
        scale_y_continuous(limits = c(axis_min_pad, axis_max_pad)) +
        annotate(
          "text",
          x = ann_x, y = ann_y,
          hjust = 0, vjust = 1,
          label = stats_label,
          size = 2.7
        )+
        labs(
          x = expression("Observed ICP-MS [Ln(mg kg"^{-1}*")]"),
          y = expression("Predicted [XRF-CS Ln(el/inc)]"),
          title = paste(el, mname, "- Pred vs Obs Training")
        ) +
        theme_bw() + theme_small +
        coord_fixed()   # <-- square panel
      
      # Save new version
      ggsave(
        file.path(element_dir, paste0(el, "_", mname, "_PredVsObs_no_err.pdf")),
        p_noerr,
        width = 12, height = 12, units = "cm"
      )
      
      # Store no-error plot for summary panel in Patch B
      if (!exists("p_po_list_noerr")) p_po_list_noerr <- list()
      p_po_list_noerr[[mname]] <- p_noerr
      
      
      # -----------------------------------------------------------      
      # Residuals
      # -----------------------------------------------------------
      res_df   <- data.frame(resid = resid_vals)
      hist_obj <- hist(resid_vals, plot = FALSE)
      
      p_res <- ggplot(res_df, aes(resid)) +
        geom_histogram(bins = 30, fill = "grey80") +
        geom_density(colour = "red") +
        annotate(
          "text",
          x = -Inf, y = Inf,
          hjust = -0.05, vjust = 1.1,
          label = stats_label,
          size = 2.7
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
    # 8.4 MODEL RANKING & DIAGNOSTIC OUTPUT (PATCHED)
    # ===============================================================
    if (nrow(ranking_df) > 0) {
      
      # 1. Compute ranking ONCE
      ranking_df <- ranking_df %>%
        dplyr::arrange(dplyr::desc(R2), RMSEP, RMSE) %>%
        dplyr::mutate(Rank = dplyr::row_number()) %>%
        dplyr::select(Rank, dplyr::everything())
      
      # 2. Add explanations (Rank preserved)
      ranking_expl_df <- ranking_df %>%
        dplyr::mutate(
          Rank_explanation  = expl_rank,
          R2_explanation    = expl_r2,
          RMSE_explanation  = expl_rmse,
          RMSEP_explanation = expl_rmsep
        ) %>%
        dplyr::select(Rank, dplyr::everything())
      
      # 3. Write CSV WITH Rank included
      write.csv(
        ranking_expl_df,
        file.path(element_dir, paste0(el, "_Model_Ranking.csv")),
        row.names = FALSE
      )
    }
    
# 8.4.1 DIAGNOSTICS SUMMARY OUTPUT (PATCH — REINSTATED)
    
    # Ensure diag_df always exists and has rows
    if (!exists("diag_df")) {
      warning("diag_df missing for element ", el, " — creating empty diagnostics table")
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
    }
    
    # Write diagnostics summary CSV per element
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
      "Additional uncertainty columns:",
      " *_L95_log_pc, *_U95_log_pc : % error in log-space.",
      " *_L95_pc_ppm, *_U95_pc_ppm : % error bounds in ppm.",
      "",
      "Models: OLS, WLS, OLS(wt), WLS(wt), Bayes, RF, PLS(LOO), PLS(k).",
      "Ranking: highest R², then lowest RMSEP, then lowest RMSE.",
      "Plots show predicted concentrations (Y-axis) versus observed ICP-MS values (X-axis); horizontal error bars represent ICP analytical uncertainty, vertical error bars represent model prediction uncertainty (95% CI).",
      "Model robustness was evaluated using a composite signal-to-noise framework incorporating prediction uncertainty, downcore smoothness, and calibration performance. Signal-to-noise ratios were calculated from prediction confidence intervals and profile roughness, scaled within each element, and combined with R² into a weighted robustness score. Models exhibiting unstable downcore behaviour were excluded from selection. For each element, the highest-ranked stable model was designated as the production model and used for subsequent interpretation.
      "
    )
    
    writeLines(readme_txt, file.path(element_dir, "README.txt"))
    
    # Check these plots are stored correctly
    message("Stored plots for ", el, ":")
    print(names(p_po_list))
    print(names(p_po_list_noerr))
    
    # ===============================================================
    # 8.7  SUMMARY PANELS (2×2) OF ERROR (SQUARE) ON 2 A4 PAGES
    # (UNCHANGED LAYOUT; USING p_po_list CALIBRATION PLOTS)
    # ===============================================================
    
    linear_models   <- c("OLS", "OLS_weighted", "WLS", "WLS_weighted")
    multivar_models <- c("Bayes", "RF", "PLS_kfold", "PLS_LOO")
    
    get_plot_or_blank <- function(model_name) {
      if (!is.null(p_po_list[[model_name]])) {
        p_po_list[[model_name]]
      } else {
        ggplot() + theme_void() + coord_fixed()
      }
    }
    
    # -------------------------
    # PAGE 1 — Linear models
    # -------------------------
    p1_lin <- get_plot_or_blank("OLS")
    p2_lin <- get_plot_or_blank("OLS_weighted")
    p3_lin <- get_plot_or_blank("WLS")
    p4_lin <- get_plot_or_blank("WLS_weighted")
    
    page1 <- (
      (p1_lin | p2_lin) /
        (p3_lin | p4_lin)
    ) +
      patchwork::plot_layout(guides = "collect") &
      theme(
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.box       = "horizontal"
      ) &
      guides(colour = guide_legend(nrow = 1))
    
    # -------------------------
    # PAGE 2 — Multivariate models
    # -------------------------
    p1_mv <- get_plot_or_blank("Bayes")
    p2_mv <- get_plot_or_blank("RF")
    p3_mv <- get_plot_or_blank("PLS_kfold")
    p4_mv <- get_plot_or_blank("PLS_LOO")
    
    page2 <- (
      (p1_mv | p2_mv) /
        (p3_mv | p4_mv)
    ) +
      patchwork::plot_layout(guides = "collect") &
      theme(
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.box       = "horizontal"
      ) &
      guides(colour = guide_legend(nrow = 1))
    
    pdf(
      file   = file.path(element_dir, paste0(el, "_PredObs_Model_Summary_err_2pp.pdf")),
      width  = 21 / 2.54,
      height = 29.7 / 2.54
    )
    print(page1)
    print(page2)
    dev.off()
    
    # ===============================================================
    # 8.9  SUMMARY PANEL — PORTRAIT 2×4 OF ERROR PLOTS (SQUARE)
    # ===============================================================
    
    all_models_order <- c(
      "OLS",
      "OLS_weighted",
      "WLS",
      "WLS_weighted",
      "Bayes",
      "RF",
      "PLS_kfold",
      "PLS_LOO"
    )
    
    make_square <- function(p) {
      p + theme(
        aspect.ratio        = 1,
        plot.title.position = "panel"
      )
    }
    
    get_plot_or_blank <- function(model_name) {
      if (!is.null(p_po_list[[model_name]])) {
        make_square(p_po_list[[model_name]])
      } else {
        make_square(ggplot() + theme_void())
      }
    }
    
    all_plots <- lapply(all_models_order, get_plot_or_blank)
    
    summary_panel_all <- wrap_plots(
      all_plots,
      ncol = 2,
      nrow = 4,
      guides = "collect"
    ) &
      theme(
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.box       = "horizontal"
      ) &
      guides(colour = guide_legend(nrow = 1))
    
    ggsave(
      filename = paste0(el, "_PredObs_Model_Summary_err.pdf"),
      path     = element_dir,
      plot     = summary_panel_all,
      width    = 21,
      height   = 29.7,
      units    = "cm"
    )
    
    # ===============================================================
    # NEW: 8.10  SUMMARY PANEL — 2×4 GRID OF NO-ERROR PLOTS (SQUARE)
    # ===============================================================
    
    all_models_order <- c(
      "OLS",
      "OLS_weighted",
      "WLS",
      "WLS_weighted",
      "Bayes",
      "RF",
      "PLS_kfold",
      "PLS_LOO"
    )
    
    get_plot_noerr <- function(model_name) {
      if (!is.null(p_po_list_noerr[[model_name]])) {
        p_po_list_noerr[[model_name]]
      } else {
        ggplot() + theme_void() + coord_fixed()
      }
    }
    
    noerr_plots <- lapply(all_models_order, get_plot_noerr)
    
    summary_panel_noerr <- wrap_plots(
      noerr_plots,
      ncol = 2, nrow = 4,
      guides = "collect"
    ) &
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal"
      ) &
      guides(colour = guide_legend(nrow = 1))
    
    # Save square A4-style portrait PDF
    ggsave(
      filename = paste0(el, "_PredObs_Summary_NoErr.pdf"),
      path     = element_dir,
      plot     = summary_panel_noerr,
      width    = 21,      # A4 width in cm
      height   = 29.7,    # A4 height in cm
      units    = "cm"
    )
    
  } # end element loop
  
  # ===============================================================
  # 9. SITE-LEVEL A4 PLOTS (DEPTH + AGE)
  # ===============================================================
  
  dynamic_layout <- function(n_panels, max_cols = 3) {
    ncol <- min(max_cols, max(1, ceiling(sqrt(n_panels))))
    nrow <- ceiling(n_panels / ncol)
    list(nrow = nrow, ncol = ncol)
  }
  
  if (!is.null(XRF_new) && length(preds_store) > 0) {
    
    # Ensure ICP-MS (Observed) has a defined colour (blue)
    if (!"ICP-MS (Observed)" %in% names(cols_all)) {
      cols_all["ICP-MS (Observed)"] <- "blue"
    }
    
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
            
            # ---- All models depth plot ----
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
              scale_color_manual(values = cols_all) +
              theme_bw() + theme_small
            
            if (!is.null(ace_overlay)) {
              p_depth_all <- p_depth_all +
                geom_path(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth, colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth, colour = "ICP-MS (Observed)"),
                  size = 1.5,
                  shape = 21,
                  fill = "white",
                  inherit.aes = FALSE
                )
            }
            
            # ---- Multivariate depth plot ----
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
              p_depth_best <- p_depth_best +
                geom_path(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth, colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay,
                  aes(x = ICP_ppm, y = depth, colour = "ICP-MS (Observed)"),
                  size = 1.5,
                  shape = 21,
                  fill = "white",
                  inherit.aes = FALSE
                )
            }
            
            # ✅ store plots for wrapping & saving
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
            
            # ---- All models age plot ----
            p_age_all <- ggplot(df_long,
                                aes(x = Pred_ppm, y = SH20_age,
                                    colour = Model, group = Model)) +
              geom_path(linewidth = 0.4) +
              labs(
                x = bquote(.(el)~"(mg"~kg^{-1}*")"),
                y = "Age (cal a BP)",
                title = paste("Site", s, "-", el, "(age, all models)")
              ) +
              scale_color_manual(values = cols_all) +
              theme_bw() + theme_small
            
            if (!is.null(ace_overlay_age)) {
              p_age_all <- p_age_all +
                geom_path(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age, colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age, colour = "ICP-MS (Observed)"),
                  size = 1.5,
                  shape = 21,
                  fill = "white",
                  inherit.aes = FALSE
                )
            }
            
            # ---- Multivariate age plot ----
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
              p_age_best <- p_age_best +
                geom_path(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age, colour = "ICP-MS (Observed)"),
                  linewidth = 0.6,
                  inherit.aes = FALSE
                ) +
                geom_point(
                  data = ace_overlay_age,
                  aes(x = ICP_ppm, y = SH20_age, colour = "ICP-MS (Observed)"),
                  size = 1.5,
                  shape = 21,
                  fill = "white",
                  inherit.aes = FALSE
                )
            }
            
            # ✅ store plots
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
    
    # Ensure ICP colour exists (blue)
    if (!"ICP-MS (Observed)" %in% names(cols_all)) {
      cols_all["ICP-MS (Observed)"] <- "blue"
    }
    
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
        dplyr::filter(!is.na(SH20_age), !is.na(.data[[obs_col]])) %>%
        dplyr::mutate(
          ICP_ppm    = .data[[obs_col]],
          ICP_sd_ppm = if (has_sd_el) .data[[sd_col]] else NA_real_
        )
      
      # Global y-range for log-equal
      y_pred_vals <- df_long_all$Pred_ppm
      y_obs_vals  <- ace_all$ICP_ppm
      y_all_vals  <- c(y_pred_vals, y_obs_vals)
      y_all_vals  <- y_all_vals[is.finite(y_all_vals)]
      
      if (!length(y_all_vals) || min(y_all_vals) <= 0) {
        y_min_el <- NA_real_; y_max_el <- NA_real_
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
      pls_l95_col  <- paste0("PLS_kfold_", el, "_L95_pc_ppm")
      pls_u95_col  <- paste0("PLS_kfold_", el, "_U95_pc_ppm")
      has_pls_cols <- all(c(pls_pred_col, pls_l95_col, pls_u95_col) %in% names(df_el))
      
      plot_list_all          <- list()
      plot_list_best         <- list()
      plot_list_all_y_equal  <- list()
      plot_list_best_y_equal <- list()
      
      for (s in site_levels_all) {
        
        df_pred_site_all  <- df_long_all  %>% dplyr::filter(Site == s)
        df_pred_site_best <- df_long_best %>% dplyr::filter(Site == s)
        ace_site          <- ace_all      %>% dplyr::filter(Site == s)
        
        # Sample counts
        df_raw_xrf_site <- df_el   %>% dplyr::filter(Site == s)
        icp_site_raw    <- ICP_ppm %>% dplyr::filter(Site == s)
        
        n_xrf <- sum(!is.na(df_raw_xrf_site$SH20_age))
        n_icp <- sum(!is.na(icp_site_raw[[obs_col]]) &
                       !is.na(icp_site_raw$SH20_age))
        n_label <- paste0("n = ", n_icp, " (ICP), ", n_xrf, " (XRF)")
        
        if (nrow(ace_site) > 0) {
          site_median <- median(ace_site$ICP_ppm, na.rm = TRUE)
          icp_upper_max <- max(ace_site$ICP_ppm + ace_site$ICP_sd_ppm, na.rm = TRUE)
        } else {
          site_median <- NA_real_; icp_upper_max <- NA_real_
        }
        
        # PLS ribbons
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
        } else df_pls_site <- NULL
        
        # Exceedance
        ast_df <- NULL
        if (!is.null(df_pls_site) && !is.na(icp_upper_max)) {
          thresh <- icp_upper_max + 2 * max(ace_site$ICP_sd_ppm, na.rm = TRUE)
          exceed_idx <- with(df_pls_site,
                             !is.na(U95_ppm) & U95_ppm > thresh)
          if (any(exceed_idx)) ast_df <- df_pls_site[exceed_idx, , drop = FALSE]
        }
        
        y_top_free <- suppressWarnings(
          max(c(df_pred_site_all$Pred_ppm, ace_site$ICP_ppm), na.rm = TRUE)
        )
        if (!is.finite(y_top_free)) y_top_free <- NA_real_
        
        ##################################################################
        # 10.1 ALL MODELS (free y)
        ##################################################################
        
        if (nrow(df_pred_site_all) == 0 && nrow(ace_site) == 0) {
          p_all <- ggplot() + theme_void()
        } else {
          
          p_all <- ggplot(df_pred_site_all,
                          aes(x = SH20_age, y = Pred_ppm,
                              colour = Model, group = Model)) +
            theme_bw() + theme_small +
            scale_color_manual(values = cols_all) +
            coord_cartesian(xlim = c(age_min, age_max)) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(all models)")
            )
          
          # ---- ribbon first (if exists) ----
          if (!is.null(df_pls_site)) {
            p_all <- p_all +
              geom_ribbon(
                data = df_pls_site,
                aes(x = SH20_age, ymin = L95_ppm, ymax = U95_ppm),
                inherit.aes = FALSE,
                fill = "grey80", alpha = 0.5
              )
          }
          
          # ---- model lines ----
          p_all <- p_all + geom_path(linewidth = 0.4)
          
          # ---- ICP ----
          if (nrow(ace_site) > 0) {
            p_all <- p_all +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white",
                stroke = 0.4, size = 1.4, inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              )
          }
          
          # ---- exceedance markers ----
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_all <- p_all +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red",
                size = 5, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          # ---- annotation ----
          p_all <- p_all +
            annotate("text", x = Inf, y = Inf,
                     label = n_label,
                     hjust = 1.05, vjust = 1.1, size = 2.7)
        }
        
        ##################################################################
        # 10.2 ALL MODELS — LOG10 Y (ribbon behind models)
        ##################################################################
        
        if (inherits(p_all$theme, "theme_void") ||
            is.na(y_min_el) || is.na(y_max_el)) {
          p_all_y <- p_all
        } else {
          
          p_all_y <- ggplot(df_pred_site_all,
                            aes(x = SH20_age, y = Pred_ppm,
                                colour = Model, group = Model)) +
            theme_bw() + theme_small +
            scale_color_manual(values = cols_all) +
            coord_cartesian(xlim = c(age_min, age_max)) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(all models, y-equal)")
            )
          
          # ---- ribbon FIRST ----
          if (!is.null(df_pls_site)) {
            p_all_y <- p_all_y +
              geom_ribbon(
                data = df_pls_site,
                aes(x = SH20_age, ymin = L95_ppm, ymax = U95_ppm),
                inherit.aes = FALSE,
                fill = "grey80", alpha = 0.5
              )
          }
          
          # ---- model lines ----
          p_all_y <- p_all_y + geom_path(linewidth = 0.4)
          
          # ---- ICP ----
          if (nrow(ace_site) > 0) {
            p_all_y <- p_all_y +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white",
                stroke = 0.4, size = 1.4, inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              )
          }
          
          # ---- exceedance markers ----
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_all_y <- p_all_y +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red",
                size = 3, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          # ---- annotation + log scale ----
          p_all_y <- p_all_y +
            annotate("text", x = Inf, y = Inf,
                     label = n_label,
                     hjust = 1, vjust = 1, size = 2.7) +
            scale_y_log10(
              limits = c(y_min_el, y_max_el),
              breaks = y_log_breaks,
              labels = y_log_labels,
              oob = scales::oob_squish
            )
        }
        
        ##################################################################
        # 10.3 MULTIVARIATE — free y (ribbon behind)
        ##################################################################
        
        if (nrow(df_pred_site_best) == 0 && nrow(ace_site) == 0) {
          p_best <- ggplot() + theme_void()
        } else {
          
          p_best <- ggplot(df_pred_site_best,
                           aes(x = SH20_age, y = Pred_ppm,
                               colour = Model, group = Model)) +
            theme_bw() + theme_small +
            scale_color_manual(values = cols_all) +
            coord_cartesian(xlim = c(age_min, age_max)) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(multivariate)")
            )
          
          # -- ribbon first --
          if (!is.null(df_pls_site)) {
            p_best <- p_best +
              geom_ribbon(
                data = df_pls_site,
                aes(x = SH20_age, ymin = L95_ppm, ymax = U95_ppm),
                inherit.aes = FALSE,
                fill = "grey80", alpha = 0.5
              )
          }
          
          # -- model lines --
          p_best <- p_best + geom_path(linewidth = 0.4)
          
          # -- ICP --
          if (nrow(ace_site) > 0) {
            p_best <- p_best +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white",
                stroke = 0.4, size = 1.4, inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              )
          }
          
          # -- exceedance --
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_best <- p_best +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red",
                size = 3, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          # -- annotation --
          p_best <- p_best +
            annotate("text", x = age_max, y = y_top_free,
                     label = n_label, hjust = 1.05, vjust = 1.1, size = 2.7)
        }
        
        ##################################################################
        # 10.4 MULTIVARIATE — y-equal log10 (ribbon behind)
        ##################################################################
        
        if (inherits(p_best$theme, "theme_void") ||
            is.na(y_min_el) || is.na(y_max_el)) {
          p_best_y <- p_best
        } else {
          
          p_best_y <- ggplot(df_pred_site_best,
                             aes(x = SH20_age, y = Pred_ppm,
                                 colour = Model, group = Model)) +
            theme_bw() + theme_small +
            scale_color_manual(values = cols_all) +
            coord_cartesian(xlim = c(age_min, age_max)) +
            labs(
              x = "Age (cal a BP)",
              y = bquote(.(el)~"(mg"~kg^{-1}*")"),
              title = paste("Site", s, "(multivariate, y-equal)")
            )
          
          # -- ribbon first --
          if (!is.null(df_pls_site)) {
            p_best_y <- p_best_y +
              geom_ribbon(
                data = df_pls_site,
                aes(x = SH20_age, ymin = L95_ppm, ymax = U95_ppm),
                inherit.aes = FALSE,
                fill = "grey80", alpha = 0.5
              )
          }
          
          # -- model lines --
          p_best_y <- p_best_y + geom_path(linewidth = 0.4)
          
          # -- ICP --
          if (nrow(ace_site) > 0) {
            p_best_y <- p_best_y +
              geom_path(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                linewidth = 0.6, inherit.aes = FALSE
              ) +
              geom_point(
                data = ace_site,
                aes(x = SH20_age, y = ICP_ppm,
                    colour = "ICP-MS (Observed)"),
                shape = 21, fill = "white",
                stroke = 0.4, size = 1.4,
                inherit.aes = FALSE
              ) +
              geom_hline(
                yintercept = site_median, linetype = "dashed",
                colour = "grey40", linewidth = 0.4
              )
          }
          
          # -- exceedance --
          if (!is.null(ast_df) && !is.na(y_top_free)) {
            p_best_y <- p_best_y +
              geom_text(
                data = transform(ast_df, y_ast = y_top_free),
                aes(x = SH20_age, y = y_ast),
                label = ".", colour = "red",
                size = 3, vjust = -0.2,
                inherit.aes = FALSE
              )
          }
          
          # -- annotate + log scale --
          p_best_y <- p_best_y +
            annotate("text", x = Inf, y = Inf,
                     label = n_label,
                     hjust = 1.05, vjust = 1.1, size = 2.7) +
            scale_y_log10(
              limits = c(y_min_el, y_max_el),
              breaks = y_log_breaks,
              labels = y_log_labels,
              oob = scales::oob_squish
            )
        }
        
        ##################################################################
        # Store all four plots
        ##################################################################
        
        plot_list_all[[s]]          <- p_all
        plot_list_all_y_equal[[s]]  <- p_all_y
        plot_list_best[[s]]         <- p_best
        plot_list_best_y_equal[[s]] <- p_best_y
      }
      
      ##################################################################
      # Output assemblers
      ##################################################################
      
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
      
      ##################################################################
      # Save outputs
      ##################################################################
      
      combined_all <- combine_1col(plot_list_all)
      if (!is.null(combined_all)) {
        ggsave(
          file.path(elements_dir_all,
                    paste0(el, "_model_age_comparison_all.pdf")),
          combined_all, width = 21, height = 29.7, units = "cm"
        )
      }
      
      combined_all_y <- combine_1col(plot_list_all_y_equal)
      if (!is.null(combined_all_y)) {
        ggsave(
          file.path(elements_dir_all,
                    paste0(el, "_model_age_comparison_all_y_equal.pdf")),
          combined_all_y, width = 21, height = 29.7, units = "cm"
        )
      }
      
      combined_best <- combine_1col(plot_list_best)
      if (!is.null(combined_best)) {
        ggsave(
          file.path(elements_dir_best,
                    paste0(el, "_model_age_comparison_multivariate.pdf")),
          combined_best, width = 21, height = 29.7, units = "cm"
        )
      }
      
      combined_best_y <- combine_1col(plot_list_best_y_equal)
      if (!is.null(combined_best_y)) {
        ggsave(
          file.path(elements_dir_best,
                    paste0(el, "_model_age_comparison_multivariate_y_equal.pdf")),
          combined_best_y, width = 21, height = 29.7, units = "cm"
        )
      }
      
    } # element loop
  } # main if
  # ===============================================================
  # 11. GLOBAL SUMMARY OF ALL ELEMENTS AND ALL MODELS (PATCHED)
  # ===============================================================
  
  message("=== Building global calibration summary (patched Section 11) ===")
  
  # global_summary is filled inside Section 8 model loop.
  # Now we compute ppm-scale RMSE and percentage metrics for all rows.
  
  if (!exists("global_summary") || nrow(global_summary) == 0) {
    warning("global_summary does not exist or is empty.")
  } else {
    
    # Expand global summary with ppm-scale RMSE values.
    # NOTE: We must recompute RMSE_ppm, RMSE_log_pc, RMSE_ppm_pc
    # because Section 8 had log-scale RMSE only.
    
    global_summary <- global_summary %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        # Recompute ppm RMSE from cached preds_store
        mean_obs_ppm = {
          el <- Element
          y_col <- paste0(el, "_ICP")
          if (y_col %in% names(data))
            mean(exp(data[[y_col]]), na.rm = TRUE)
          else NA_real_
        },
        
        RMSE_log_pc = {
          # log-space percentage error
          mean_pred_log <- {
            el <- Element
            mname <- Model
            m <- model_store[[el]][[mname]]
            if (is.null(m)) NA_real_ else {
              if (inherits(m, "mvr"))
                mean(as.numeric(predict(m, ncomp = attr(m, "best_ncomp"))), na.rm = TRUE)
              else
                mean(as.numeric(predict(m)), na.rm = TRUE)
            }
          }
          if (is.na(mean_pred_log) || abs(mean_pred_log) < 1e-12) NA_real_
          else (RMSE / abs(mean_pred_log)) * 100
        },
        
        RMSE_ppm_pc = {
          if (is.na(mean_obs_ppm) || mean_obs_ppm == 0 || is.na(RMSE_ppm))
            NA_real_
          else
            (RMSE_ppm * (RMSE_log_pc / 100))
        },
        
        RMSEP_ppm_pc = {
          if (is.na(mean_obs_ppm) || mean_obs_ppm == 0 || is.na(RMSEP_ppm))
            NA_real_
          else
            (RMSEP_ppm / mean_obs_ppm) * 100
        }
      ) %>%
      dplyr::ungroup()
    
    # -----------------------------------------------------------
    # NEW: Compute per-element model ranking (matches Section 17)
    # -----------------------------------------------------------
    
    global_summary <- global_summary %>%
      dplyr::group_by(Element) %>%
      dplyr::arrange(
        dplyr::desc(R2),
        RMSEP,
        RMSE,
        .by_group = TRUE
      ) %>%
      dplyr::mutate(Rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    # Write updated global CSV
    global_file <- file.path(all_base, "AllElements_ModelSummary.csv")
    write.csv(global_summary, global_file, row.names = FALSE)
    
    # Determine best model per element (Rank = 1)
    best_models <- global_summary %>%
      dplyr::group_by(Element) %>%
      dplyr::arrange(desc(R2), RMSE, RMSEP) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    
    best_file <- file.path(all_base, "BestModels_PerElement.csv")
    write.csv(best_models, best_file, row.names = FALSE)
    
    message("=== Section 11 summaries written successfully ===")
  }
  
  # ===============================================================  
  # SECTION 11A: ROBUSTNESS of XRF_pred predictions file
  # XRF_pred SNR, robustness, stability flags, production model
  # Ranked by Robustness score using Global Rank Priority:
  #   1) Robustness_score (overall) - 0.4x SNR & 0.3x Smoothness & 0.3x R2
  #   2) R2 (overall)
  #   3) RMSE (log-space) (overall)
  #   4) RMSEP (RMSE_ppm) (overall)
  # Flagged as unstable using Roughness scores <0.7 sd(diff(pred) normalised
  # ===============================================================
  
  # CREATE ROBUSTNESS OUTPUT DIRECTORY
  
  robust_root <- file.path(all_base, "Robustness")
  
  if (!dir.exists(robust_root)) {
    dir.create(robust_root, recursive = TRUE, showWarnings = FALSE)
  }
  
  stopifnot(dir.exists(robust_root))
  
  # Set up SNR calucaltions 
  
  scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(diff(rng)) || diff(rng) == 0)
      return(rep(NA_real_, length(x)))
    (x - rng[1]) / diff(rng)
  }
  
  snr_rows <- list()
  
  # 1. Compute SNR + stability metrics (XRF_pred ONLY)
  
  for (el in elements) {
    
    df <- preds_store[[el]]
    if (is.null(df)) next
    
    for (mname in names(model_store[[el]])) {
      
      pred_col <- paste0(mname, "_", el, "_Pred_ppm")
      l_col    <- paste0(mname, "_", el, "_L95_ppm")
      u_col    <- paste0(mname, "_", el, "_U95_ppm")
      
      if (!all(c(pred_col, l_col, u_col) %in% names(df))) next
      
      pred <- df[[pred_col]]
      l95  <- df[[l_col]]
      u95  <- df[[u_col]]
      
      ok <- !is.na(pred) & !is.na(l95) & !is.na(u95)
      if (sum(ok) < 6) next
      
      pred <- pred[ok]
      ci_w <- u95[ok] - l95[ok]
      
      rmse_ppm <- global_summary %>%
        dplyr::filter(Element == el, Model == mname) %>%
        dplyr::pull(RMSE_ppm)
      
      roughness <- sd(diff(pred), na.rm = TRUE)
      
      snr_rows[[length(snr_rows) + 1]] <- data.frame(
        Element        = el,
        Model          = mname,
        Mean_ppm       = mean(pred, na.rm = TRUE),
        RMSE_ppm       = rmse_ppm,
        Mean_CI_width  = mean(ci_w, na.rm = TRUE),
        SD_pred        = sd(pred, na.rm = TRUE),
        Roughness      = roughness,
        SNR_model      = mean(pred, na.rm = TRUE) / rmse_ppm,
        SNR_CI         = sd(pred, na.rm = TRUE) / mean(ci_w, na.rm = TRUE),
        SNR_smooth     = sd(pred, na.rm = TRUE) / roughness
      )
    }
  }
  
  snr_df <- dplyr::bind_rows(snr_rows)
  
  # 2. Merge SNR metrics into global_summary
  
  global_summary <- global_summary %>%
    dplyr::left_join(
      snr_df %>%
        dplyr::select(
          Element, Model,
          Mean_ppm, Mean_CI_width,
          SNR_model, SNR_CI, SNR_smooth,
          Roughness
        ),
      by = c("Element", "Model")
    )
  
  # 3. Robustness score (scaled SNR + R²)
  
  global_summary <- global_summary %>%
    dplyr::group_by(Element) %>%
    dplyr::mutate(
      SNR_CI_s      = scale01(SNR_CI),
      SNR_smooth_s = scale01(SNR_smooth),
      R2_s         = scale01(R2),
      
      Robustness_Score =
        0.4 * SNR_CI_s +
        0.3 * SNR_smooth_s +
        0.3 * R2_s
    ) %>%
    dplyr::arrange(dplyr::desc(Robustness_Score), .by_group = TRUE) %>%
    dplyr::mutate(
      Robust_Rank = dplyr::row_number()
    ) %>%
    dplyr::ungroup()
  
  # PATCH 3.1: GLOBAL_RANK using same priority as Section 11
  global_summary <- global_summary %>%
    dplyr::group_by(Element) %>%
    dplyr::arrange(
      dplyr::desc(Robustness_Score),
      dplyr::desc(R2),
      RMSE_ppm,
      RMSEP,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      Global_Rank = dplyr::row_number()
    ) %>%
    dplyr::ungroup()
  
  # 4. Flag unstable downcore behaviour
  
  global_summary <- global_summary %>%
    dplyr::group_by(Element) %>%
    dplyr::mutate(
      Roughness_norm = scale01(Roughness),
      Unstable_flag  = ifelse(
        is.na(Roughness_norm), NA,
        Roughness_norm > 0.7
      )
    ) %>%
    dplyr::ungroup()
  
  # 5. Auto-select ONE production model per element
  
  production_models <- global_summary %>%
    dplyr::filter(
      Robust_Rank == 1,
      Unstable_flag == FALSE | is.na(Unstable_flag)
    ) %>%
    dplyr::group_by(Element) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Production_Model = Model
    )
  
  global_summary <- global_summary %>%
    dplyr::left_join(
      production_models %>%
        dplyr::select(Element, Production_Model),
      by = "Element"
    ) %>%
    dplyr::mutate(
      Is_Production_Model = Model == Production_Model
    )
  
  # 6. Robustness class (for colour-coding)
  
  global_summary <- global_summary %>%
    dplyr::mutate(
      Robustness_Class = dplyr::case_when(
        Global_Rank == 1 & !Unstable_flag ~ "Preferred 1",
        Global_Rank <= 2 & !Unstable_flag ~ "Acceptable 2",
        Global_Rank <= 3 & !Unstable_flag ~ "Acceptable 3",
        Global_Rank <= 4 & !Unstable_flag ~ "Acceptable 4",
        TRUE                              ~ "Unstable"
      )
    )
  
  # 7. Write updated global summary
  
  write.csv(
    global_summary,
    file.path(all_base, "AllElements_ModelSummary.csv"),
    row.names = FALSE
  )
  
  # 8. Supplementary audit table
  
  
  # Ensure Robustness output directory exists
  robust_root <- file.path(all_base, "Robustness")
  if (!dir.exists(robust_root)) {
    dir.create(robust_root, recursive = TRUE, showWarnings = FALSE)
  }
  
  audit_table <- global_summary %>%
    dplyr::select(
      Element, Model,
      R2, RMSE_ppm,
      Mean_ppm, Mean_CI_width,
      SNR_model, SNR_CI, SNR_smooth,
      Roughness,
      Robustness_Score, Robust_Rank,
      Robustness_Class,
      Unstable_flag,
      Is_Production_Model
    ) %>%
    dplyr::arrange(Element, Robust_Rank)
  
  write.csv(
    audit_table,
    file.path(robust_root, "AllElementsModel_Robustness_Stability.csv"),
    row.names = FALSE
  )
  
  # 9. CONFIDENCE FLAGS — UNCERTAINTY VS SIGNAL DOMINANCE
  
  global_summary <- global_summary %>%
    dplyr::mutate(
      Confidence_Flag = dplyr::case_when(
        
        # Strong model: high robustness, not unstable
        Robustness_Score >= 0.75 &
          !Unstable_flag ~
          "High confidence (signal-dominated)",
        
        # Uncertainty-dominated: wide CI relative to variance
        SNR_CI_s < 0.4 &
          Robustness_Score < 0.6 ~
          "Low confidence (uncertainty-dominated)",
        
        # Noise-dominated: rough downcore behaviour
        Unstable_flag == TRUE ~
          "Low confidence (noise-dominated)",
        
        # Intermediate but usable
        Robustness_Score >= 0.6 &
          Robustness_Score < 0.75 ~
          "Moderate confidence",
        
        TRUE ~
          "Low confidence"
      )
    )

  # 10. UPDATED AUDIT TABLES (WITH CONFIDENCE FLAGS)
  
  # 10.1 ALL-SITES AUDIT 
  
  audit_all <- global_summary %>%
    dplyr::select(
      Element,
      Model,
      Global_Rank,
      Robust_Rank,
      Robustness_Score,
      Robustness_Class,
      Confidence_Flag,
      R2,
      RMSE,
      RMSE_ppm,
      Mean_ppm,
      Mean_CI_width,
      SNR_model,
      SNR_CI,
      SNR_smooth,
      Roughness,
      Unstable_flag,
      Is_Production_Model
    ) %>%
    dplyr::arrange(Element, Global_Rank)
  
  write.csv(
    audit_all,
    file.path(
      robust_root,
      "AllElementsModel_Robustness_Stability_ALL.csv"
    ),
    row.names = FALSE
  )
  
  # 10.2 PER-SITE AUDIT TABLES 
  
  # Uses preds_store to infer which models were applied at each site
  
  audit_site_list <- list()
  
  sites_all <- sort(unique(unlist(
    lapply(preds_store, function(x)
      if ("Site" %in% names(x)) unique(x$Site))
  )))
  
  for (s in sites_all) {
    
    audit_site <- audit_all %>%
      dplyr::mutate(Site = s) %>%
      dplyr::select(
        Site,
        dplyr::everything()
      )
    
    write.csv(
      audit_site,
      file.path(
        robust_root,
        paste0(
          "AllElementsModel_Robustness_Stability_Site_",
          s,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    audit_site_list[[s]] <- audit_site
  }
  
  message("✓ Robustness & stability audit tables written (ALL + per site)")
  

  # 11. MASTER ROBUSTNESS AUDIT TABLE
  
  audit_master <- dplyr::bind_rows(
    audit_all %>% dplyr::mutate(Site = "ALL"),
    dplyr::bind_rows(audit_site_list)
  ) %>%
    dplyr::select(
      Site,
      Element,
      Model,
      Global_Rank,
      Robust_Rank,
      Robustness_Score,
      Robustness_Class,
      Confidence_Flag,
      R2,
      RMSE,
      RMSE_ppm,
      Mean_ppm,
      Mean_CI_width,
      SNR_model,
      SNR_CI,
      SNR_smooth,
      Roughness,
      Unstable_flag,
      Is_Production_Model
    ) %>%
    dplyr::arrange(Site, Element, Global_Rank)
  
  write.csv(
    audit_master,
    file.path(
      all_base,
      "AllElementsModel_Robustness_Stability_MASTER.csv"
    ),
    row.names = FALSE
  )
  
  message("✓ Confidence flags added and MASTER robustness audit table written")
  
  # ===============================================================
  # SECTION 11B: ROBUSTNESS VISUAL SUMMARY - OVERALL & PER SITE
  # Heatmap (all elements × models) + Radar plots (per element)
  # ===============================================================
  
  pdf(
    file = file.path(robust_root, "Robustness.pdf"),
    width = 11,
    height = 8.5
  )
  
  # 1. SNR HEATMAP (white → red)
  
  p_heat <- ggplot(
    global_summary,
    aes(x = Model, y = Element, fill = SNR_CI_s)
  ) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low  = "white",
      high = "red",
      limits = c(0, 1),
      name = "Scaled SNR\n(CI-based)",
      na.value = "grey90"
    ) +
    labs(
      title = "Model signal-to-noise robustness (ppm scale)",
      x = "Model",
      y = "Element"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(face = "bold")
    )
  
  print(p_heat)
  
  
  # # 2. STABILITY CLASS HEATMAP (white → red)
  # 
  # p_heat_stability <- ggplot(
  #   global_summary,
  #   aes(
  #     x = Model,
  #     y = Element,
  #     fill = dplyr::case_when(
  #       Robustness_Class == "Unstable"        ~ 0,
  #       grepl("Acceptable", Robustness_Class) ~ 0.5,
  #       grepl("Preferred", Robustness_Class)  ~ 1,
  #       TRUE                                  ~ NA_real_
  #     )
  #   )
  # ) +
  #   geom_tile(colour = "white", linewidth = 0.3) +
  #   scale_fill_gradient(
  #     low  = "white",
  #     high = "red",
  #     limits = c(0, 1),
  #     name = "Stability class\n(Unstable → Preferred)",
  #     na.value = "grey90"
  #   ) +
  #   labs(
  #     title = "Model stability classification",
  #     x = "Model",
  #     y = "Element"
  #   ) +
  #   theme_minimal(base_size = 9) +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     plot.title  = element_text(face = "bold")
  #   )
  # 
  # print(p_heat_stability)
  # 
  # 
  # # 3. CONFIDENCE CLASS HEATMAP (white → red)
  # 
  # p_heat_confidence <- ggplot(
  #   global_summary,
  #   aes(
  #     x = Model,
  #     y = Element,
  #     fill = dplyr::case_when(
  #       grepl("High confidence", Confidence_Flag)     ~ 1,
  #       grepl("Moderate confidence", Confidence_Flag) ~ 0.5,
  #       grepl("Low confidence", Confidence_Flag)      ~ 0,
  #       TRUE                                          ~ NA_real_
  #     )
  #   )
  # ) +
  #   geom_tile(colour = "white", linewidth = 0.3) +
  #   scale_fill_gradient(
  #     low  = "white",
  #     high = "red",
  #     limits = c(0, 1),
  #     name = "Confidence class\n(Low → High)",
  #     na.value = "grey90"
  #   ) +
  #   labs(
  #     title = "Model confidence classification",
  #     x = "Model",
  #     y = "Element"
  #   ) +
  #   theme_minimal(base_size = 9) +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     plot.title  = element_text(face = "bold")
  #   )
  # 
  # print(p_heat_confidence)
  
  
  # 3. RADAR PLOTS — ONE PDF PER ELEMENT (3×3 GRID OF SITES)
  
  sites_order <- c("ALL", "BI10", "HERPB42", "KER1", "KER3", "PB1")
  sites_order <- sites_order[sites_order %in% snr_df$Site]
  
  for (el in elements) {
    
    radar_pdf <- file.path(
      robust_root,
      paste0("Robustness_Radar_", el, ".pdf")
    )
    
    pdf(radar_pdf, width = 11, height = 11)
    
    layout(matrix(1:9, nrow = 3, byrow = TRUE))
    par(mar = c(2, 2, 3, 2))
    
    plot_count <- 0
    
    for (s in sites_order) {
      
      df_r <- snr_df %>%
        dplyr::filter(Element == el, Site == s) %>%
        dplyr::select(
          Model,
          SNR_CI_s,
          SNR_smooth_s,
          Robustness_score
        ) %>%
        tidyr::drop_na()
      
      if (nrow(df_r) < 2) {
        plot.new()
        title(paste(el, "-", s, "\n(no data)"))
        plot_count <- plot_count + 1
        next
      }
      
      radar_vals <- as.data.frame(df_r[, -1])
      rownames(radar_vals) <- df_r$Model
      
      radar_mat <- rbind(
        max = rep(1, ncol(radar_vals)),
        min = rep(0, ncol(radar_vals)),
        radar_vals
      )
      
      radar_cols <- scales::hue_pal()(nrow(radar_vals))
      
      plot.new()
      pdf(file.path(robust_root, paste0("Robustness_Radars_", el, ".pdf")),
          width = 11, height = 11)
      par(bg = "white", fg = "black")
      
      fmsb::radarchart(
        radar_mat,
        axistype = 1,
        pcol  = cols,
        pfcol = scales::alpha(cols, 0.35),
        plwd  = 2,
        cglcol = "grey80",
        cglty  = 1,
        axislabcol = "grey40",
        vlcex = 0.7,
        title = paste(el, "-", s)
      )
      
      plot_count <- plot_count + 1
    }
    
    # pad remaining panels to 3×3
    while (plot_count < 9) {
      plot.new()
      plot_count <- plot_count + 1
    }
    
    dev.off()
    
    message("✓ Radar written: ", radar_pdf)
  }
  
  
  # 4. One radar per MODEL  across ELEMENTS), 3×3 sites grid
  
  message("=== Creating model-centric robustness radar summaries ===")
  
  model_radar_dir <- file.path(robust_root, "Model_Radars")
  dir.create(model_radar_dir, showWarnings = FALSE, recursive = TRUE)
  
  sites_order <- c("ALL", "BI10", "HERPB42", "KER1", "KER3", "PB1")
  sites_order <- sites_order[sites_order %in% snr_df$Site]
  
  for (mname in unique(snr_df$Model)) {
    
    radar_grobs <- list()
    
    for (s in sites_order) {
      
      df_r <- snr_df %>%
        dplyr::filter(Model == mname, Site == s) %>%
        dplyr::select(
          Element,
          SNR_CI_s,
          SNR_smooth_s,
          Robustness_score
        ) %>%
        tidyr::drop_na()
      
      if (nrow(df_r) < 2) {
        radar_grobs[[s]] <- grid::nullGrob()
        next
      }
      
      radar_vals <- as.data.frame(df_r[, -1])
      rownames(radar_vals) <- df_r$Element
      
      radar_mat <- rbind(
        max = rep(1, ncol(radar_vals)),
        min = rep(0, ncol(radar_vals)),
        radar_vals
      )
      
      radar_cols <- scales::hue_pal()(nrow(radar_vals))
      
      radar_grobs[[s]] <- grid::grid.grabExpr({
        fmsb::radarchart(
          radar_mat,
          axistype   = 1,
          pcol       = radar_cols,
          pfcol      = scales::alpha(radar_cols, 0.35),
          plwd       = 2,
          cglcol     = "grey80",
          cglty      = 1,
          axislabcol = "grey40",
          vlcex      = 0.7,
          title      = paste("Model:", mname, "| Site:", s)
        )
      })
    }
    
    # Pad to 3×3
    while (length(radar_grobs) < 9) {
      radar_grobs[[paste0("blank_", length(radar_grobs))]] <-
        grid::nullGrob()
    }
    
    radar_file <- file.path(
      model_radar_dir,
      paste0("Robustness_Radar_Model_", mname, ".pdf")
    )
    
    pdf(radar_file, width = 11, height = 11)
    grid::grid.newpage()
    gridExtra::grid.arrange(
      grobs = radar_grobs,
      ncol  = 3,
      nrow  = 3
    )
    dev.off()
    
    message("✓ Model radar written: ", basename(radar_file))
  }
  
  message("=== Model-centric radar summaries complete ===")
  
  message("✓ Robustness heatmap (white→red) + radar plots written to Robustness.pdf")
  
  # ===============================================================
  # SECTION 11C:ROBUSTNESS & SNR ANALYSIS (OVERALL + PER SITE)
  # OUTPUT: Robustness folder
  # ===============================================================

  # Output directory
  
  robust_root <- file.path(all_base, "Robustness")
  
  if (!dir.exists(robust_root)) {
    dir.create(robust_root, recursive = TRUE, showWarnings = FALSE)
  }
  stopifnot(dir.exists(robust_root))
  
  message("Robustness outputs will be written to: ",
          normalizePath(robust_root, winslash = "/"))
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(fmsb)
    library(grid)
    library(patchwork)
  })
  
  # Helper functions

  scale01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (!is.finite(diff(r)) || diff(r) == 0)
      return(rep(NA_real_, length(x)))
    (x - r[1]) / diff(r)
  }
  
  calc_snr_ci <- function(pred, lwr, upr) {
    signal <- sd(pred, na.rm = TRUE)
    noise  <- mean(upr - lwr, na.rm = TRUE)
    ifelse(is.na(signal) | is.na(noise) | noise == 0,
           NA_real_, signal / noise)
  }
  
  calc_snr_smooth <- function(pred, depth = NULL) {
    if (length(pred) < 6) return(NA_real_)
    if (!is.null(depth)) pred <- pred[order(depth)]
    d <- diff(pred)
    if (all(is.na(d))) return(NA_real_)
    1 / sd(d, na.rm = TRUE)
  }
  
# for building heatmaps
  class_to_numeric <- function(x) {
    if (is.numeric(x)) return(x)
    f <- factor(x)
    as.numeric(f) / max(as.numeric(f), na.rm = TRUE)
  }
  
  # Build robustness table

  rows <- list()
  
  for (el in elements) {
    
    df0 <- preds_store[[el]]
    if (is.null(df0) || !"Site" %in% names(df0)) next
    
    sites <- c("ALL", sort(unique(df0$Site)))
    
    for (s in sites) {
      
      df_s <- if (s == "ALL") df0 else df0 %>% filter(Site == s)
      if (nrow(df_s) < 6) next
      
      for (mname in names(model_store[[el]])) {
        
        pred_col <- paste0(mname, "_", el, "_Pred_ppm")
        lwr_col  <- paste0(mname, "_", el, "_L95_ppm")
        upr_col  <- paste0(mname, "_", el, "_U95_ppm")
        
        if (!all(c(pred_col, lwr_col, upr_col) %in% names(df_s))) next
        
        ok <- is.finite(df_s[[pred_col]]) &
          is.finite(df_s[[lwr_col]])  &
          is.finite(df_s[[upr_col]])
        
        if (sum(ok) < 6) next
        
        pred <- df_s[[pred_col]][ok]
        lwr  <- df_s[[lwr_col]][ok]
        upr  <- df_s[[upr_col]][ok]
        
        depth_vec <- if ("depth" %in% names(df_s)) df_s$depth[ok] else NULL
        
        rows[[length(rows) + 1]] <- data.frame(
          Element    = el,
          Site       = s,
          Model      = mname,
          n          = length(pred),
          Mean_ppm   = mean(pred, na.rm = TRUE),
          SD_pred    = sd(pred, na.rm = TRUE),
          Mean_CI_w  = mean(upr - lwr, na.rm = TRUE),
          SNR_CI     = calc_snr_ci(pred, lwr, upr),
          SNR_smooth = calc_snr_smooth(pred, depth_vec),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  snr_df <- bind_rows(rows)
  
  if (nrow(snr_df) == 0) {
    stop("No robustness results produced — check preds_store content.")
  }
  
  # Scale + composite robustness

  snr_df <- snr_df %>%
    group_by(Element, Site) %>%
    mutate(
      SNR_CI_s         = scale01(SNR_CI),
      SNR_smooth_s     = scale01(SNR_smooth),
      Robustness_score = 0.6 * SNR_CI_s + 0.4 * SNR_smooth_s
    ) %>%
    arrange(desc(Robustness_score), .by_group = TRUE) %>%
    mutate(Robustness_Rank = row_number()) %>%
    ungroup()
  
  # Write CSV outputs

  write.csv(
    snr_df,
    file.path(robust_root, "Robustness_ALL_SITES.csv"),
    row.names = FALSE
  )
  
  for (s in unique(snr_df$Site)) {
    write.csv(
      snr_df %>% filter(Site == s),
      file.path(robust_root, paste0("Robustness_Site_", s, ".csv")),
      row.names = FALSE
    )
  }
  
  message("✓ Robustness CSVs written")
  
  
# HEATMAPS — ROBUSTNESS - one combined 3×3 page

  sites_order <- c("ALL", "BI10", "HER42PB", "KER1", "KER3", "PB1")
  sites_order <- sites_order[sites_order %in% snr_df$Site]
  
  heatmap_plots <- list()
  
  for (s in sites_order) {
    
    df_h <- snr_df %>% filter(Site == s)
    if (nrow(df_h) == 0) next
    
    heatmap_plots[[s]] <- ggplot(
      df_h,
      aes(x = Model, y = Element, fill = Robustness_score)
    ) +
      geom_tile(color = "grey85", linewidth = 0.3) +
      scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
      labs(title = paste("Robustness: ", s)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  while (length(heatmap_plots) < 9) {
    heatmap_plots[[paste0("blank_", length(heatmap_plots))]] <-
      ggplot() + theme_void()
  }
  
  pdf(file.path(robust_root, "Heatmaps_Robustness_AllSites.pdf"), 11, 11)
  
  print(
    wrap_plots(
      heatmap_plots,
      ncol   = 3,
      nrow   = 3,
      guides = "collect"   # <<< COLLECT LEGENDS
    ) &
      theme(
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.box       = "horizontal"
      )
  )
  
  dev.off()
  
  
# STABILITY CLASS HEATMAPS — ALL + PER SITE
  
  # DEFINE STABILITY CLASS (EXPLICIT) - for heatmapping
  
  snr_df <- snr_df %>%
    dplyr::mutate(
      Stability_class = dplyr::case_when(
        Robustness_score >= 0.75                    ~ "Stable",
        Robustness_score >= 0.50 & Robustness_score < 0.75 ~ "Acceptable",
        Robustness_score >= 0.30 & Robustness_score < 0.50 ~ "Marginal",
        TRUE                                        ~ "Unstable"
      )
    )
  
  # FIXED CLASS → NUMERIC MAPPINGS FOR HEATMAPS 
  stability_to_score <- function(x) {
    dplyr::case_when(
      x == "Stable"      ~ 1.00,
      x == "Acceptable" ~ 0.67,
      x == "Marginal"   ~ 0.33,
      x == "Unstable"   ~ 0.00,
      TRUE              ~ NA_real_
    )
  }
  
  
  heatmap_plots <- list()
  
  for (s in sites_order) {
    
    df_h <- snr_df %>%
      dplyr::filter(Site == s)
    
    if (nrow(df_h) == 0) next
    
    df_h <- df_h %>%
      dplyr::mutate(
        Stability_score_hm = stability_to_score(Stability_class)
      )
    
    heatmap_plots[[s]] <- ggplot(
      df_h,
      aes(x = Model, y = Element, fill = Stability_score_hm)
    ) +
      geom_tile(color = "grey85", linewidth = 0.3) +
      scale_fill_gradient(
        low    = "white",
        high   = "red",
        limits = c(0, 1),
        breaks = c(0, 0.33, 0.67, 1),
        labels = c("Unstable", "Marginal", "Acceptable", "Stable"),
        name   = "Stability class",
        na.value = "grey90"
      ) +
      labs(
        title = paste("Model stability:", s),
        x = NULL,
        y = NULL
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
  
  while (length(heatmap_plots) < 9) {
    heatmap_plots[[paste0("blank_", length(heatmap_plots))]] <-
      ggplot() + theme_void()
  }
  
  pdf(file.path(robust_root, "Heatmaps_Stability_AllSites.pdf"), 11, 11)
  
  print(
    patchwork::wrap_plots(
      heatmap_plots,
      ncol   = 3,
      nrow   = 3,
      guides = "collect"
    ) &
      theme(
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.box       = "horizontal"
      )
  )
  
  dev.off()
  
# CONFIDENCE CLASS HEATMAPS — ALL + PER SITE
  
  # DEFINE CONFIDENCE CLASS (EXPLICIT) FOR HEATMAPS
  
  snr_df <- snr_df %>%
    dplyr::mutate(
      Confidence_class = dplyr::case_when(
        SNR_CI_s >= 0.67 & SNR_smooth_s >= 0.67 ~ "Signal-dominated",
        SNR_CI_s <  0.33 | SNR_smooth_s <  0.33 ~ "Uncertainty-dominated",
        TRUE                                   ~ "Balanced"
      )
    ) 
  
  confidence_to_score <- function(x) {
    dplyr::case_when(
      x == "Signal-dominated"       ~ 1.00,
      x == "Balanced"               ~ 0.50,
      x == "Uncertainty-dominated"  ~ 0.00,
      TRUE                          ~ NA_real_
    )
  }
  
  heatmap_plots <- list()
  
  for (s in sites_order) {
    
    df_h <- snr_df %>%
      dplyr::filter(Site == s)
    
    if (nrow(df_h) == 0) next
    
    df_h <- df_h %>%
      dplyr::mutate(
        Confidence_score_hm = confidence_to_score(Confidence_class)
      )
    
    heatmap_plots[[s]] <- ggplot(
      df_h,
      aes(x = Model, y = Element, fill = Confidence_score_hm)
    ) +
      geom_tile(color = "grey85", linewidth = 0.3) +
      scale_fill_gradient(
        low    = "white",
        high   = "red",
        limits = c(0, 1),
        breaks = c(0, 0.5, 1),
        labels = c(
          "Uncertainty-dominated",
          "Balanced",
          "Signal-dominated"
        ),
        name   = "Confidence class",
        na.value = "grey90"
      ) +
      labs(
        title = paste("Prediction confidence:", s),
        x = NULL,
        y = NULL
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
  
  while (length(heatmap_plots) < 9) {
    heatmap_plots[[paste0("blank_", length(heatmap_plots))]] <-
      ggplot() + theme_void()
  }
  
  pdf(file.path(robust_root, "Heatmaps_Confidence_AllSites.pdf"), 11, 11)
  
  print(
    patchwork::wrap_plots(
      heatmap_plots,
      ncol   = 3,
      nrow   = 3,
      guides = "collect"
    ) &
      theme(
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.box       = "horizontal"
      )
  )
  
  dev.off()
  
  # RADAR PLOTS — ONE PDF PER ELEMENT (3×3 grid of sites)
  
  # cb_palette <- c(
  #   "#E69F00", "#56B4E9", "#009E73",
  #   "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  # )
  # 
  # # for (el in elements) {
  # #   
  # #   pdf(
  # #     file.path(robust_root, paste0("Robustness_Radars_El_", el, ".pdf")),
  # #     width = 11, height = 11
  # #   )
  # #   
  # #   layout(matrix(1:9, nrow = 3, byrow = TRUE))
  # #   par(mar = c(2, 2, 3, 2))
  # #   
  # #   plot_count <- 0
  #   
  #   for (s in sites_order) {
  #     
  #     df_r <- snr_df %>% filter(Element == el, Site == s)
  #     
  #     if (nrow(df_r) < 2) {
  #       plot.new()
  #       title(paste(el, "-", s, "\n(no data)"))
  #       plot_count <- plot_count + 1
  #       next
  #     }
  #     
  #     radar_vals <- df_r %>%
  #       dplyr::select(SNR_CI_s, SNR_smooth_s, Robustness_score) %>%
  #       as.data.frame()
  #     
  #     rownames(radar_vals) <- df_r$Model
  #     
  #     radar_mat <- rbind(
  #       max = rep(1, ncol(radar_vals)),
  #       min = rep(0, ncol(radar_vals)),
  #       radar_vals
  #     )
  #     
  #     cols <- cb_palette[seq_len(nrow(radar_vals))]
  #     
  #     plot.new()
  #     par(bg = "white", fg = "black")
  #     
  #     fmsb::radarchart(
  #       radar_mat,
  #       axistype = 1,
  #       pcol = cols,
  #       pfcol = scales::alpha(cols, 0.35),
  #       plwd = 2,
  #       cglcol = "grey80",
  #       title = paste(el, "-", s)
  #     )
  #     
  #     plot_count <- plot_count + 1
  #     if (plot_count >= 9) break
  #   }
  #   
  #   while (plot_count < 9) {
  #     plot.new()
  #     plot_count <- plot_count + 1
  #   }
  #   
  #   dev.off()
  # }
  # 
  message("DONE — all robustness outputs saved to Robustness/")
  
  # ===============================================================
  # SECTION 11D: SITE SUMMARY CSVs — RANK MODELS PER SITE × ELEMENT
  # Priority:
  #   1) Robustness_score (site)
  #   2) R2_overall
  #   3) RMSE_overall
  #   4) RMSEP_overall
  # ===============================================================
  
  if (!exists("snr_df") || nrow(snr_df) == 0) {
    warning("snr_df missing/empty — skipping per-site robustness ranking summaries.")
  } else if (!exists("global_summary") || nrow(global_summary) == 0) {
    warning("global_summary missing/empty — skipping per-site robustness ranking summaries.")
  } else {
    
    # Ensure Rank exists in global_summary
    if (!"Rank" %in% names(global_summary)) {
      global_summary <- global_summary %>%
        dplyr::group_by(Element) %>%
        dplyr::arrange(dplyr::desc(R2), RMSEP, RMSE, RMSE_ppm) %>%
        dplyr::mutate(Rank = dplyr::row_number()) %>%
        dplyr::ungroup()
    }
    
    # Pull overall metrics used for ranking
    global_metrics <- global_summary %>%
      dplyr::select(
        Element,
        Model,
        R2_overall    = R2,
        RMSE_overall  = RMSE,
        RMSEP_overall = RMSE_ppm
      )
    
    # Summarise robustness safely
    robust_site_model <- snr_df %>%
      dplyr::group_by(Site, Element, Model) %>%
      dplyr::summarise(
        Robustness_score = mean(Robustness_score, na.rm = TRUE),
        SNR_CI_s         = mean(SNR_CI_s, na.rm = TRUE),
        SNR_smooth_s     = mean(SNR_smooth_s, na.rm = TRUE),
        n_mean           = round(mean(n, na.rm = TRUE), 0),
        .groups = "drop"
      ) %>%
      dplyr::left_join(global_metrics, by = c("Element", "Model"))
    
    # Rank models per Site × Element
    site_ranked <- robust_site_model %>%
      dplyr::group_by(Site, Element) %>%
      dplyr::arrange(
        dplyr::desc(Robustness_score),
        dplyr::desc(R2_overall),
        RMSE_overall,
        RMSEP_overall
      ) %>%
      dplyr::mutate(Rank_site = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    # Write one CSV per site (including ALL)
    sites_out <- sort(unique(site_ranked$Site))
    
    for (s in sites_out) {
      
      out_df <- site_ranked %>%
        dplyr::filter(Site == s) %>%
        dplyr::select(
          Site,
          Element,
          Model,
          Rank_site,
          Robustness_score,
          R2_overall,
          RMSE_overall,
          RMSEP_overall,
          SNR_CI_s,
          SNR_smooth_s,
          n_mean
        ) %>%
        dplyr::arrange(Element, Rank_site)
      
      write.csv(
        out_df,
        file.path(robust_root, paste0("SiteModelRank_", s, ".csv")),
        row.names = FALSE
      )
    }
    
    # Combined all-sites table
    write.csv(
      site_ranked %>%
        dplyr::select(
          Site,
          Element,
          Model,
          Rank_site,
          Robustness_score,
          R2_overall,
          RMSE_overall,
          RMSEP_overall,
          SNR_CI_s,
          SNR_smooth_s,
          n_mean
        ) %>%
        dplyr::arrange(Site, Element, Rank_site),
      file.path(robust_root, "SiteModelRank_ALLSITES.csv"),
      row.names = FALSE
    )
    
    message("✓ Per-site model ranking CSVs written with RMSEP included")
  }
  
  # ===============================================================
  # 12. FINISH
  # ===============================================================
  
  message("=== COMPLETED SUCCESSFULLY ===")
  message("All results saved in: ", out_dir)
  
  # ===============================================================
  # Sections 13–20: USER OPTION FOR ACE PREDICTIONS (TRUE/FALSE) 
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
  # 13. ACE_dataset PREDICTION - XRF_pred as input 
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
      
      # Access fitted models + RMSE from calibration step (Section 8)
      models_el <- model_store[[el]]
      rmse_el   <- calibration_rmse[[el]]
      
      # Compute predictions, store per model
      for (mname in names(models_el)) {
        
        m <- models_el[[mname]]
        
        # Predictions in log space
        pred_log <- tryCatch({
          if (inherits(m, "mvr")) {
            as.numeric(predict(m,
                               ncomp   = attr(m, "best_ncomp"),
                               newdata = ACE_dataset))
          } else {
            as.numeric(predict(m, newdata = ACE_dataset))
          }
        }, error = function(e) rep(NA_real_, nrow(ACE_dataset)))
        
        # Compute CI from calibration RMSE
        rmse_val <- rmse_el[[mname]]
        if (is.null(rmse_val) || is.na(rmse_val)) rmse_val <- NA_real_
        
        L95_log <- pred_log - 1.96 * rmse_val
        U95_log <- pred_log + 1.96 * rmse_val
        
        # ppm-scale values (keep original names for backward compatibility)
        base_nm <- paste0(mname, "_", el)
        
        ace_df[[paste0(base_nm, "_Pred_ppm")]] <- exp(pred_log)
        ace_df[[paste0(base_nm, "_L95_ppm")]]  <- exp(L95_log)
        ace_df[[paste0(base_nm, "_U95_ppm")]]  <- exp(U95_log)
        
        # Duplicate into *_XRF_pred* naming for error-bar use
        # (these are the same ppm values, just different column names)
        ace_df[[paste0(base_nm, "_XRF_pred")]]        <- ace_df[[paste0(base_nm, "_Pred_ppm")]]
        ace_df[[paste0(base_nm, "_XRF_pred_lower")]]  <- ace_df[[paste0(base_nm, "_L95_ppm")]]
        ace_df[[paste0(base_nm, "_XRF_pred_upper")]]  <- ace_df[[paste0(base_nm, "_U95_ppm")]]
        
        # ============================================================
        # LOG-SPACE % ERROR + PPM-SCALE % ERROR INTERVALS (ACE)
        # ============================================================
        
        safe_pred_log <- ifelse(is.na(pred_log) | pred_log == 0, NA, pred_log)
        
        # log-space percentage errors
        L95_log_pc <- (safe_pred_log - L95_log) / safe_pred_log * 100
        U95_log_pc <- (U95_log - safe_pred_log) / safe_pred_log * 100
        
        # store log-space percentages
        ace_df[[paste0(mname, "_", el, "_L95_log_pc")]] <- L95_log_pc
        ace_df[[paste0(mname, "_", el, "_U95_log_pc")]] <- U95_log_pc
        
        # ppm-scale percentage equivalents
        pred_ppm <- exp(pred_log)
        
        L95_pc_ppm <- pred_ppm - (pred_ppm * (L95_log_pc / 100))
        U95_pc_ppm <- pred_ppm + (pred_ppm * (U95_log_pc / 100))
        
        ace_df[[paste0(mname, "_", el, "_L95_pc_ppm")]] <- L95_pc_ppm
        ace_df[[paste0(mname, "_", el, "_U95_pc_ppm")]] <- U95_pc_ppm
        
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
  # 14. ACE DIAGNOSTIC EVALUATION vs ICP_ppm  (FIXED FOR PPM METRICS)
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
      
      # Merge ACE predictions with observed ICP ppm
      joined <- ace_df %>%
        dplyr::left_join(
          ICP_ppm %>% dplyr::select(Site, depth, SH20_age, !!sym(obs_col), !!sym(sd_col)),
          by = c("Site", "depth", "SH20_age")
        )
      
      diag_out <- data.frame()
      
      for (mname in names(model_store[[el]])) {
        
        pred_nm <- paste0(mname, "_", el, "_Pred_ppm")
        if (!pred_nm %in% names(joined)) next
        
        pred <- joined[[pred_nm]]
        obs  <- joined[[obs_col]]
        
        ok <- !is.na(pred) & !is.na(obs)
        if (!any(ok)) next
        
        # -------------------------------
        # PPM-BASED METRICS (CORRECT)
        # -------------------------------
        SS_res <- sum((obs[ok] - pred[ok])^2, na.rm = TRUE)
        SS_tot <- sum((obs[ok] - mean(obs[ok], na.rm = TRUE))^2, na.rm = TRUE)
        
        R2_ppm   <- ifelse(SS_tot == 0, NA_real_, 1 - SS_res / SS_tot)
        RMSE_ppm <- sqrt(mean((obs[ok] - pred[ok])^2, na.rm = TRUE))
        
        diag_out <- rbind(diag_out,
                          data.frame(
                            Model = mname,
                            R2    = R2_ppm,
                            RMSE  = RMSE_ppm
                          ))
        
      }
      
      write.csv(
        diag_out,
        file.path(el_dir, paste0(el, "_ACE_diagnostics.csv")),
        row.names = FALSE
      )
    }
  }
  # ===============================================================
  # 15. ACE PREDICTION DIAGNOSTIC PLOTS (FINAL SCIENTIFIC + EQUAL AXES)
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== Creating ACE diagnostic plots (FINAL Section 15) ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    
    if (!exists("p_po_list_ace"))
      p_po_list_ace <- list()
    
    for (el in elements) {
      
      el_dir  <- file.path(ace_root, el)
      pred_fl <- file.path(el_dir, paste0(el, "_ACE_predictions.csv"))
      if (!file.exists(pred_fl)) next
      
      ace_df <- read.csv(pred_fl)
      
      # Observed ICP column + SD column (used for horizontal error bars)
      obs_col <- paste0(el, "_ICP")
      sd_col  <- paste0(el, "_ICP_sd")
      
      if (!(obs_col %in% names(ICP_ppm))) next
      
      joined <- ace_df %>%
        dplyr::left_join(
          ICP_ppm %>% dplyr::select(
            Site, depth, SH20_age,
            !!sym(obs_col),
            !!sym(sd_col)
          ),
          by = c("Site", "depth", "SH20_age")
        )
      
      for (mname in names(model_store[[el]])) {
        
        pred_nm <- paste0(mname, "_", el, "_Pred_ppm")
        Lnm     <- paste0(mname, "_", el, "_L95_ppm")
        Unm     <- paste0(mname, "_", el, "_U95_ppm")
        Lpc_nm  <- paste0(mname, "_", el, "_L95_pc_ppm")
        Upc_nm  <- paste0(mname, "_", el, "_U95_pc_ppm")
        
        if (!pred_nm %in% names(joined)) next
        
        df <- joined %>%
          dplyr::filter(!is.na(.data[[obs_col]]),
                        !is.na(.data[[pred_nm]]))
        if (nrow(df) == 0) next
        
        obs  <- df[[obs_col]]
        pred <- df[[pred_nm]]
        
        SS_res <- sum((obs - pred)^2, na.rm = TRUE)
        SS_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
        
        R2_ppm   <- ifelse(SS_tot == 0, NA_real_, 1 - SS_res / SS_tot)
        RMSE_ppm <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
        
        stats_label <- paste0(
          "R² pred = ", ifelse(is.na(R2_ppm), "NA", signif(R2_ppm, 3)),
          "\nRMSE pred = ", signif(RMSE_ppm, 4), " ppm"
        )
        
        ann_x <- -Inf
        ann_y <-  Inf
        
        if (!el %in% names(p_po_list_ace)) p_po_list_ace[[el]] <- list()
        if (!mname %in% names(p_po_list_ace[[el]]))
          p_po_list_ace[[el]][[mname]] <- list()
        
        # ===========================================================
        # 1. NO ERROR
        # ===========================================================
        
        x_min <- min(df[[obs_col]],  na.rm = TRUE)
        x_max <- max(df[[obs_col]],  na.rm = TRUE)
        y_min <- min(df[[pred_nm]], na.rm = TRUE)
        y_max <- max(df[[pred_nm]], na.rm = TRUE)
        
        axis_min <- min(x_min, y_min)
        axis_max <- max(x_max, y_max)
        
        pad <- 0.03 * (axis_max - axis_min)
        axis_min_pad <- axis_min - pad
        axis_max_pad <- axis_max + pad
        
        # --- annotation coordinates: top-left INSIDE the padded limits ---
        ann_x <- axis_min_pad + 0.02 * (axis_max_pad - axis_min_pad)
        ann_y <- axis_max_pad - 0.02 * (axis_max_pad - axis_min_pad)
        
        p_no_err <- ggplot(df, aes(x = .data[[obs_col]],
                                   y = .data[[pred_nm]],
                                   colour = Site)) +
          geom_point(size = 1.4) +
          scale_x_continuous(labels = scales::label_scientific(),
                             limits = c(axis_min_pad, axis_max_pad)) +
          scale_y_continuous(labels = scales::label_scientific(),
                             limits = c(axis_min_pad, axis_max_pad)) +
          ggsci::scale_color_jco() +
          geom_abline(slope = 1, intercept = 0) +
          annotate(
            "text",
            x = ann_x, y = ann_y,
            hjust = 0, vjust = 1,
            label = stats_label,
            size = 2.7
          )  +
          labs(
            x = expression("Observed ICP-MS (mg kg"^{-1}*")"),
            y = expression("Predicted XRF-calibrated (mg kg"^{-1}*")"),
            title = paste(el, mname, "- External Validation")
          ) +
          theme_bw() + theme_small +
          coord_fixed()
        
        ggsave(
          file.path(el_dir, paste0(el, "_", mname, "_ACE_PredVsObs_no_err.pdf")),
          p_no_err, width = 12, height = 12, units = "cm"
        )
        
        p_po_list_ace[[el]][[mname]]$no_err <- p_no_err
        
        # ===========================================================
        # 2. ABSOLUTE ERROR (±2×ICP_sd + model CI)
        # ===========================================================
        
        # <<< PATCH >>> convert CI bounds to centred asymmetric errors
        df$Pred_err_low  <- df[[pred_nm]] - df[[Lnm]]
        df$Pred_err_high <- df[[Unm]]     - df[[pred_nm]]
        
        x_min <- min(df[[obs_col]] - 2 * df[[sd_col]], na.rm = TRUE)
        x_max <- max(df[[obs_col]] + 2 * df[[sd_col]], na.rm = TRUE)
        y_min <- min(df[[pred_nm]] - df$Pred_err_low,  na.rm = TRUE)
        y_max <- max(df[[pred_nm]] + df$Pred_err_high, na.rm = TRUE)
        
        axis_min <- min(x_min, y_min)
        axis_max <- max(x_max, y_max)
        
        pad <- 0.03 * (axis_max - axis_min)
        axis_min_pad <- axis_min - pad
        axis_max_pad <- axis_max + pad
        
        # --- annotation coordinates: top-left INSIDE the padded limits ---
        ann_x <- axis_min_pad + 0.02 * (axis_max_pad - axis_min_pad)
        ann_y <- axis_max_pad - 0.02 * (axis_max_pad - axis_min_pad)
        
        if (all(c(sd_col, Lnm, Unm) %in% names(df))) {
          
          p_abs_err <- ggplot(df, aes(x = .data[[obs_col]],
                                      y = .data[[pred_nm]],
                                      colour = Site)) +
            geom_segment(aes(
              x = .data[[obs_col]] - 2 * .data[[sd_col]],
              xend = .data[[obs_col]] + 2 * .data[[sd_col]],
              y = .data[[pred_nm]],
              yend = .data[[pred_nm]]
            ), colour = "grey80") +
            geom_segment(aes(
              x = .data[[obs_col]],
              xend = .data[[obs_col]],
              y = .data[[pred_nm]] - Pred_err_low,
              yend = .data[[pred_nm]] + Pred_err_high
            ), colour = "grey80") +
            geom_point(size = 1.4) +
            scale_x_continuous(labels = scales::label_scientific(),
                               limits = c(axis_min_pad, axis_max_pad)) +
            scale_y_continuous(labels = scales::label_scientific(),
                               limits = c(axis_min_pad, axis_max_pad)) +
            ggsci::scale_color_jco() +
            geom_abline(slope = 1, intercept = 0) +
            annotate(
              "text",
              x = ann_x, y = ann_y,
              hjust = 0, vjust = 1,
              label = stats_label,
              size = 2.7
            )  +
            labs(
              x = expression("Observed ICP-MS (mg kg"^{-1}*")"),
              y = expression("Predicted XRF-calibrated (mg kg"^{-1}*")"),
              title = paste(el, mname, "- External Validation")
            ) +
            theme_bw() + theme_small +
            coord_fixed()
          
          ggsave(
            file.path(el_dir, paste0(el, "_", mname, "_ACE_PredVsObs_abs_err.pdf")),
            p_abs_err, width = 12, height = 12, units = "cm"
          )
          
          p_po_list_ace[[el]][[mname]]$abs_err <- p_abs_err
        }
        
        # ===========================================================
        # 3. PERCENT ERROR (±2×ICP_sd + % error ppm)
        # ===========================================================
        
        x_min <- min(df[[obs_col]] - 2 * df[[sd_col]], na.rm = TRUE)
        x_max <- max(df[[obs_col]] + 2 * df[[sd_col]], na.rm = TRUE)
        y_min <- min(df[[pred_nm]] - df[[Lpc_nm]], na.rm = TRUE)
        y_max <- max(df[[pred_nm]] + df[[Upc_nm]], na.rm = TRUE)
        
        axis_min <- min(x_min, y_min)
        axis_max <- max(x_max, y_max)
        
        pad <- 0.03 * (axis_max - axis_min)
        axis_min_pad <- axis_min - pad
        axis_max_pad <- axis_max + pad
        
        # --- annotation coordinates: top-left INSIDE the padded limits ---
        ann_x <- axis_min_pad + 0.02 * (axis_max_pad - axis_min_pad)
        ann_y <- axis_max_pad - 0.02 * (axis_max_pad - axis_min_pad)
        
        if (all(c(Lpc_nm, Upc_nm) %in% names(df))) {
          
          p_pc_err <- ggplot(df, aes(x = .data[[obs_col]],
                                     y = .data[[pred_nm]],
                                     colour = Site)) +
            geom_segment(aes(
              x = .data[[obs_col]] - 2 * .data[[sd_col]],
              xend = .data[[obs_col]] + 2 * .data[[sd_col]],
              y = .data[[pred_nm]],
              yend = .data[[pred_nm]]
            ), colour = "grey80") +
            geom_segment(aes(
              x = .data[[obs_col]],
              xend = .data[[obs_col]],
              y = .data[[pred_nm]] - .data[[Lpc_nm]],
              yend = .data[[pred_nm]] + .data[[Upc_nm]]
            ), colour = "grey80") +
            geom_point(size = 1.4) +
            scale_x_continuous(labels = scales::label_scientific(),
                               limits = c(axis_min_pad, axis_max_pad)) +
            scale_y_continuous(labels = scales::label_scientific(),
                               limits = c(axis_min_pad, axis_max_pad)) +
            ggsci::scale_color_jco() +
            geom_abline(slope = 1, intercept = 0) +
            annotate(
              "text",
              x = ann_x, y = ann_y,
              hjust = 0, vjust = 1,
              label = stats_label,
              size = 2.7
            ) +
            labs(
              x = expression("Observed ICP-MS (mg kg"^{-1}*")"),
              y = expression("Predicted XRF-calibrated (mg kg"^{-1}*")"),
              title = paste(el, mname, "- External Validation")
            ) +
            theme_bw() + theme_small +
            coord_fixed()
          
          ggsave(
            file.path(el_dir, paste0(el, "_", mname, "_ACE_PredVsObs_pc_err.pdf")),
            p_pc_err, width = 12, height = 12, units = "cm"
          )
          
          p_po_list_ace[[el]][[mname]]$pc_err <- p_pc_err
        }
      }
    }
  }
  # ===============================================================
  # 16. ACE Summary Pages (Square Panels, Scientific Axes, FINAL)
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== Creating ACE summary pages (FINAL Section 16) ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    
    if (!exists("p_po_list_ace")) {
      warning("p_po_list_ace missing — cannot create summary pages.")
    } else {
      
      for (el in elements) {
        
        el_dir <- file.path(ace_root, el)
        if (!dir.exists(el_dir)) next
        if (!el %in% names(p_po_list_ace)) next
        
        model_plots <- p_po_list_ace[[el]]
        if (length(model_plots) == 0) next
        
        # -----------------------------------------------------------
        # Helper: enforce square panels without touching geoms
        # -----------------------------------------------------------
        sq <- function(p) {
          if (is.null(p)) {
            ggplot() + theme_void() + coord_fixed()
          } else {
            p + coord_fixed()
          }
        }
        
        models <- names(model_plots)
        
        no_err_list  <- lapply(models, function(m) sq(model_plots[[m]]$no_err))
        abs_err_list <- lapply(models, function(m) sq(model_plots[[m]]$abs_err))
        pc_err_list  <- lapply(models, function(m) sq(model_plots[[m]]$pc_err))
        
        grid <- function(lst) {
          wrap_plots(lst, ncol = 2, nrow = 4, guides = "collect") &
            theme(
              legend.position  = "bottom",
              legend.direction = "horizontal",
              legend.box       = "horizontal"
            ) &
            guides(colour = guide_legend(nrow = 1))
        }
        
        # -----------------------------------------------------------
        # Save summary pages (A4 portrait, square panels)
        # -----------------------------------------------------------
        pdf(file.path(el_dir, paste0(el, "_ACE_PredObs_Model_Summary_no_err.pdf")),
            width = 21/2.54, height = 29.7/2.54)
        print(grid(no_err_list))
        dev.off()
        
        pdf(file.path(el_dir, paste0(el, "_ACE_PredObs_Model_Summary_abs_err.pdf")),
            width = 21/2.54, height = 29.7/2.54)
        print(grid(abs_err_list))
        dev.off()
        
        pdf(file.path(el_dir, paste0(el, "_ACE_PredObs_Model_Summary_pc_err.pdf")),
            width = 21/2.54, height = 29.7/2.54)
        print(grid(pc_err_list))
        dev.off()
      }
    }
  }
  
  
  # ===============================================================
  # 17. ACE GLOBAL SUMMARY (patched for correct RMSE%, ppm logic)
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== ACE Global Summary ===")
    
    ace_root     <- file.path(all_base, "Calibration_pred")
    summary_file <- file.path(ace_root, "AllElements_ModelSummary_pred.csv")
    best_file    <- file.path(ace_root, "BestModels_PerElement_pred.csv")
    
    ace_global <- data.frame()
    
    for (el in elements) {
      
      el_dir    <- file.path(ace_root, el)
      diag_file <- file.path(el_dir, paste0(el, "_ACE_diagnostics.csv"))
      pred_file <- file.path(el_dir, paste0(el, "_ACE_predictions.csv"))
      
      if (!file.exists(diag_file) || !file.exists(pred_file)) next
      
      d      <- read.csv(diag_file)
      ace_df <- read.csv(pred_file)
      
      if (nrow(d) == 0) next
      
      obs_col <- paste0(el, "_ICP")
      if (!obs_col %in% names(ICP_ppm)) next
      
      ACE_obs <- ICP_ppm[[obs_col]]
      mean_obs_log <- mean(log(ACE_obs), na.rm = TRUE)
      
      # -----------------------------------------------
      # Log-space %RMSE (same logic as Section 11)
      # -----------------------------------------------
      d$RMSE_log_pc <- ifelse(
        is.na(mean_obs_log) || abs(mean_obs_log) < 1e-12,
        NA_real_,
        (d$RMSE / abs(mean_obs_log)) * 100
      )
      
      # -----------------------------------------------
      # RMSE_ppm is already ppm for ACE
      # -----------------------------------------------
      d$RMSE_ppm <- d$RMSE
      
      d$RMSE_ppm_pc <- ifelse(
        is.na(d$RMSE_ppm) | is.na(d$RMSE_log_pc),
        NA_real_,
        d$RMSE_ppm * (d$RMSE_log_pc / 100)
      )
      
      # -----------------------------------------------
      # NEW: RMSEP metrics for ACE
      # -----------------------------------------------
      
      # RMSEP_ppm should equal RMSEP (ACE diagnostics RMSE = ppm units)
      if ("RMSEP" %in% names(d)) {
        d$RMSEP_ppm <- d$RMSEP
      } else {
        d$RMSEP_ppm <- NA_real_
      }
      
      # percentage RMSEP in ppm space
      d$RMSEP_ppm_pc <- ifelse(
        is.na(d$RMSEP_ppm) | is.na(d$RMSE_log_pc),
        NA_real_,
        d$RMSEP_ppm * (d$RMSE_log_pc / 100)
      )
      
      # Add element label & ranking
      d$Element <- el
      
      d <- d %>%
        dplyr::arrange(desc(R2), RMSE) %>%
        dplyr::mutate(Rank = dplyr::row_number())
      
      ace_global <- rbind(ace_global, d)
    }
    
    # -----------------------------------------------------------
    # PATCH: Compute Rank per element for ACE global summary
    # -----------------------------------------------------------
    ace_global <- ace_global %>%
      dplyr::group_by(Element) %>%
      dplyr::arrange(dplyr::desc(R2), RMSE, .by_group = TRUE) %>%
      dplyr::mutate(Rank = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(Rank, dplyr::everything())
    
    if (nrow(ace_global) > 0) {
      
      write.csv(ace_global, summary_file, row.names = FALSE)
      
      best_models <- ace_global %>% dplyr::filter(Rank == 1)
      write.csv(best_models, best_file, row.names = FALSE)
      
      message("=== ACE global summaries written ===")
    }
  }
  # ===============================================================
  # 18. ACE ppm-scale summary (UPDATED / PATCHED)
  # ===============================================================
  
  if (do_ACE_predictions && !is.null(ACE_dataset)) {
    
    message("=== ACE ppm-scale summary (patched Section 18) ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    summary_file <- file.path(ace_root, "AllElements_ModelSummary_ppm_pred.csv")
    
    # Must have ace_global from Section 17
    if (!exists("ace_global") || nrow(ace_global) == 0) {
      warning("ace_global is missing or empty — skipping ACE ppm summary.")
    } else {
      
      # Build mean ppm for each element
      elem_means_list <- lapply(elements, function(el) {
        obs_col <- paste0(el, "_ICP")
        if (!(obs_col %in% names(ICP_ppm))) {
          return(data.frame(Element = el, Mean_ppm = NA_real_))
        }
        vals <- ICP_ppm[[obs_col]]
        vals <- vals[!is.na(vals)]
        data.frame(
          Element  = el,
          Mean_ppm = ifelse(length(vals) == 0, NA_real_, mean(vals))
        )
      })
      
      elem_means_df <- do.call(rbind, elem_means_list)
      
      # Merge means into ACE global table
      ace_ppm_summary <- ace_global %>%
        dplyr::left_join(elem_means_df, by = "Element") %>%
        dplyr::mutate(
          
          # RMSE_ppm is already included in ace_global
          RMSE_ppm_pc = ifelse(
            is.na(Mean_ppm) | Mean_ppm == 0 | is.na(RMSE_ppm),
            NA_real_,
            (RMSE_ppm / Mean_ppm) * 100
          ),
          
          RMSEP_ppm_pc = ifelse(
            is.na(Mean_ppm) | Mean_ppm == 0 | is.na(RMSEP_ppm),
            NA_real_,
            (RMSEP_ppm / Mean_ppm) * 100
          )
        )
      
      write.csv(ace_ppm_summary, summary_file, row.names = FALSE)
      
      message("=== ACE ppm-scale summary written successfully ===")
    }
  }
  # ===============================================================
  # 19. ACE DIAGNOSTICS DASHBOARD & VALIDATOR 
  # ===============================================================
  
  if (do_ACE_predictions && exists("ace_global")) {
    
    message("=== Creating ACE diagnostics dashboard ===")
    
    ace_root <- file.path(all_base, "Calibration_pred")
    dashboard_file <- file.path(ace_root, "ACE_Diagnostics_Dashboard.pdf")
    
    suppressMessages({
      pdf(dashboard_file, width=29.7/2.54, height=21/2.54)
      
      # 1. Title page
      plot.new()
      title(main="ACE Diagnostics Dashboard",
            sub=paste("Generated on", Sys.Date()),
            cex.main=2)
      
      # 2. R2 ranking plot
      ace_global_clean <- ace_global %>%
        dplyr::filter(!is.na(R2))
      
      p_r2 <- ggplot(ace_global_clean,
                     aes(x=Element, y=R2, fill=Model)) +
        geom_col(position="dodge") +
        theme_bw() +
        ggtitle("ACE R² by Element and Model")
      
      print(p_r2)
      
      # 3. RMSE plot
      p_rmse <- ggplot(ace_global_clean,
                       aes(x=Element, y=RMSE, fill=Model)) +
        geom_col(position="dodge") +
        theme_bw() +
        ggtitle("ACE RMSE by Element and Model")
      
      print(p_rmse)
      
      # 4. Best model table
      best_models <- ace_global_clean %>% dplyr::filter(Rank==1)
      gridExtra::grid.table(best_models)
      
      dev.off()
    })
  }
  
  # ===============================================================
  # 19.1 VALIDATOR: Check all expected output files for Sections 8 & 13–18
  # ===============================================================
  
  validate_full_pipeline <- function(
    elements,
    all_base,
    calib_base = file.path(all_base, "Calibration_outputs"),
    ace_base   = file.path(all_base, "Calibration_pred"),
    p_po_list  = NULL,
    p_po_list_ace_orig = NULL,
    p_po_list_ace_pc   = NULL
  ) {
    
    message("\n=== RUNNING PIPELINE VALIDATION ===\n")
    
    results <- list()
    
    add_result <- function(category, file, exists) {
      results[[length(results) + 1]] <<- data.frame(
        Category = category,
        File     = file,
        Exists   = exists
      )
    }
    
    # -----------------------------------------------------------
    # 1. Calibration per-element outputs (Section 8)
    # -----------------------------------------------------------
    message("Checking Section 8 calibration outputs...")
    
    for (el in elements) {
      el_dir <- file.path(calib_base, el)
      
      # Expected files
      expected_files <- c(
        paste0(el, "_Model_Ranking.csv"),
        paste0(el, "_Diagnostics_Summary.csv"),
        paste0(el, "_preds.csv"),
        paste0(el, "_PredObs_Model_Summary_1page.pdf"),
        paste0(el, "_PredObs_Model_Summary_2page.pdf")
      )
      
      for (f in expected_files) {
        add_result(
          category = paste("Calibration:", el),
          file     = file.path(el_dir, f),
          exists   = file.exists(file.path(el_dir, f))
        )
      }
    }
    
    # -----------------------------------------------------------
    # 2. Global calibration summaries (Section 11)
    # -----------------------------------------------------------
    message("Checking Section 11 global summaries...")
    
    global_files <- c(
      "AllElements_ModelSummary.csv",
      "BestModels_PerElement.csv"
    )
    
    for (f in global_files) {
      add_result(
        category = "Global Calibration Summary",
        file     = file.path(all_base, f),
        exists   = file.exists(file.path(all_base, f))
      )
    }
    
    # -----------------------------------------------------------
    # 3. ACE predictions (Section 13)
    # -----------------------------------------------------------
    message("Checking ACE prediction folders...")
    
    for (el in elements) {
      el_dir <- file.path(ace_base, el)
      expected <- paste0(el, "_ACE_predictions.csv")
      add_result(
        category = paste("ACE Predictions:", el),
        file     = file.path(el_dir, expected),
        exists   = file.exists(file.path(el_dir, expected))
      )
    }
    
    # -----------------------------------------------------------
    # 4. ACE diagnostics (Section 14)
    # -----------------------------------------------------------
    message("Checking ACE diagnostics...")
    
    for (el in elements) {
      f <- file.path(ace_base, el, paste0(el, "_ACE_diagnostics.csv"))
      add_result(
        category = paste("ACE Diagnostics:", el),
        file     = f,
        exists   = file.exists(f)
      )
    }
    
    # -----------------------------------------------------------
    # 5. ACE diagnostic plots (Section 15)
    # -----------------------------------------------------------
    message("Checking ACE diagnostic plots (Section 15)...")
    
    plot_suffixes <- c(
      "_ACE_PredVsObs.pdf",
      "_ACE_PredVsObs_err.pdf",
      "_ACE_PredVsObs_pc.pdf",
      "_ACE_PredVsObs_pc_err.pdf",
      "_ACE_Residuals.pdf"
    )
    
    for (el in elements) {
      for (m in names(model_store[[el]])) {
        for (suf in plot_suffixes) {
          f <- file.path(ace_base, el, paste0(el, "_", m, suf))
          add_result(
            category = paste("ACE Plots:", el, m),
            file     = f,
            exists   = file.exists(f)
          )
        }
      }
    }
    
    # -----------------------------------------------------------
    # 6. Section 16 ACE summary panels
    # -----------------------------------------------------------
    message("Checking ACE summary panels (Section 16)...")
    
    sum_suffixes <- c(
      "_ACE_Model_Summary_original.pdf",
      "_ACE_Model_Summary_percentage.pdf"
    )
    
    for (el in elements) {
      for (suf in sum_suffixes) {
        f <- file.path(ace_base, el, paste0(el, suf))
        add_result(
          category = paste("ACE Summary:", el),
          file     = f,
          exists   = file.exists(f)
        )
      }
    }
    
    # -----------------------------------------------------------
    # 7. ACE global summaries (Section 17)
    # -----------------------------------------------------------
    message("Checking ACE global summaries (Section 17)...")
    
    ace_global_files <- c(
      "AllElements_ModelSummary_pred.csv",
      "BestModels_PerElement_pred.csv"
    )
    
    for (f in ace_global_files) {
      add_result(
        category = "ACE Global Summary",
        file     = file.path(ace_base, f),
        exists   = file.exists(file.path(ace_base, f))
      )
    }
    
    # -----------------------------------------------------------
    # 8. ACE ppm summary (Section 18)
    # -----------------------------------------------------------
    message("Checking ACE ppm summary (Section 18)...")
    
    ppm_file <- file.path(ace_base, "AllElements_ModelSummary_ppm_pred.csv")
    add_result(
      category = "ACE PPM Summary",
      file     = ppm_file,
      exists   = file.exists(ppm_file)
    )
    
    # -----------------------------------------------------------
    # 9. Check in-memory plot stores (Section 15 + 16 dependencies)
    # -----------------------------------------------------------
    message("Checking in-memory ACE plot stores…")
    
    for (el in elements) {
      for (m in names(model_store[[el]])) {
        add_result(
          category = paste("ACE Plot Store (original PI):", el, m),
          file     = paste0("p_po_list_ace_orig[['", el, "']][['", m, "']]"),
          exists   = !is.null(p_po_list_ace_orig[[el]][[m]])
        )
        add_result(
          category = paste("ACE Plot Store (percentage PI):", el, m),
          file     = paste0("p_po_list_ace_pc[['", el, "']][['", m, "']]"),
          exists   = !is.null(p_po_list_ace_pc[[el]][[m]])
        )
      }
    }
    
    # -----------------------------------------------------------
    # Final report
    # -----------------------------------------------------------
    results_df <- do.call(rbind, results)
    
    missing <- results_df %>% dplyr::filter(!Exists)
    
    message("\n=== VALIDATION COMPLETE ===")
    message("Total checks: ", nrow(results_df))
    message("Missing files: ", nrow(missing))
    
    if (nrow(missing) > 0) {
      message("\n--- Missing Outputs ---")
      print(missing)
    } else {
      message("\nAll expected outputs are present.")
    }
    
    return(results_df)
  }
  
  # ===============================================================
  # 20. TESTER — CHECK ALL OUTPUTS EXIST
  # ===============================================================
  
  test_results <- list()
  
  ace_root <- file.path(all_base, "Calibration_pred")
  
  for (el in elements) {
    
    el_dir <- file.path(ace_root, el)
    
    expected_files <- c(
      paste0(el, "_ACE_predictions.csv"),
      paste0(el, "_ACE_diagnostics.csv"),
      paste0(el, "_ACE_PredObs_Model_Summary_2page.pdf"),
      paste0(el, "_ACE_PredObs_Model_Summary_1page.pdf"),
      paste0(el, "_ACE_PredObs_Model_Summary_1page_x.pdf")
    )
    
    exists_vec <- file.exists(file.path(el_dir, expected_files))
    
    test_results[[el]] <- data.frame(
      File = expected_files,
      Exists = exists_vec
    )
  }
  
  # Print test results
  for (el in names(test_results)) {
    message("Element: ", el)
    print(test_results[[el]])
  }
  
  # Check global files
  global_files <- c(
    "AllElements_ModelSummary_pred.csv",
    "BestModels_PerElement_pred.csv",
    "AllElements_ModelSummary_ppm_pred.csv",
    "ACE_Diagnostics_Dashboard.pdf"
  )
  
  message("--- Global files ---")
  print(data.frame(
    File   = global_files,
    Exists = file.exists(file.path(ace_root, global_files))
  ))
  
  
  
  
  
  # ===============================================================
  # 21. MULTI-ELEMENT CORRELATION DIAGNOSTICS
  # ===============================================================
  
  message("=== Running multi-element correlation diagnostics ===")
  
  corr_dir <- file.path(all_base, "Correlation_diagnostics")
  dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------------------------------
  # 21.1 CORRELATION MATRIX OF RAW ICP LOG VALUES
  # -------------------------------------------------------------------
  
  # extract ICP log variables
  icp_cols <- paste0(elements, "_ICP")
  icp_present <- icp_cols[icp_cols %in% names(data)]
  
  if (length(icp_present) >= 2) {
    
    icp_df <- data[, icp_present, drop = FALSE]
    
    corr_mat <- cor(icp_df, use = "pairwise.complete.obs")
    
    # save CSV
    write.csv(
      corr_mat,
      file.path(corr_dir, "ICP_logspace_correlation_matrix.csv"),
      row.names = TRUE
    )
    
    # heatmap
    melted <- reshape2::melt(corr_mat)
    p_corr <- ggplot(melted, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
      scale_fill_gradient2(
        low = "blue", high = "red", mid = "white", midpoint = 0,
        name = "r"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)
      ) +
      labs(
        title = "ICP log-space correlation matrix",
        x = "", y = ""
      )
    
    ggsave(
      file.path(corr_dir, "ICP_logspace_correlation_heatmap.pdf"),
      p_corr, width = 16, height = 12, units = "cm"
    )
  }
  
  # -------------------------------------------------------------------
  # 21.2 CORRELATION BETWEEN MODEL PREDICTIONS (CALIBRATION DATA)
  # -------------------------------------------------------------------
  
  message("=== Producing element-wise prediction correlation plots ===")
  
  for (el in elements) {
    
    element_dir <- file.path(corr_dir, el)
    dir.create(element_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (!el %in% names(preds_store)) next
    df <- preds_store[[el]]
    
    # find prediction columns (ppm back-converted)
    pred_cols <- grep(paste0("_", el, "_Pred_ppm$"), names(df), value = TRUE)
    if (length(pred_cols) < 2) next
    
    preds_only <- df[, pred_cols, drop = FALSE]
    
    # correlation matrix
    corr_pred <- cor(preds_only, use = "pairwise.complete.obs")
    
    # save CSV
    write.csv(
      corr_pred,
      file.path(element_dir, paste0(el, "_prediction_correlation_matrix.csv")),
      row.names = TRUE
    )
    
    # heatmap
    melted_p <- reshape2::melt(corr_pred)
    
    p_pred_corr <- ggplot(melted_p, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
      scale_fill_gradient2(
        low = "blue", high = "red", mid = "white", midpoint = 0,
        name = "r"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title   = element_text(size = 10, face = "bold")
      ) +
      labs(
        title = paste("Model prediction correlation –", el),
        x = "", y = ""
      )
    
    ggsave(
      file.path(element_dir, paste0(el, "_prediction_correlation_heatmap.pdf")),
      p_pred_corr,
      width = 16, height = 12, units = "cm"
    )
  }
  
  # -------------------------------------------------------------------
  # 19.3 PAIRWISE SCATTERPLOTS OF ICP ELEMENTS
  # -------------------------------------------------------------------
  
  if (length(icp_present) >= 2) {
    
    pairs_dir <- file.path(corr_dir, "pairwise_scatterplots")
    dir.create(pairs_dir, recursive = TRUE, showWarnings = FALSE)
    
    icp_df2 <- data[, c("Site", icp_present), drop = FALSE]
    
    for (i in seq_along(icp_present)) {
      for (j in seq_along(icp_present)) {
        
        if (j <= i) next
        
        e1 <- icp_present[i]
        e2 <- icp_present[j]
        
        p_sc <- ggplot(icp_df2, aes_string(x = e1, y = e2, colour = "Site")) +
          geom_point(size = 1.5) +      # no transparency
          ggsci::scale_color_jco(name = "Site") +
          theme_bw() +
          theme(
            plot.title = element_text(size = 10, face = "bold")
          ) +
          labs(
            title = paste("ICP relationship:", e1, "vs", e2),
            x = e1, y = e2
          )
        
        ggsave(
          file.path(
            pairs_dir,
            paste0("ICP_scatter_", e1, "_vs_", e2, ".pdf")
          ),
          p_sc, width = 13, height = 10, units = "cm"
        )
      }
    }
  }
  
  message("=== Correlation diagnostics completed successfully ===")

  # ===============================================================
  # 22. WRITE GLOBAL README FILE (placed in Regression_Multivariate)
  # ===============================================================
  
  readme_contents <- c(
    "REGRESSION MULTIVARIATE OUTPUT STRUCTURE",
    "========================================",
    "",
      "Overview",
      "",
    "This workflow performs complete calibration, validation, and reconstruction for multiple geochemical elements using log-transformed ICP-MS calibration data and log-transformed XRF datasets.",
    "For each element, the function fits eight statistical models (univariate and multivariate where applicable):",
    "",
    "- Ordinary Least Squares (OLS)",
    "- Weighted Least Squares (WLS)",
    "- Weighted OLS (using standard deviations)",
    "- Two-step WLS",
    "- Bayesian Generalised Linear Model (BGLM)",
    "- Random Forest (RF)",
    "- Partial Least Squares (PLS) with Leave-One-Out (LOO) CV",
    "- PLS with k-fold CV",
    "",
    "Model performance is evaluated using:",
    "- R², RMSE, 10-fold RMSEP",
    "- AIC/BIC (where applicable)",
    "- Normality and heteroscedasticity tests",
    "",
    "The workflow generates predicted vs observed plots, residual and influence diagnostics, and creates model rankings based on R² → RMSEP → RMSE.",
    "It then produces ppm-scale predictions with 95% confidence intervals, site-level and multi-site/multi-element comparison plots (depth and age), and compiles global summary tables in both log and ppm space.",
    "",
    "Random Forest OOB (Out-of-Bag)",
    "",
    "Random Forest models include OOB error estimation, providing a robust measure of performance without a separate test set.",
    "Each tree predicts observations not included in its bootstrap sample, and OOB predictions are averaged to compute unbiased error and R² estimates.",
    "",
    "Models Included",
    "",
    "Single-predictor models:",
    "- OLS",
    "- WLS (variance-weighted OLS)",
    "- OLS_weighted",
    "- WLS_weighted",
    "",
    "Multi-predictor models (using all 6 elements as predictors):",
    "- Bayes",
    "- PLS_LOO",
    "- PLS_kfold",
    "- Random Forest",
    "",
    "Metrics, Tests, and Outputs",
    "",
    "Calculated metrics:",
    "- True R² = 1 - (SS_res / SS_tot)",
    "- RMSE",
    "- 10-fold RMSEP",
    "- RF OOB MSE and R²",
    "- PLS RMSEP at optimal ncomp",
    "",
    "Statistical tests:",
    "- AIC & BIC (where applicable)",
    "- Shapiro–Wilk normality test",
    "- Breusch–Pagan heteroscedasticity test (LM/Bayes)",
    "",
    "Plots (PDF):",
    "- Residual distribution plots",
    "- Predicted vs Observed plots (coloured by Site using a JCO palette)",
    "- Influence plots (Cook’s distance + leverage; LM models only)",
    "- PLS RMSEP validation plots",
    "- RF variable importance plots",
    "",
    "Per-element outputs:",
    "- Diagnostics CSV (AIC, BIC, tests, CV metrics, explanations)",
    "- Performance summary CSV (R², RMSE, RMSEP, AIC, BIC)",
    "- Best-model table",
    "",
    "Ranking Logic",
    "",
    "Models are ranked by:",
    "1. Highest R² (best explanatory fit)",
    "2. Lowest RMSEP (best predictive CV performance)",
    "3. Lowest RMSE",
    "",
    "No PLS priority and no AIC/BIC in the ranking criteria.",
    "",
    "Directory Output Structure Per Element:",
    "- *_Model_Ranking.csv",
    "- *_Diagnostics_Summary.csv",
    "- *_PredVsObs.pdf",
    "- *_Residuals.pdf",
    "- *_Influence.pdf (LM models only)",
    "- *_PLS_RMSEP.pdf",
    "- README.txt",
    "",
    "Global outputs in All_models & Best_models folders:",
    "- AllElements_ModelSummary.csv",
    "- BestModels_PerElement.csv",
    "",
    "Acknowledgements",
    "",
    "Based on original OLS, WLS, PLS regression code by Steve Roberts.",
    "Refined and made into a function with assistance from ChatGPT (OpenAI, version 5.1).",
    "We acknowledge the use of ChatGPT for assistance in constructing functions based on manual code written by the authors and generating some preliminary code.",
    "The final implementations of the code and functions were produced and verified independently by the authors.",
    "",    
    "This file describes all folders and files generated by run_full_regressions().",
    "",
    "-------------------------------------------------------------",
    "1. All_models_ppm/",
    "-------------------------------------------------------------",
    "Contains ALL model outputs for all 8 calibration models:",
    "  OLS, OLS_weighted, WLS, WLS_weighted, Bayes, RF, PLS(LOO), PLS(k).",
    "",
    "Subfolders:",
    "",
    "  Calibration_outputs/<Element>/",
    "    * <el>_Model_Ranking.csv",
    "    * <el>_Diagnostics_Summary.csv",
    "    * <el>_<model>_PredVsObs.pdf",
    "    * <el>_<model>_Residuals.pdf",
    "    * <el>_<model>_Influence.pdf",
    "    * <el>_<model>_PLS_RMSEP.pdf",
    "    * <el>_PredObs_Model_Summary_2Page.pdf    (8-model A4 2 pages, 2×4 layout)",
    "    * <el>_PredObs_Model_Summary_1Page.pdf    (8-model A4 1 page, 2×8 layout)",
    "    * <el>_preds.csv                    (ppm predictions for XRF_new)",
    "    * README.txt",
    "",
    "  Sites/",
    "    Depth and age profile reconstructions per site:",
    "    * Site_<site>_Profiles_depth_all.pdf",
    "    * Site_<site>_Profiles_depth_multivariate.pdf",
    "    * Site_<site>_Profiles_age_all.pdf",
    "    * Site_<site>_Profiles_age_multivariate.pdf",
    "",
    "  elements/",
    "    Cross-site age–concentration comparisons:",
    "    * <el>_model_age_comparison_all.pdf",
    "    * <el>_model_age_comparison_all_y_equal.pdf",
    "    * <el>_model_age_comparison_multivariate.pdf",
    "    * <el>_model_age_comparison_multivariate_y_equal.pdf",
    "",
    "Global model summaries:",
    "  * AllElements_ModelSummary.csv",
    "  * BestModels_PerElement.csv",
    "  * AllElements_ModelSummary_ppm.csv",
    "",
    "-------------------------------------------------------------",
    "2. Best_models_ppm/",
    "-------------------------------------------------------------",
    "Contains ONLY the best or multivariate models.",
    "Includes 'Sites/' and 'elements/' outputs using only chosen models.",
    "",
    "-------------------------------------------------------------",
    "3. Calibration_pred/",
    "-------------------------------------------------------------",
    "Predictions applied to ACE/ICP dataset (external predictions).",
    "",
    "Structure:",
    "  Calibration_pred/<Element>/",
    "    * <el>_ICP_pred.csv",
    "    * <el>_<model>_PredVsObs_ACE.pdf",
    "    * <el>_PredObs_Model_Summary_pred.pdf",
    "    * <el>_Model_Ranking_pred.csv",
    "    * <el>_Diagnostics_Summary_pred.csv",
    "Plots show predicted concentrations (Y-axis) versus observed ICP-MS values (X-axis); horizontal error bars represent ICP analytical uncertainty, vertical error bars represent model prediction uncertainty (95% CI).",
    "",
    "Global summaries:",
    "  * AllElements_ModelSummary_pred.csv",
    "  * BestModels_PerElement_pred.csv",
    "  * AllElements_ModelSummary_ppm_pred.csv",
    " Model robustness was evaluated using a composite signal-to-noise framework 
      incorporating prediction uncertainty, downcore smoothness, and calibration 
      performance. Signal-to-noise ratios were calculated from prediction confidence 
      intervals and profile roughness, scaled within each element, and combined with
      R² into a weighted robustness score overall and per site.
      Overall robustness across the whole calibration dataset was useful, but models
      exhibiting unstable downcore behaviour were assessed on a per site basis and
      excluded from selection if classified as unstable (robustness <0.75).
      For each element, and at each site, the highest-ranked
      stable model was designated as the production model and used for subsequent
      interpretation.

    Metrics used
    
    Quantity	Definition
    Signal 	  Variability of predicted concentrations
    Noise 	  Model uncertainty + downcore roughness
    Scale	    ppm (mg kg⁻¹), not log
    
    Metrics computed per element × model:
    1.	SNR_model = mean predicted ppm / RMSE_ppm
    2.	SNR_CI (preferred) = sd(predicted ppm) / mean(U95 - L95)
    3.	SNR_smooth (downcore stability) = sd(predicted ppm) / sd(diff(predicted ppm)
    
    For each element × model, robustness is derived from:
    1.	Accuracy
    •	R² (calibration or validation performance)
    2.	Noise / uncertainty
    •	SNR_CI
    Signal-to-noise ratio based on prediction uncertainty
    (signal variability ÷ mean CI width)
    3.	Downcore stability
    •	SNR_smooth
    Penalises high point-to-point roughness (spiky profiles)
    
    Robustness score
    
    Each metric is scaled 0–1 within each element, then combined:
    
    Robustness Score = 0.4 × SNR_CI_scaled + 0.3 × SNR_smooth_scaled + 0.3 × R²_scaled
    
    Models are ranked per element by:
    1.	Highest SNR_CI
    2.	Highest SNR_smooth
    3.	Highest R²
    
    This ensures:
    •	Models with narrow CIs,
    •	smooth downcore behaviour, and
    •	high explanatory power
    are favoured.
    
    Class	Definition
    Preferred:	Highest robustness score and not flagged unstable
    Acceptable:	Rank 2–3 robustness, stable but slightly noisier
    Unstable:	High roughness or poor SNR despite reasonable R²
    
    Unstable models are never selected as production models, even if R² is high.
    
    Generates a summary Heatmap matrix for visualising robustness quickly for all
    models vs elements and radar plots for each element ",
    "",
    "-------------------------------------------------------------",
    "4. Multi-element correlation diagnostics (Section 21)",
    "-------------------------------------------------------------",
    "May include correlation heatmaps, PCA, and pairwise plots.",
    "Stored in: All_models_ppm/Diagnostics/ or Calibration_pred/Diagnostics/",
    "",
    "-------------------------------------------------------------",
    "5. Logging",
    "-------------------------------------------------------------",
    "Warnings/messages optionally stored in logs.",
    "",
    "-------------------------------------------------------------",
    "6. Completion Message",
    "-------------------------------------------------------------",
    "=== COMPLETED SUCCESSFULLY ===",
    "",
    "All results saved in: <Regression_Multivariate>"
  )
  
  # Write to Regression_Multivariate/ReadMe.txt
  writeLines(
    readme_contents,
    file.path(out_dir, "ReadMe.txt")
  )
  
  message("README written to: ", file.path(out_dir, 'ReadMe.txt'))  
  
} # end run_full_regressions()
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Run it ...
# 6-elements "Ca", "Ti", "Fe", "Mn", "Sr", "Zr" ---------------------------
run_full_regressions(
  elements = c( "Ca", "Ti", "Fe", "Mn", "Sr", "Zr" ), # choose elements to include 
  data     = ACE_dataset,   # dataset for calibration model (log-space as log(element) ICP-MS & XRF-CS log(element/inc))
  XRF_new  = XRF_pred,      # predictors (log-space) for conversion to ppm across 8 models and for elements defined above
  ICP_ppm  = ICP_obs,       # matched ICP data in ppm: Ti_ICP, Ti_ICP_sd, Site, SH20_age, etc., (XRF in cps not used as an input)
  save_dir = "/Users/sjro/Desktop"  # function will create /Desktop/Regression_Multivariate - copy this to output folder elsewhere
)

# 4-elements "Ca", "Ti", "Sr", "Zr" ---------------------------------------
run_full_regressions(
  elements = c( "Ca", "Ti", "Sr", "Zr" ), # choose elements to include 
  data     = ACE_dataset,   # dataset for calibration model (log-space as log(element) ICP-MS & XRF-CS log(element/inc))
  XRF_new  = XRF_pred,      # predictors (log-space) for conversion to ppm across 8 models and for elements defined above
  ICP_ppm  = ICP_obs,       # matched ICP data in ppm: Ti_ICP, Ti_ICP_sd, Site, SH20_age, etc., (XRF in cps not used as an input)
  save_dir = "/Users/sjro/Desktop"  # function will create /Desktop/Regression_Multivariate - copy this to output folder elsewhere
)
# -------------------------------------------------------------------------


