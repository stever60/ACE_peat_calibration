# Figure 3b - Linear models, ordinary and weighted least squares regression analysis (OLS, WLS)

# clr 

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
              'directlabels', 'performance', 'lmtest', 'ggpmisc', 'cowplot', 'Hmisc')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

# Set up ------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()

# Define elements to use ----------------------------------------------------------

# ICP elements as defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")

# Below were defined by autocorrelation (acf) analysis of ACE09 composite ITRAX dataframe analysis file - see other R files for mor info

# XRF elements defined by ITRAX acf and matched to Francois ICPMS element list above
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

# ICP elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# Import existing ACE cps matched file from Fig 2 - make clr dataset ----------------------------------------------

ACE_xrf_icp_matched <-read_csv("Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/clr/ACE_matched_xrf_icp_cps.csv")
is.na(ACE_xrf_icp_matched)<-sapply(ACE_xrf_icp_matched, is.infinite) # replace any infinite values with NA

# Replace zeros with half minimum value for each element 
ACE_all0 <- ACE_xrf_icp_matched %>% 
  mutate(across(all_of(c(icp_Elements_min, icp_Elements_min_sd)), ~ ifelse(.x < 0, 0, .x))) %>%  # replace any ICPMS values <0 with zero
  #zeroreplace() %>% # replace zeros as outlined in the compositions package
  mutate_at(vars(all_of(c(icp_Elements_min, acf_icp_Elements_min))), ## replace zeros with half minimum value to allow linear modelling to work
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% # Recommended procedure from Bertrand et al. (submitted) - retains dataframe structure
  mutate_at(vars(all_of(c(icp_Elements_min_sd, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% 
  select(Location:MSE, all_of(acf_icp_Elements_min), coh_inc, all_of(acf_icp_Elements_min_sd), coh_inc_sd,
         all_of(icp_Elements_min), dry_mass_pc, all_of(icp_Elements_min_sd), dry_mass_err) %>% 
  na.omit() %>% #remove rows with NAs - in this case there is only one at the end of CAT1-S1-1F
  print()

# Make clr file of while dataset after replacing zeros with half min value and using using composition package 
ACE_all_clr_itrax <- ACE_all0 %>% 
  select(all_of(acf_icp_Elements_min)) %>%
  select(-c(Mo_inc, Mo_coh)) %>% 
  clr() %>% # make centred log ratio dataframe
  as_tibble()

ACE_all_clr_itrax_sd <- ACE_all_clr_itrax %>% 
  mutate(across(where(is.numeric), ~ .x * 0.05)) %>% 
  rename_with(.fn = function(.x){paste0(.x,"_sd")})

ACE_all_clr_icp <- ACE_all0 %>% 
  select(all_of(icp_Elements_min)) %>%
  clr() %>% # make centred log ratio dataframe
  as_tibble()

ACE_all_clr_icp_sd <- ACE_all_clr_icp %>% 
  mutate(across(where(is.numeric), ~ .x * 0.05)) %>% 
  rename_with(.fn = function(.x){paste0(.x,"_sd")})

# Make other parameters / variables 
ACE_all_text <- ACE_all0 %>% 
  select(Location:MSE)
ACE_all_rest <- ACE_all0 %>% 
  select(Mo_inc, Mo_coh, Mo_inc_sd, Mo_coh_sd, coh_inc, coh_inc_sd,
         dry_mass_pc, dry_mass_err)

# Bind into single file & define dataset to use for linear modelling
ACE_LM1 <- bind_cols(ACE_all_text, ACE_all_clr_itrax, ACE_all_clr_itrax_sd, 
                     ACE_all_clr_icp, ACE_all_clr_icp_sd, ACE_all_rest) %>%
  filter(!Site =="POB4") %>% #remove POB4 data
  #calculate error in clr ICPMS
  rename(K_ICP_sd_FDV2023 = K_ICP_sd) %>% 
  mutate(K_ICP_sd = K_ICP*0.18) %>%
  relocate(K_ICP_sd, .after = K_ICP) %>%
  rename(Ca_ICP_sd_FDV2023 = Ca_ICP_sd) %>% 
  mutate(Ca_ICP_sd = Ca_ICP*0.13) %>% 
  relocate(Ca_ICP_sd, .after = Ca_ICP) %>%
  rename(Ti_ICP_sd_FDV2023 = Ti_ICP_sd) %>% 
  mutate(Ti_ICP_sd = Ti_ICP*0.14) %>%
  relocate(Ti_ICP_sd, .after = Ti_ICP) %>%
  rename(Mn_ICP_sd_FDV2023 = Mn_ICP_sd) %>% 
  mutate(Mn_ICP_sd = Mn_ICP*0.12) %>% 
  relocate(Mn_ICP_sd, .after = Mn_ICP) %>%
  rename(Fe_ICP_sd_FDV2023 = Fe_ICP_sd) %>% 
  mutate(Fe_ICP_sd = Fe_ICP*0.18) %>%
  relocate(Fe_ICP_sd, .after = Fe_ICP) %>%
  rename(Co_ICP_sd_FDV2023 = Co_ICP_sd) %>% 
  mutate(Co_ICP_sd = Co_ICP*0.07) %>%
  relocate(Co_ICP_sd, .after = Co_ICP) %>%
  rename(Ni_ICP_sd_FDV2023 = Ni_ICP_sd) %>% 
  mutate(Ni_ICP_sd = Ni_ICP*0.54) %>%
  relocate(Ni_ICP_sd, .after = Ni_ICP) %>%
  rename(Cu_ICP_sd_FDV2023 = Cu_ICP_sd) %>% 
  mutate(Cu_ICP_sd = Cu_ICP*0.14) %>%
  relocate(Cu_ICP_sd, .after = Cu_ICP) %>%
  rename(Zn_ICP_sd_FDV2023 = Zn_ICP_sd) %>% 
  mutate(Zn_ICP_sd = Zn_ICP*0.09) %>% 
  relocate(Zn_ICP_sd, .after = Zn_ICP) %>%
  rename(Rb_ICP_sd_FDV2023 = Rb_ICP_sd) %>% 
  mutate(Rb_ICP_sd = Rb_ICP*0.16) %>%
  relocate(Rb_ICP_sd, .after = Rb_ICP) %>%
  rename(Sr_ICP_sd_FDV2023 = Sr_ICP_sd) %>% 
  mutate(Sr_ICP_sd = Sr_ICP*0.02) %>% 
  relocate(Sr_ICP_sd, .after = Sr_ICP) %>% 
  rename(Zr_ICP_sd_FDV2023 = Zr_ICP_sd) %>% 
  mutate(Zr_ICP_sd = Zr_ICP*0.05) %>% 
  relocate(Zr_ICP_sd, .after = Zr_ICP) %>% 
  #calculate error in clr ITRAX data
  mutate(K_sd = K*0.1) %>%
  relocate(K_sd, .after = K) %>%
  mutate(Ca_sd = Ca*0.1) %>% 
  relocate(Ca_sd, .after = Ca) %>%
  mutate(Ti_sd = Ti*0.1) %>%
  relocate(Ti_sd, .after = Ti) %>%
  mutate(Mn_sd = Mn*0.1) %>% 
  relocate(Mn_sd, .after = Mn) %>%
  mutate(Fe_sd = Fe*0.1) %>%
  relocate(Fe_sd, .after = Fe) %>%
  mutate(Co_sd = Co*0.1) %>%
  relocate(Co_sd, .after = Co) %>%
  mutate(Ni_sd = Ni*0.1) %>%
  relocate(Ni_sd, .after = Ni) %>%
  mutate(Cu_sd = Cu*0.1) %>%
  relocate(Cu_sd, .after = Cu) %>%
  mutate(Zn_sd = Zn*0.1) %>% 
  relocate(Zn_sd, .after = Zn) %>%
  mutate(Rb_sd = Rb*0.1) %>%
  relocate(Rb_sd, .after = Rb) %>%
  mutate(Sr_sd = Sr*0.1) %>% 
  relocate(Sr_sd, .after = Sr) %>% 
  mutate(Zr_sd = Zr*0.1) %>% 
  relocate(Zr_sd, .after = Zr)
write.csv(ACE_LM1,"Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/clr/ACE_xrf_icp_matched_noPOB4_clr.csv", row.names = FALSE)
ACE_LM1

ACE_LM1_stats <- ACE_LM1 %>%
  select(kcps: dry_mass_err) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_LM1_stats,"Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/clr/ACE_LM1_clr_stats.csv", row.names = FALSE)

# Convert Site to use as a grouping variable
ACE_LM1$Site <- as.factor(ACE_LM1$Site)

# OLS & WLS Linear modelling -----------------------------------------------------------

# Load required packages
library(performance) # linear model testing and graphical outputs
library(ggpmisc) # this package is used for simple/quick for labelling (stat_poly_eq)
library(lmtest) # linear model testing 
library(car) # Leverage

# Linear models & performance assessment /stats -----------------------------------------------

# After first run, change input to e.g., ACE_K_OLS_wt_no_outliers.csv 
# or e.g., ACE_K_WLS_wt_no_outliers.csv & skip to Save stats & 
# stats tests plots to file section & rerun each element & rerun

# Set up labels
ACE_dataset <- ACE_LM1 
site_title <- "ACE"
itrax_dataset <- " clr"
icp_dataset <- " [clr]"
ACE_dataset$Site <- as.factor(ACE_dataset$Site) # Convert Site as a grouping variable

# Primary Matched Elements shown in Figure 3 ----------------------------------------------------
# ACE_Ca LM -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ca_lm <- lm(Ca_ICP ~ Ca, data = ACE_dataset)
summary(ACE_Ca_lm)
glance(ACE_Ca_lm)
model_performance(ACE_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Ca_wlm <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/(Ca_ICP_sd)^2)
summary(ACE_Ca_wlm)
glance(ACE_Ca_wlm)
model_performance(ACE_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Ca_model <- lm(Ca_ICP ~ Ca, data = ACE_dataset) # define model
ACE_Ca_wt <- 1 / lm(abs(ACE_Ca_model$residuals) ~ ACE_Ca_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ca_wls <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ca_wls) # summary stats
glance(ACE_Ca_wls) # summary stats including AIC
model_performance(ACE_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Ca_model_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/Ca_ICP_sd^2) # define model
ACE_Ca_wt_wt <- 1 / lm(abs(ACE_Ca_model_wt$residuals) ~ ACE_Ca_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ca_wls_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ca_wls_wt) # summary stats
glance(ACE_Ca_wls_wt) # summary stats including AIC
model_performance(ACE_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_lm_hats <- as.data.frame(hatvalues(ACE_Ca_lm))
ACE_Ca_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_lm_cooksD <- cooks.distance(ACE_Ca_lm)
ACE_Ca_lm_influential <- ACE_Ca_lm_cooksD[(ACE_Ca_lm_cooksD > (3 * mean(ACE_Ca_lm_cooksD, na.rm = TRUE)))]
ACE_Ca_lm_influential
ACE_Ca_lm_influential_names <- names(ACE_Ca_lm_influential)
ACE_Ca_lm_outliers <- ACE_dataset[ACE_Ca_lm_influential_names,] # outliers only using of index values
ACE_Ca_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ca_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_wlm_hats <- as.data.frame(hatvalues(ACE_Ca_wlm))
ACE_Ca_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_wlm_cooksD <- cooks.distance(ACE_Ca_wlm)
ACE_Ca_wlm_influential <- ACE_Ca_wlm_cooksD[(ACE_Ca_wlm_cooksD > (3 * mean(ACE_Ca_wlm_cooksD, na.rm = TRUE)))]
ACE_Ca_wlm_influential
ACE_Ca_wlm_influential_names <- names(ACE_Ca_wlm_influential)
ACE_Ca_wlm_outliers <- ACE_dataset[ACE_Ca_wlm_influential_names,] # outliers only using of index values
ACE_Ca_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ca_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_wls_hats <- as.data.frame(hatvalues(ACE_Ca_wls))
ACE_Ca_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_wls_cooksD <- cooks.distance(ACE_Ca_wls)
ACE_Ca_wls_influential <- ACE_Ca_wls_cooksD[(ACE_Ca_wls_cooksD > (3 * mean(ACE_Ca_wls_cooksD, na.rm = TRUE)))]
ACE_Ca_wls_influential
ACE_Ca_wls_influential_names <- names(ACE_Ca_wls_influential)
ACE_Ca_wls_outliers <- ACE_dataset[ACE_Ca_wls_influential_names,] # outliers only using of index values
ACE_Ca_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ca_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_wls_wt_hats <- as.data.frame(hatvalues(ACE_Ca_wls_wt))
ACE_Ca_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_wls_wt_cooksD <- cooks.distance(ACE_Ca_wls_wt)
ACE_Ca_wls_wt_influential <- ACE_Ca_wls_wt_cooksD[(ACE_Ca_wls_wt_cooksD > (3 * mean(ACE_Ca_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Ca_wls_wt_influential
ACE_Ca_wls_wt_influential_names <- names(ACE_Ca_wls_wt_influential)
ACE_Ca_wls_wt_outliers <- ACE_dataset[ACE_Ca_wls_wt_influential_names,] # outliers only using of index values
ACE_Ca_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ca_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_OLS_summary.txt")
summary(ACE_Ca_lm)
glance(ACE_Ca_lm)
model_performance(ACE_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_lm_lev.pdf")
leveragePlots(ACE_Ca_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_OLS_wt_summary.txt")
summary(ACE_Ca_wlm)
glance(ACE_Ca_wlm)
model_performance(ACE_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wlm_lev.pdf")
leveragePlots(ACE_Ca_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_WLS_summary.txt")
summary(ACE_Ca_wls)
glance(ACE_Ca_wls)
model_performance(ACE_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wls_lev.pdf")
leveragePlots(ACE_Ca_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_WLS_wt_summary.txt")
summary(ACE_Ca_wls_wt)
glance(ACE_Ca_wls_wt)
model_performance(ACE_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wls_wt_lev.pdf")
leveragePlots(ACE_Ca_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/ACE_Ca_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Ca"
Ca_WLS_wt = ACE_Ca_wt
Ca_WLS_err_wt = ACE_Ca_wt_wt
theme_set(theme_classic(10))

ACE_Ca <- ggplot(ACE_dataset, aes(x = Ca, y = Ca_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Ca_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Ca_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Ca_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Ca



# Define p value, equation & R2 as a string to add to plots

ACE_Ca_lm_p <- function(ACE_Ca_lm) {
  f <- summary(ACE_Ca_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_lm_p(ACE_Ca_lm)

ACE_Ca_lm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_lm_p(ACE_Ca_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Ca_wlm_p <- function(ACE_Ca_wlm) {
  f <- summary(ACE_Ca_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wlm_p(ACE_Ca_wlm)

ACE_Ca_wlm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wlm_p(ACE_Ca_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Ca_wlm_p <- function(ACE_Ca_wlm) {
  f <- summary(ACE_Ca_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wlm_p(ACE_Ca_wlm)

ACE_Ca_wlm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/(Ca_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wlm_p(ACE_Ca_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Ca_wls_p <- function(ACE_Ca_wls) {
  f <- summary(ACE_Ca_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wls_p(ACE_Ca_wls)

ACE_Ca_wls_eqn <- function(df){
  m <- ACE_Ca_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wls_p(ACE_Ca_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Ca_wls_wt_p <- function(ACE_Ca_wls_wt) {
  f <- summary(ACE_Ca_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wls_wt_p(ACE_Ca_wls_wt)

ACE_Ca_wls_wt_eqn <- function(df){
  m <- ACE_Ca_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wls_wt_p(ACE_Ca_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Ca_final <- ACE_Ca + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ca_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ca_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ca_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ca_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ca_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Ca_x.reg <- ACE_dataset$Ca
Ca_y.reg <- ACE_dataset$Ca_ICP
# Build linear regression model
Ca_model1 <- lm(Ca_y.reg ~ Ca, data = ACE_dataset) # OSL unweighted
Ca_model2 <- lm(Ca_y.reg ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt) # WLS unweighted
# Model summary statistics
Ca_model1
summary(Ca_model1)
confint(Ca_model1)
Ca_model2
summary(Ca_model2)
confint(Ca_model2)
# Create predicted values and upper/lower CI to check model is working
Ca_new1.y <- data.frame(Ca_x.reg = c(1, 2, 3))
predict(Ca_model1, Ca_newdata1 = Ca_new1.y)
predict(Ca_model1, Ca_newdata1 = Ca_new1.y, Ca_interval1 = "confidence")
Ca_new2.y <- data.frame(Ca_x.reg = c(1, 2, 3))
predict(Ca_model2, Ca_newdata2 = Ca_new2.y)
predict(Ca_model2, Ca_newdata2 = Ca_new2.y, Ca_interval2 = "confidence")
# Add prediction intervals to model data frame
Ca_pred.int1 <- predict(Ca_model1, interval = "prediction")
Ca_data_1_out <- bind_cols(ACE_dataset, Ca_pred.int1) %>%
  rename(Ca_fit_OLS = fit, Ca_lwr_OLS = lwr, Ca_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Ca, Ca_sd, Ca_ICP, Ca_ICP_sd, Ca_fit_OLS, Ca_lwr_OLS, Ca_upr_OLS))
Ca_pred.int2 <- predict(Ca_model2, interval = "prediction")
Ca_data_2_out <- bind_cols(Ca_data_1_out, Ca_pred.int2) %>%
  rename(Ca_fit_WLS = fit, Ca_lwr_WLS = lwr, Ca_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Ca, Ca_sd, Ca_ICP, Ca_ICP_sd, Ca_fit_OLS, Ca_lwr_OLS, Ca_upr_OLS, Ca_fit_WLS, Ca_lwr_WLS, Ca_upr_WLS)) %>% 
  print()
write.csv(Ca_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Ca_predict <- ACE_Ca_final + 
  geom_line(data = Ca_data_1_out, aes(y = Ca_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Ca_data_1_out, aes(y = Ca_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Ca_data_2_out, aes(y = Ca_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Ca_data_2_out, aes(y = Ca_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Ca_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ca/Ca_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")



 
# ACE_Ti LM -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ti_lm <- lm(Ti_ICP ~ Ti, data = ACE_dataset)
summary(ACE_Ti_lm)
glance(ACE_Ti_lm)
model_performance(ACE_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Ti_wlm <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/(Ti_ICP_sd)^2)
summary(ACE_Ti_wlm)
glance(ACE_Ti_wlm)
model_performance(ACE_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Ti_model <- lm(Ti_ICP ~ Ti, data = ACE_dataset) # define model
ACE_Ti_wt <- 1 / lm(abs(ACE_Ti_model$residuals) ~ ACE_Ti_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ti_wls <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ti_wls) # summary stats
glance(ACE_Ti_wls) # summary stats including AIC
model_performance(ACE_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Ti_model_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/Ti_ICP_sd^2) # define model
ACE_Ti_wt_wt <- 1 / lm(abs(ACE_Ti_model_wt$residuals) ~ ACE_Ti_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ti_wls_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ti_wls_wt) # summary stats
glance(ACE_Ti_wls_wt) # summary stats including AIC
model_performance(ACE_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_lm_hats <- as.data.frame(hatvalues(ACE_Ti_lm))
ACE_Ti_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_lm_cooksD <- cooks.distance(ACE_Ti_lm)
ACE_Ti_lm_influential <- ACE_Ti_lm_cooksD[(ACE_Ti_lm_cooksD > (3 * mean(ACE_Ti_lm_cooksD, na.rm = TRUE)))]
ACE_Ti_lm_influential
ACE_Ti_lm_influential_names <- names(ACE_Ti_lm_influential)
ACE_Ti_lm_outliers <- ACE_dataset[ACE_Ti_lm_influential_names,] # outliers only using of index values
ACE_Ti_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ti_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_wlm_hats <- as.data.frame(hatvalues(ACE_Ti_wlm))
ACE_Ti_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_wlm_cooksD <- cooks.distance(ACE_Ti_wlm)
ACE_Ti_wlm_influential <- ACE_Ti_wlm_cooksD[(ACE_Ti_wlm_cooksD > (3 * mean(ACE_Ti_wlm_cooksD, na.rm = TRUE)))]
ACE_Ti_wlm_influential
ACE_Ti_wlm_influential_names <- names(ACE_Ti_wlm_influential)
ACE_Ti_wlm_outliers <- ACE_dataset[ACE_Ti_wlm_influential_names,] # outliers only using of index values
ACE_Ti_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ti_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_wls_hats <- as.data.frame(hatvalues(ACE_Ti_wls))
ACE_Ti_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_wls_cooksD <- cooks.distance(ACE_Ti_wls)
ACE_Ti_wls_influential <- ACE_Ti_wls_cooksD[(ACE_Ti_wls_cooksD > (3 * mean(ACE_Ti_wls_cooksD, na.rm = TRUE)))]
ACE_Ti_wls_influential
ACE_Ti_wls_influential_names <- names(ACE_Ti_wls_influential)
ACE_Ti_wls_outliers <- ACE_dataset[ACE_Ti_wls_influential_names,] # outliers only using of index values
ACE_Ti_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ti_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_wls_wt_hats <- as.data.frame(hatvalues(ACE_Ti_wls_wt))
ACE_Ti_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_wls_wt_cooksD <- cooks.distance(ACE_Ti_wls_wt)
ACE_Ti_wls_wt_influential <- ACE_Ti_wls_wt_cooksD[(ACE_Ti_wls_wt_cooksD > (3 * mean(ACE_Ti_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Ti_wls_wt_influential
ACE_Ti_wls_wt_influential_names <- names(ACE_Ti_wls_wt_influential)
ACE_Ti_wls_wt_outliers <- ACE_dataset[ACE_Ti_wls_wt_influential_names,] # outliers only using of index values
ACE_Ti_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ti_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_OLS_summary.txt")
summary(ACE_Ti_lm)
glance(ACE_Ti_lm)
model_performance(ACE_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_lm_lev.pdf")
leveragePlots(ACE_Ti_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_OLS_wt_summary.txt")
summary(ACE_Ti_wlm)
glance(ACE_Ti_wlm)
model_performance(ACE_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wlm_lev.pdf")
leveragePlots(ACE_Ti_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_WLS_summary.txt")
summary(ACE_Ti_wls)
glance(ACE_Ti_wls)
model_performance(ACE_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wls_lev.pdf")
leveragePlots(ACE_Ti_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_WLS_wt_summary.txt")
summary(ACE_Ti_wls_wt)
glance(ACE_Ti_wls_wt)
model_performance(ACE_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wls_wt_lev.pdf")
leveragePlots(ACE_Ti_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/ACE_Ti_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Ti"
Ti_WLS_wt = ACE_Ti_wt
Ti_WLS_err_wt = ACE_Ti_wt_wt
theme_set(theme_classic(10))

ACE_Ti <- ggplot(ACE_dataset, aes(x = Ti, y = Ti_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Ti_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Ti_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Ti_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Ti

# Define p value, equation & R2 as a string to add to plots
 
ACE_Ti_lm_p <- function(ACE_Ti_lm) {
  f <- summary(ACE_Ti_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_lm_p(ACE_Ti_lm)

ACE_Ti_lm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_lm_p(ACE_Ti_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Ti_wlm_p <- function(ACE_Ti_wlm) {
  f <- summary(ACE_Ti_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wlm_p(ACE_Ti_wlm)

ACE_Ti_wlm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wlm_p(ACE_Ti_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Ti_wlm_p <- function(ACE_Ti_wlm) {
  f <- summary(ACE_Ti_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wlm_p(ACE_Ti_wlm)

ACE_Ti_wlm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/(Ti_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wlm_p(ACE_Ti_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Ti_wls_p <- function(ACE_Ti_wls) {
  f <- summary(ACE_Ti_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wls_p(ACE_Ti_wls)

ACE_Ti_wls_eqn <- function(df){
  m <- ACE_Ti_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wls_p(ACE_Ti_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Ti_wls_wt_p <- function(ACE_Ti_wls_wt) {
  f <- summary(ACE_Ti_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wls_wt_p(ACE_Ti_wls_wt)

ACE_Ti_wls_wt_eqn <- function(df){
  m <- ACE_Ti_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wls_wt_p(ACE_Ti_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Ti_final <- ACE_Ti + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ti_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ti_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ti_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ti_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ti_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Ti_x.reg <- ACE_dataset$Ti
Ti_y.reg <- ACE_dataset$Ti_ICP
# Build linear regression model
Ti_model1 <- lm(Ti_y.reg ~ Ti, data = ACE_dataset) # OSL unweighted
Ti_model2 <- lm(Ti_y.reg ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt) # WLS unweighted
# Model summary statistics
Ti_model1
summary(Ti_model1)
confint(Ti_model1)
Ti_model2
summary(Ti_model2)
confint(Ti_model2)
# Create predicted values and upper/lower CI to check model is working
Ti_new1.y <- data.frame(Ti_x.reg = c(1, 2, 3))
predict(Ti_model1, Ti_newdata1 = Ti_new1.y)
predict(Ti_model1, Ti_newdata1 = Ti_new1.y, Ti_interval1 = "confidence")
Ti_new2.y <- data.frame(Ti_x.reg = c(1, 2, 3))
predict(Ti_model2, Ti_newdata2 = Ti_new2.y)
predict(Ti_model2, Ti_newdata2 = Ti_new2.y, Ti_interval2 = "confidence")
# Add prediction intervals to model data frame
Ti_pred.int1 <- predict(Ti_model1, interval = "prediction")
Ti_data_1_out <- bind_cols(ACE_dataset, Ti_pred.int1) %>%
  rename(Ti_fit_OLS = fit, Ti_lwr_OLS = lwr, Ti_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Ti, Ti_sd, Ti_ICP, Ti_ICP_sd, Ti_fit_OLS, Ti_lwr_OLS, Ti_upr_OLS))
Ti_pred.int2 <- predict(Ti_model2, interval = "prediction")
Ti_data_2_out <- bind_cols(Ti_data_1_out, Ti_pred.int2) %>%
  rename(Ti_fit_WLS = fit, Ti_lwr_WLS = lwr, Ti_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Ti, Ti_sd, Ti_ICP, Ti_ICP_sd, Ti_fit_OLS, Ti_lwr_OLS, Ti_upr_OLS, Ti_fit_WLS, Ti_lwr_WLS, Ti_upr_WLS)) %>% 
  print()
write.csv(Ti_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Ti_predict <- ACE_Ti_final + 
  geom_line(data = Ti_data_1_out, aes(y = Ti_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Ti_data_1_out, aes(y = Ti_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Ti_data_2_out, aes(y = Ti_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Ti_data_2_out, aes(y = Ti_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Ti_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ti/Ti_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# ACE_Mn LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Mn_lm <- lm(Mn_ICP ~ Mn, data = ACE_dataset)
summary(ACE_Mn_lm)
glance(ACE_Mn_lm)
model_performance(ACE_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Mn_wlm <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/(Mn_ICP_sd)^2)
summary(ACE_Mn_wlm)
glance(ACE_Mn_wlm)
model_performance(ACE_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Mn_model <- lm(Mn_ICP ~ Mn, data = ACE_dataset) # define model
ACE_Mn_wt <- 1 / lm(abs(ACE_Mn_model$residuals) ~ ACE_Mn_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Mn_wls <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt) #perform weighted least squares regression
# Checks
summary(ACE_Mn_wls) # summary stats
glance(ACE_Mn_wls) # summary stats including AIC
model_performance(ACE_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Mn_model_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/Mn_ICP_sd^2) # define model
ACE_Mn_wt_wt <- 1 / lm(abs(ACE_Mn_model_wt$residuals) ~ ACE_Mn_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Mn_wls_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Mn_wls_wt) # summary stats
glance(ACE_Mn_wls_wt) # summary stats including AIC
model_performance(ACE_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_lm_hats <- as.data.frame(hatvalues(ACE_Mn_lm))
ACE_Mn_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_lm_cooksD <- cooks.distance(ACE_Mn_lm)
ACE_Mn_lm_influential <- ACE_Mn_lm_cooksD[(ACE_Mn_lm_cooksD > (3 * mean(ACE_Mn_lm_cooksD, na.rm = TRUE)))]
ACE_Mn_lm_influential
ACE_Mn_lm_influential_names <- names(ACE_Mn_lm_influential)
ACE_Mn_lm_outliers <- ACE_dataset[ACE_Mn_lm_influential_names,] # outliers only using of index values
ACE_Mn_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Mn_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_wlm_hats <- as.data.frame(hatvalues(ACE_Mn_wlm))
ACE_Mn_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_wlm_cooksD <- cooks.distance(ACE_Mn_wlm)
ACE_Mn_wlm_influential <- ACE_Mn_wlm_cooksD[(ACE_Mn_wlm_cooksD > (3 * mean(ACE_Mn_wlm_cooksD, na.rm = TRUE)))]
ACE_Mn_wlm_influential
ACE_Mn_wlm_influential_names <- names(ACE_Mn_wlm_influential)
ACE_Mn_wlm_outliers <- ACE_dataset[ACE_Mn_wlm_influential_names,] # outliers only using of index values
ACE_Mn_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Mn_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_wls_hats <- as.data.frame(hatvalues(ACE_Mn_wls))
ACE_Mn_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_wls_cooksD <- cooks.distance(ACE_Mn_wls)
ACE_Mn_wls_influential <- ACE_Mn_wls_cooksD[(ACE_Mn_wls_cooksD > (3 * mean(ACE_Mn_wls_cooksD, na.rm = TRUE)))]
ACE_Mn_wls_influential
ACE_Mn_wls_influential_names <- names(ACE_Mn_wls_influential)
ACE_Mn_wls_outliers <- ACE_dataset[ACE_Mn_wls_influential_names,] # outliers only using of index values
ACE_Mn_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Mn_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_wls_wt_hats <- as.data.frame(hatvalues(ACE_Mn_wls_wt))
ACE_Mn_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_wls_wt_cooksD <- cooks.distance(ACE_Mn_wls_wt)
ACE_Mn_wls_wt_influential <- ACE_Mn_wls_wt_cooksD[(ACE_Mn_wls_wt_cooksD > (3 * mean(ACE_Mn_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Mn_wls_wt_influential
ACE_Mn_wls_wt_influential_names <- names(ACE_Mn_wls_wt_influential)
ACE_Mn_wls_wt_outliers <- ACE_dataset[ACE_Mn_wls_wt_influential_names,] # outliers only using of index values
ACE_Mn_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Mn_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_OLS_summary.txt")
summary(ACE_Mn_lm)
glance(ACE_Mn_lm)
model_performance(ACE_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_lm_lev.pdf")
leveragePlots(ACE_Mn_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_OLS_wt_summary.txt")
summary(ACE_Mn_wlm)
glance(ACE_Mn_wlm)
model_performance(ACE_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wlm_lev.pdf")
leveragePlots(ACE_Mn_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_WLS_summary.txt")
summary(ACE_Mn_wls)
glance(ACE_Mn_wls)
model_performance(ACE_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wls_lev.pdf")
leveragePlots(ACE_Mn_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_WLS_wt_summary.txt")
summary(ACE_Mn_wls_wt)
glance(ACE_Mn_wls_wt)
model_performance(ACE_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wls_wt_lev.pdf")
leveragePlots(ACE_Mn_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/ACE_Mn_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Mn"
Mn_WLS_wt = ACE_Mn_wt
Mn_WLS_err_wt = ACE_Mn_wt_wt
theme_set(theme_classic(10))

ACE_Mn <- ggplot(ACE_dataset, aes(x = Mn, y = Mn_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Mn_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Mn_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Mn_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Mn

# Define p value, equation & R2 as a string to add to plots
 
ACE_Mn_lm_p <- function(ACE_Mn_lm) {
  f <- summary(ACE_Mn_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_lm_p(ACE_Mn_lm)

ACE_Mn_lm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_lm_p(ACE_Mn_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Mn_wlm_p <- function(ACE_Mn_wlm) {
  f <- summary(ACE_Mn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wlm_p(ACE_Mn_wlm)

ACE_Mn_wlm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wlm_p(ACE_Mn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Mn_wlm_p <- function(ACE_Mn_wlm) {
  f <- summary(ACE_Mn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wlm_p(ACE_Mn_wlm)

ACE_Mn_wlm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/(Mn_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wlm_p(ACE_Mn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Mn_wls_p <- function(ACE_Mn_wls) {
  f <- summary(ACE_Mn_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wls_p(ACE_Mn_wls)

ACE_Mn_wls_eqn <- function(df){
  m <- ACE_Mn_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wls_p(ACE_Mn_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Mn_wls_wt_p <- function(ACE_Mn_wls_wt) {
  f <- summary(ACE_Mn_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wls_wt_p(ACE_Mn_wls_wt)

ACE_Mn_wls_wt_eqn <- function(df){
  m <- ACE_Mn_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wls_wt_p(ACE_Mn_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Mn_final <- ACE_Mn + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Mn_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Mn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Mn_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Mn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Mn_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Mn_x.reg <- ACE_dataset$Mn
Mn_y.reg <- ACE_dataset$Mn_ICP
# Build linear regression model
Mn_model1 <- lm(Mn_y.reg ~ Mn, data = ACE_dataset) # OSL unweighted
Mn_model2 <- lm(Mn_y.reg ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt) # WLS unweighted
# Model summary statistics
Mn_model1
summary(Mn_model1)
confint(Mn_model1)
Mn_model2
summary(Mn_model2)
confint(Mn_model2)
# Create predicted values and upper/lower CI to check model is working
Mn_new1.y <- data.frame(Mn_x.reg = c(1, 2, 3))
predict(Mn_model1, Mn_newdata1 = Mn_new1.y)
predict(Mn_model1, Mn_newdata1 = Mn_new1.y, Mn_interval1 = "confidence")
Mn_new2.y <- data.frame(Mn_x.reg = c(1, 2, 3))
predict(Mn_model2, Mn_newdata2 = Mn_new2.y)
predict(Mn_model2, Mn_newdata2 = Mn_new2.y, Mn_interval2 = "confidence")
# Add prediction intervals to model data frame
Mn_pred.int1 <- predict(Mn_model1, interval = "prediction")
Mn_data_1_out <- bind_cols(ACE_dataset, Mn_pred.int1) %>%
  rename(Mn_fit_OLS = fit, Mn_lwr_OLS = lwr, Mn_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Mn, Mn_sd, Mn_ICP, Mn_ICP_sd, Mn_fit_OLS, Mn_lwr_OLS, Mn_upr_OLS))
Mn_pred.int2 <- predict(Mn_model2, interval = "prediction")
Mn_data_2_out <- bind_cols(Mn_data_1_out, Mn_pred.int2) %>%
  rename(Mn_fit_WLS = fit, Mn_lwr_WLS = lwr, Mn_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Mn, Mn_sd, Mn_ICP, Mn_ICP_sd, Mn_fit_OLS, Mn_lwr_OLS, Mn_upr_OLS, Mn_fit_WLS, Mn_lwr_WLS, Mn_upr_WLS)) %>% 
  print()
write.csv(Mn_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Mn_predict <- ACE_Mn_final + 
  geom_line(data = Mn_data_1_out, aes(y = Mn_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Mn_data_1_out, aes(y = Mn_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Mn_data_2_out, aes(y = Mn_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Mn_data_2_out, aes(y = Mn_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Mn_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Mn/Mn_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")


# ACE_Fe LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Fe_lm <- lm(Fe_ICP ~ Fe, data = ACE_dataset)
summary(ACE_Fe_lm)
glance(ACE_Fe_lm)
model_performance(ACE_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Fe_wlm <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/(Fe_ICP_sd)^2)
summary(ACE_Fe_wlm)
glance(ACE_Fe_wlm)
model_performance(ACE_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Fe_model <- lm(Fe_ICP ~ Fe, data = ACE_dataset) # define model
ACE_Fe_wt <- 1 / lm(abs(ACE_Fe_model$residuals) ~ ACE_Fe_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Fe_wls <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt) #perform weighted least squares regression
# Checks
summary(ACE_Fe_wls) # summary stats
glance(ACE_Fe_wls) # summary stats including AIC
model_performance(ACE_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Fe_model_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/Fe_ICP_sd^2) # define model
ACE_Fe_wt_wt <- 1 / lm(abs(ACE_Fe_model_wt$residuals) ~ ACE_Fe_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Fe_wls_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Fe_wls_wt) # summary stats
glance(ACE_Fe_wls_wt) # summary stats including AIC
model_performance(ACE_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_lm_hats <- as.data.frame(hatvalues(ACE_Fe_lm))
ACE_Fe_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_lm_cooksD <- cooks.distance(ACE_Fe_lm)
ACE_Fe_lm_influential <- ACE_Fe_lm_cooksD[(ACE_Fe_lm_cooksD > (3 * mean(ACE_Fe_lm_cooksD, na.rm = TRUE)))]
ACE_Fe_lm_influential
ACE_Fe_lm_influential_names <- names(ACE_Fe_lm_influential)
ACE_Fe_lm_outliers <- ACE_dataset[ACE_Fe_lm_influential_names,] # outliers only using of index values
ACE_Fe_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Fe_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_wlm_hats <- as.data.frame(hatvalues(ACE_Fe_wlm))
ACE_Fe_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_wlm_cooksD <- cooks.distance(ACE_Fe_wlm)
ACE_Fe_wlm_influential <- ACE_Fe_wlm_cooksD[(ACE_Fe_wlm_cooksD > (3 * mean(ACE_Fe_wlm_cooksD, na.rm = TRUE)))]
ACE_Fe_wlm_influential
ACE_Fe_wlm_influential_names <- names(ACE_Fe_wlm_influential)
ACE_Fe_wlm_outliers <- ACE_dataset[ACE_Fe_wlm_influential_names,] # outliers only using of index values
ACE_Fe_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Fe_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_wls_hats <- as.data.frame(hatvalues(ACE_Fe_wls))
ACE_Fe_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_wls_cooksD <- cooks.distance(ACE_Fe_wls)
ACE_Fe_wls_influential <- ACE_Fe_wls_cooksD[(ACE_Fe_wls_cooksD > (3 * mean(ACE_Fe_wls_cooksD, na.rm = TRUE)))]
ACE_Fe_wls_influential
ACE_Fe_wls_influential_names <- names(ACE_Fe_wls_influential)
ACE_Fe_wls_outliers <- ACE_dataset[ACE_Fe_wls_influential_names,] # outliers only using of index values
ACE_Fe_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Fe_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_wls_wt_hats <- as.data.frame(hatvalues(ACE_Fe_wls_wt))
ACE_Fe_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_wls_wt_cooksD <- cooks.distance(ACE_Fe_wls_wt)
ACE_Fe_wls_wt_influential <- ACE_Fe_wls_wt_cooksD[(ACE_Fe_wls_wt_cooksD > (3 * mean(ACE_Fe_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Fe_wls_wt_influential
ACE_Fe_wls_wt_influential_names <- names(ACE_Fe_wls_wt_influential)
ACE_Fe_wls_wt_outliers <- ACE_dataset[ACE_Fe_wls_wt_influential_names,] # outliers only using of index values
ACE_Fe_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Fe_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_OLS_summary.txt")
summary(ACE_Fe_lm)
glance(ACE_Fe_lm)
model_performance(ACE_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_lm_lev.pdf")
leveragePlots(ACE_Fe_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_OLS_wt_summary.txt")
summary(ACE_Fe_wlm)
glance(ACE_Fe_wlm)
model_performance(ACE_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wlm_lev.pdf")
leveragePlots(ACE_Fe_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_WLS_summary.txt")
summary(ACE_Fe_wls)
glance(ACE_Fe_wls)
model_performance(ACE_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wls_lev.pdf")
leveragePlots(ACE_Fe_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_WLS_wt_summary.txt")
summary(ACE_Fe_wls_wt)
glance(ACE_Fe_wls_wt)
model_performance(ACE_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wls_wt_lev.pdf")
leveragePlots(ACE_Fe_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/ACE_Fe_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Fe"
Fe_WLS_wt = ACE_Fe_wt
Fe_WLS_err_wt = ACE_Fe_wt_wt
theme_set(theme_classic(10))

ACE_Fe <- ggplot(ACE_dataset, aes(x = Fe, y = Fe_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Fe_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Fe_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Fe_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Fe

# Define p value, equation & R2 as a string to add to plots
 
ACE_Fe_lm_p <- function(ACE_Fe_lm) {
  f <- summary(ACE_Fe_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_lm_p(ACE_Fe_lm)

ACE_Fe_lm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_lm_p(ACE_Fe_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Fe_wlm_p <- function(ACE_Fe_wlm) {
  f <- summary(ACE_Fe_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wlm_p(ACE_Fe_wlm)

ACE_Fe_wlm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wlm_p(ACE_Fe_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Fe_wlm_p <- function(ACE_Fe_wlm) {
  f <- summary(ACE_Fe_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wlm_p(ACE_Fe_wlm)

ACE_Fe_wlm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/(Fe_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wlm_p(ACE_Fe_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Fe_wls_p <- function(ACE_Fe_wls) {
  f <- summary(ACE_Fe_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wls_p(ACE_Fe_wls)

ACE_Fe_wls_eqn <- function(df){
  m <- ACE_Fe_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wls_p(ACE_Fe_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Fe_wls_wt_p <- function(ACE_Fe_wls_wt) {
  f <- summary(ACE_Fe_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wls_wt_p(ACE_Fe_wls_wt)

ACE_Fe_wls_wt_eqn <- function(df){
  m <- ACE_Fe_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wls_wt_p(ACE_Fe_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Fe_final <- ACE_Fe + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Fe_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Fe_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Fe_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Fe_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Fe_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Fe_x.reg <- ACE_dataset$Fe
Fe_y.reg <- ACE_dataset$Fe_ICP
# Build linear regression model
Fe_model1 <- lm(Fe_y.reg ~ Fe, data = ACE_dataset) # OSL unweighted
Fe_model2 <- lm(Fe_y.reg ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt) # WLS unweighted
# Model summary statistics
Fe_model1
summary(Fe_model1)
confint(Fe_model1)
Fe_model2
summary(Fe_model2)
confint(Fe_model2)
# Create predicted values and upper/lower CI to check model is working
Fe_new1.y <- data.frame(Fe_x.reg = c(1, 2, 3))
predict(Fe_model1, Fe_newdata1 = Fe_new1.y)
predict(Fe_model1, Fe_newdata1 = Fe_new1.y, Fe_interval1 = "confidence")
Fe_new2.y <- data.frame(Fe_x.reg = c(1, 2, 3))
predict(Fe_model2, Fe_newdata2 = Fe_new2.y)
predict(Fe_model2, Fe_newdata2 = Fe_new2.y, Fe_interval2 = "confidence")
# Add prediction intervals to model data frame
Fe_pred.int1 <- predict(Fe_model1, interval = "prediction")
Fe_data_1_out <- bind_cols(ACE_dataset, Fe_pred.int1) %>%
  rename(Fe_fit_OLS = fit, Fe_lwr_OLS = lwr, Fe_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Fe, Fe_sd, Fe_ICP, Fe_ICP_sd, Fe_fit_OLS, Fe_lwr_OLS, Fe_upr_OLS))
Fe_pred.int2 <- predict(Fe_model2, interval = "prediction")
Fe_data_2_out <- bind_cols(Fe_data_1_out, Fe_pred.int2) %>%
  rename(Fe_fit_WLS = fit, Fe_lwr_WLS = lwr, Fe_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Fe, Fe_sd, Fe_ICP, Fe_ICP_sd, Fe_fit_OLS, Fe_lwr_OLS, Fe_upr_OLS, Fe_fit_WLS, Fe_lwr_WLS, Fe_upr_WLS)) %>% 
  print()
write.csv(Fe_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Fe_predict <- ACE_Fe_final + 
  geom_line(data = Fe_data_1_out, aes(y = Fe_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Fe_data_1_out, aes(y = Fe_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Fe_data_2_out, aes(y = Fe_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Fe_data_2_out, aes(y = Fe_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Fe_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Fe/Fe_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Sr LM -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Sr_lm <- lm(Sr_ICP ~ Sr, data = ACE_dataset)
summary(ACE_Sr_lm)
glance(ACE_Sr_lm)
model_performance(ACE_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Sr_wlm <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/(Sr_ICP_sd)^2)
summary(ACE_Sr_wlm)
glance(ACE_Sr_wlm)
model_performance(ACE_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Sr_model <- lm(Sr_ICP ~ Sr, data = ACE_dataset) # define model
ACE_Sr_wt <- 1 / lm(abs(ACE_Sr_model$residuals) ~ ACE_Sr_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Sr_wls <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt) #perform weighted least squares regression
# Checks
summary(ACE_Sr_wls) # summary stats
glance(ACE_Sr_wls) # summary stats including AIC
model_performance(ACE_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Sr_model_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/Sr_ICP_sd^2) # define model
ACE_Sr_wt_wt <- 1 / lm(abs(ACE_Sr_model_wt$residuals) ~ ACE_Sr_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Sr_wls_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Sr_wls_wt) # summary stats
glance(ACE_Sr_wls_wt) # summary stats including AIC
model_performance(ACE_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_lm_hats <- as.data.frame(hatvalues(ACE_Sr_lm))
ACE_Sr_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_lm_cooksD <- cooks.distance(ACE_Sr_lm)
ACE_Sr_lm_influential <- ACE_Sr_lm_cooksD[(ACE_Sr_lm_cooksD > (3 * mean(ACE_Sr_lm_cooksD, na.rm = TRUE)))]
ACE_Sr_lm_influential
ACE_Sr_lm_influential_names <- names(ACE_Sr_lm_influential)
ACE_Sr_lm_outliers <- ACE_dataset[ACE_Sr_lm_influential_names,] # outliers only using of index values
ACE_Sr_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Sr_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_wlm_hats <- as.data.frame(hatvalues(ACE_Sr_wlm))
ACE_Sr_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_wlm_cooksD <- cooks.distance(ACE_Sr_wlm)
ACE_Sr_wlm_influential <- ACE_Sr_wlm_cooksD[(ACE_Sr_wlm_cooksD > (3 * mean(ACE_Sr_wlm_cooksD, na.rm = TRUE)))]
ACE_Sr_wlm_influential
ACE_Sr_wlm_influential_names <- names(ACE_Sr_wlm_influential)
ACE_Sr_wlm_outliers <- ACE_dataset[ACE_Sr_wlm_influential_names,] # outliers only using of index values
ACE_Sr_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Sr_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_wls_hats <- as.data.frame(hatvalues(ACE_Sr_wls))
ACE_Sr_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_wls_cooksD <- cooks.distance(ACE_Sr_wls)
ACE_Sr_wls_influential <- ACE_Sr_wls_cooksD[(ACE_Sr_wls_cooksD > (3 * mean(ACE_Sr_wls_cooksD, na.rm = TRUE)))]
ACE_Sr_wls_influential
ACE_Sr_wls_influential_names <- names(ACE_Sr_wls_influential)
ACE_Sr_wls_outliers <- ACE_dataset[ACE_Sr_wls_influential_names,] # outliers only using of index values
ACE_Sr_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Sr_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_wls_wt_hats <- as.data.frame(hatvalues(ACE_Sr_wls_wt))
ACE_Sr_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_wls_wt_cooksD <- cooks.distance(ACE_Sr_wls_wt)
ACE_Sr_wls_wt_influential <- ACE_Sr_wls_wt_cooksD[(ACE_Sr_wls_wt_cooksD > (3 * mean(ACE_Sr_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Sr_wls_wt_influential
ACE_Sr_wls_wt_influential_names <- names(ACE_Sr_wls_wt_influential)
ACE_Sr_wls_wt_outliers <- ACE_dataset[ACE_Sr_wls_wt_influential_names,] # outliers only using of index values
ACE_Sr_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Sr_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_OLS_summary.txt")
summary(ACE_Sr_lm)
glance(ACE_Sr_lm)
model_performance(ACE_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_lm_lev.pdf")
leveragePlots(ACE_Sr_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_OLS_wt_summary.txt")
summary(ACE_Sr_wlm)
glance(ACE_Sr_wlm)
model_performance(ACE_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wlm_lev.pdf")
leveragePlots(ACE_Sr_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_WLS_summary.txt")
summary(ACE_Sr_wls)
glance(ACE_Sr_wls)
model_performance(ACE_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wls_lev.pdf")
leveragePlots(ACE_Sr_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_WLS_wt_summary.txt")
summary(ACE_Sr_wls_wt)
glance(ACE_Sr_wls_wt)
model_performance(ACE_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wls_wt_lev.pdf")
leveragePlots(ACE_Sr_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/ACE_Sr_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Sr"
Sr_WLS_wt = ACE_Sr_wt
Sr_WLS_err_wt = ACE_Sr_wt_wt
theme_set(theme_classic(10))

ACE_Sr <- ggplot(ACE_dataset, aes(x = Sr, y = Sr_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Sr_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Sr_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Sr_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Sr

# Define p value, equation & R2 as a string to add to plots
 
ACE_Sr_lm_p <- function(ACE_Sr_lm) {
  f <- summary(ACE_Sr_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_lm_p(ACE_Sr_lm)

ACE_Sr_lm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_lm_p(ACE_Sr_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Sr_wlm_p <- function(ACE_Sr_wlm) {
  f <- summary(ACE_Sr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wlm_p(ACE_Sr_wlm)

ACE_Sr_wlm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wlm_p(ACE_Sr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Sr_wlm_p <- function(ACE_Sr_wlm) {
  f <- summary(ACE_Sr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wlm_p(ACE_Sr_wlm)

ACE_Sr_wlm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/(Sr_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wlm_p(ACE_Sr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Sr_wls_p <- function(ACE_Sr_wls) {
  f <- summary(ACE_Sr_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wls_p(ACE_Sr_wls)

ACE_Sr_wls_eqn <- function(df){
  m <- ACE_Sr_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wls_p(ACE_Sr_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Sr_wls_wt_p <- function(ACE_Sr_wls_wt) {
  f <- summary(ACE_Sr_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wls_wt_p(ACE_Sr_wls_wt)

ACE_Sr_wls_wt_eqn <- function(df){
  m <- ACE_Sr_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wls_wt_p(ACE_Sr_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Sr_final <- ACE_Sr + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Sr_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Sr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Sr_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Sr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Sr_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Sr_x.reg <- ACE_dataset$Sr
Sr_y.reg <- ACE_dataset$Sr_ICP
# Build linear regression model
Sr_model1 <- lm(Sr_y.reg ~ Sr, data = ACE_dataset) # OSL unweighted
Sr_model2 <- lm(Sr_y.reg ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt) # WLS unweighted
# Model summary statistics
Sr_model1
summary(Sr_model1)
confint(Sr_model1)
Sr_model2
summary(Sr_model2)
confint(Sr_model2)
# Create predicted values and upper/lower CI to check model is working
Sr_new1.y <- data.frame(Sr_x.reg = c(1, 2, 3))
predict(Sr_model1, Sr_newdata1 = Sr_new1.y)
predict(Sr_model1, Sr_newdata1 = Sr_new1.y, Sr_interval1 = "confidence")
Sr_new2.y <- data.frame(Sr_x.reg = c(1, 2, 3))
predict(Sr_model2, Sr_newdata2 = Sr_new2.y)
predict(Sr_model2, Sr_newdata2 = Sr_new2.y, Sr_interval2 = "confidence")
# Add prediction intervals to model data frame
Sr_pred.int1 <- predict(Sr_model1, interval = "prediction")
Sr_data_1_out <- bind_cols(ACE_dataset, Sr_pred.int1) %>%
  rename(Sr_fit_OLS = fit, Sr_lwr_OLS = lwr, Sr_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Sr, Sr_sd, Sr_ICP, Sr_ICP_sd, Sr_fit_OLS, Sr_lwr_OLS, Sr_upr_OLS))
Sr_pred.int2 <- predict(Sr_model2, interval = "prediction")
Sr_data_2_out <- bind_cols(Sr_data_1_out, Sr_pred.int2) %>%
  rename(Sr_fit_WLS = fit, Sr_lwr_WLS = lwr, Sr_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Sr, Sr_sd, Sr_ICP, Sr_ICP_sd, Sr_fit_OLS, Sr_lwr_OLS, Sr_upr_OLS, Sr_fit_WLS, Sr_lwr_WLS, Sr_upr_WLS)) %>% 
  print()
write.csv(Sr_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Sr_predict <- ACE_Sr_final + 
  geom_line(data = Sr_data_1_out, aes(y = Sr_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Sr_data_1_out, aes(y = Sr_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Sr_data_2_out, aes(y = Sr_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Sr_data_2_out, aes(y = Sr_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Sr_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Sr/Sr_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")


# ACE_Zr LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Zr_lm <- lm(Zr_ICP ~ Zr, data = ACE_dataset)
summary(ACE_Zr_lm)
glance(ACE_Zr_lm)
model_performance(ACE_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Zr_wlm <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/(Zr_ICP_sd)^2)
summary(ACE_Zr_wlm)
glance(ACE_Zr_wlm)
model_performance(ACE_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Zr_model <- lm(Zr_ICP ~ Zr, data = ACE_dataset) # define model
ACE_Zr_wt <- 1 / lm(abs(ACE_Zr_model$residuals) ~ ACE_Zr_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Zr_wls <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zr_wls) # summary stats
glance(ACE_Zr_wls) # summary stats including AIC
model_performance(ACE_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Zr_model_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/Zr_ICP_sd^2) # define model
ACE_Zr_wt_wt <- 1 / lm(abs(ACE_Zr_model_wt$residuals) ~ ACE_Zr_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Zr_wls_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zr_wls_wt) # summary stats
glance(ACE_Zr_wls_wt) # summary stats including AIC
model_performance(ACE_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_lm_hats <- as.data.frame(hatvalues(ACE_Zr_lm))
ACE_Zr_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_lm_cooksD <- cooks.distance(ACE_Zr_lm)
ACE_Zr_lm_influential <- ACE_Zr_lm_cooksD[(ACE_Zr_lm_cooksD > (3 * mean(ACE_Zr_lm_cooksD, na.rm = TRUE)))]
ACE_Zr_lm_influential
ACE_Zr_lm_influential_names <- names(ACE_Zr_lm_influential)
ACE_Zr_lm_outliers <- ACE_dataset[ACE_Zr_lm_influential_names,] # outliers only using of index values
ACE_Zr_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zr_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_wlm_hats <- as.data.frame(hatvalues(ACE_Zr_wlm))
ACE_Zr_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_wlm_cooksD <- cooks.distance(ACE_Zr_wlm)
ACE_Zr_wlm_influential <- ACE_Zr_wlm_cooksD[(ACE_Zr_wlm_cooksD > (3 * mean(ACE_Zr_wlm_cooksD, na.rm = TRUE)))]
ACE_Zr_wlm_influential
ACE_Zr_wlm_influential_names <- names(ACE_Zr_wlm_influential)
ACE_Zr_wlm_outliers <- ACE_dataset[ACE_Zr_wlm_influential_names,] # outliers only using of index values
ACE_Zr_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zr_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_wls_hats <- as.data.frame(hatvalues(ACE_Zr_wls))
ACE_Zr_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_wls_cooksD <- cooks.distance(ACE_Zr_wls)
ACE_Zr_wls_influential <- ACE_Zr_wls_cooksD[(ACE_Zr_wls_cooksD > (3 * mean(ACE_Zr_wls_cooksD, na.rm = TRUE)))]
ACE_Zr_wls_influential
ACE_Zr_wls_influential_names <- names(ACE_Zr_wls_influential)
ACE_Zr_wls_outliers <- ACE_dataset[ACE_Zr_wls_influential_names,] # outliers only using of index values
ACE_Zr_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zr_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_wls_wt_hats <- as.data.frame(hatvalues(ACE_Zr_wls_wt))
ACE_Zr_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_wls_wt_cooksD <- cooks.distance(ACE_Zr_wls_wt)
ACE_Zr_wls_wt_influential <- ACE_Zr_wls_wt_cooksD[(ACE_Zr_wls_wt_cooksD > (3 * mean(ACE_Zr_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Zr_wls_wt_influential
ACE_Zr_wls_wt_influential_names <- names(ACE_Zr_wls_wt_influential)
ACE_Zr_wls_wt_outliers <- ACE_dataset[ACE_Zr_wls_wt_influential_names,] # outliers only using of index values
ACE_Zr_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zr_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_OLS_summary.txt")
summary(ACE_Zr_lm)
glance(ACE_Zr_lm)
model_performance(ACE_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_lm_lev.pdf")
leveragePlots(ACE_Zr_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_OLS_wt_summary.txt")
summary(ACE_Zr_wlm)
glance(ACE_Zr_wlm)
model_performance(ACE_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wlm_lev.pdf")
leveragePlots(ACE_Zr_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_WLS_summary.txt")
summary(ACE_Zr_wls)
glance(ACE_Zr_wls)
model_performance(ACE_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wls_lev.pdf")
leveragePlots(ACE_Zr_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_WLS_wt_summary.txt")
summary(ACE_Zr_wls_wt)
glance(ACE_Zr_wls_wt)
model_performance(ACE_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wls_wt_lev.pdf")
leveragePlots(ACE_Zr_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/ACE_Zr_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Zr"
Zr_WLS_wt = ACE_Zr_wt
Zr_WLS_err_wt = ACE_Zr_wt_wt
theme_set(theme_classic(10))

ACE_Zr <- ggplot(ACE_dataset, aes(x = Zr, y = Zr_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Zr_ICP-Zr_ICP_sd, ymax=Zr_ICP+Zr_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Zr-Zr_sd, xmax=Zr+Zr_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Zr_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Zr_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Zr_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Zr

# Define p value, equation & R2 as a string to add to plots
 
ACE_Zr_lm_p <- function(ACE_Zr_lm) {
  f <- summary(ACE_Zr_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_lm_p(ACE_Zr_lm)

ACE_Zr_lm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_lm_p(ACE_Zr_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Zr_wlm_p <- function(ACE_Zr_wlm) {
  f <- summary(ACE_Zr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wlm_p(ACE_Zr_wlm)

ACE_Zr_wlm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wlm_p(ACE_Zr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Zr_wlm_p <- function(ACE_Zr_wlm) {
  f <- summary(ACE_Zr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wlm_p(ACE_Zr_wlm)

ACE_Zr_wlm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/(Zr_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wlm_p(ACE_Zr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Zr_wls_p <- function(ACE_Zr_wls) {
  f <- summary(ACE_Zr_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wls_p(ACE_Zr_wls)

ACE_Zr_wls_eqn <- function(df){
  m <- ACE_Zr_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wls_p(ACE_Zr_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Zr_wls_wt_p <- function(ACE_Zr_wls_wt) {
  f <- summary(ACE_Zr_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wls_wt_p(ACE_Zr_wls_wt)

ACE_Zr_wls_wt_eqn <- function(df){
  m <- ACE_Zr_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wls_wt_p(ACE_Zr_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Zr_final <- ACE_Zr + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Zr_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Zr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Zr_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Zr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Zr_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Zr_x.reg <- ACE_dataset$Zr
Zr_y.reg <- ACE_dataset$Zr_ICP
# Build linear regression model
Zr_model1 <- lm(Zr_y.reg ~ Zr, data = ACE_dataset) # OSL unweighted
Zr_model2 <- lm(Zr_y.reg ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt) # WLS unweighted
# Model summary statistics
Zr_model1
summary(Zr_model1)
confint(Zr_model1)
Zr_model2
summary(Zr_model2)
confint(Zr_model2)
# Create predicted values and upper/lower CI to check model is working
Zr_new1.y <- data.frame(Zr_x.reg = c(1, 2, 3))
predict(Zr_model1, Zr_newdata1 = Zr_new1.y)
predict(Zr_model1, Zr_newdata1 = Zr_new1.y, Zr_interval1 = "confidence")
Zr_new2.y <- data.frame(Zr_x.reg = c(1, 2, 3))
predict(Zr_model2, Zr_newdata2 = Zr_new2.y)
predict(Zr_model2, Zr_newdata2 = Zr_new2.y, Zr_interval2 = "confidence")
# Add prediction intervals to model data frame
Zr_pred.int1 <- predict(Zr_model1, interval = "prediction")
Zr_data_1_out <- bind_cols(ACE_dataset, Zr_pred.int1) %>%
  rename(Zr_fit_OLS = fit, Zr_lwr_OLS = lwr, Zr_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Zr, Zr_sd, Zr_ICP, Zr_ICP_sd, Zr_fit_OLS, Zr_lwr_OLS, Zr_upr_OLS))
Zr_pred.int2 <- predict(Zr_model2, interval = "prediction")
Zr_data_2_out <- bind_cols(Zr_data_1_out, Zr_pred.int2) %>%
  rename(Zr_fit_WLS = fit, Zr_lwr_WLS = lwr, Zr_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Zr, Zr_sd, Zr_ICP, Zr_ICP_sd, Zr_fit_OLS, Zr_lwr_OLS, Zr_upr_OLS, Zr_fit_WLS, Zr_lwr_WLS, Zr_upr_WLS)) %>% 
  print()
write.csv(Zr_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Zr_predict <- ACE_Zr_final + 
  geom_line(data = Zr_data_1_out, aes(y = Zr_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Zr_data_1_out, aes(y = Zr_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Zr_data_2_out, aes(y = Zr_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Zr_data_2_out, aes(y = Zr_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Zr_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zr/Zr_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 3a - Primary matched elements -------------------------------------------------------------------------

# Summary 3x2 matrix plots of ITRAX-acf primary matched elements
# OLS & WLS summary - unweighted stats on plot
ggarrange(ACE_Ca_predict, ACE_Ti_predict, ACE_Mn_predict, 
          ACE_Fe_predict, ACE_Sr_predict, ACE_Zr_predict,
          ncol = 2, nrow = 3, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/ACE_Fig3b_OLS_WLS_clr_key_predict.pdf",
       height = c(40), width = c(30), dpi = 600, units = "cm")


# Secondary Matched Elements shown in Supplementary Material ----------------------------------------------
# ACE_K LM -----------------------------------------------------------------------

# 1) OLS (Ordinary Least Squares) - linear model
ACE_K_lm <- lm(K_ICP ~ K, data = ACE_dataset)
summary(ACE_K_lm)
glance(ACE_K_lm)
model_performance(ACE_K_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_lm_hats <- as.data.frame(hatvalues(ACE_K_lm))
# Cooks distance - 2-3 x difference from mean 
ACE_K_lm_cooksD <- cooks.distance(ACE_K_lm)
ACE_K_lm_influential <- ACE_K_lm_cooksD[(ACE_K_lm_cooksD > (3 * mean(ACE_K_lm_cooksD, na.rm = TRUE)))]
ACE_K_lm_influential
ACE_K_lm_influential_names <- names(ACE_K_lm_influential)
ACE_K_lm_outliers <- ACE_dataset[ACE_K_lm_influential_names,] # outliers only using of index values
ACE_K_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_K_lm_outliers) # generates a generates a new dataset with outliers removed
write.csv(ACE_K_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model & checks
ACE_K_wlm <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/(K_ICP_sd)^2)
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_wlm_hats <- as.data.frame(hatvalues(ACE_K_wlm))
# Cooks distance - 2-3 x difference from mean 
ACE_K_wlm_cooksD <- cooks.distance(ACE_K_wlm)
ACE_K_wlm_influential <- ACE_K_wlm_cooksD[(ACE_K_wlm_cooksD > (3 * mean(ACE_K_wlm_cooksD, na.rm = TRUE)))]
ACE_K_wlm_influential
ACE_K_wlm_influential_names <- names(ACE_K_wlm_influential)
ACE_K_wlm_outliers <- ACE_dataset[ACE_K_wlm_influential_names,] # outliers only using of index values
ACE_K_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_K_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_K_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) WLS (Weighted Least Squares)
ACE_K_model <- lm(K_ICP ~ K, data = ACE_dataset) # define model
ACE_K_wt <- 1 / lm(abs(ACE_K_model$residuals) ~ ACE_K_model$fitted.values)$fitted.values^2 #define weights to use
ACE_K_wls <- lm(K_ICP ~ K, data = ACE_dataset, weights=ACE_K_wt) #perform weighted least squares regression
# Checks
summary(ACE_K_wls) # summary stats
glance(ACE_K_wls) # summary stats including AIC
model_performance(ACE_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_wls_hats <- as.data.frame(hatvalues(ACE_K_wls))
# Cooks distance - 2-3 x difference from mean 
ACE_K_wls_cooksD <- cooks.distance(ACE_K_wls)
ACE_K_wls_influential <- ACE_K_wls_cooksD[(ACE_K_wls_cooksD > (3 * mean(ACE_K_wls_cooksD, na.rm = TRUE)))]
ACE_K_wls_influential
ACE_K_wls_influential_names <- names(ACE_K_wls_influential)
ACE_K_wls_outliers <- ACE_dataset[ACE_K_wls_influential_names,] # outliers only using of index values
ACE_K_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_K_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_K_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_WLS_no_outliers.csv", row.names = FALSE)

# 4) WLS (ICPMS error weighted)
ACE_K_model_wt <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/K_ICP_sd^2) # define model
ACE_K_wt_wt <- 1 / lm(abs(ACE_K_model_wt$residuals) ~ ACE_K_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_K_wls_wt <- lm(K_ICP ~ K, data = ACE_dataset, weights=ACE_K_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_K_wls_wt) # summary stats
glance(ACE_K_wls_wt) # summary stats including AIC
model_performance(ACE_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_wls_wt_hats <- as.data.frame(hatvalues(ACE_K_wls_wt))
# Cooks distance - 2-3 x difference from mean 
ACE_K_wls_wt_cooksD <- cooks.distance(ACE_K_wls_wt)
ACE_K_wls_wt_influential <- ACE_K_wls_wt_cooksD[(ACE_K_wls_wt_cooksD > (3 * mean(ACE_K_wls_wt_cooksD, na.rm = TRUE)))]
ACE_K_wls_wt_influential
ACE_K_wls_wt_influential_names <- names(ACE_K_wls_wt_influential)
ACE_K_wls_wt_outliers <- ACE_dataset[ACE_K_wls_wt_influential_names,] # outliers only using of index values
ACE_K_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_K_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_K_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_WLS_wt_no_outliers.csv", row.names = FALSE)

# Save stats & stats tests plots to file

# 1) OLS
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_OLS_wt_summary.txt")
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wlm) # Performance package summary check for heteroscedasticity
icc(ACE_K_wlm) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_lm_lev_bar.pdf")
barplot(hatvalues(ACE_K_wlm), col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_lm_lev.pdf")
leveragePlots(ACE_K_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wlm)
dev.off()

# 2) OLS weighted
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_OLS_summary.txt")
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wlm) # Performance package summary check for heteroscedasticity
icc(ACE_K_wlm) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_K_wlm), col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wlm_lev.pdf")
leveragePlots(ACE_K_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wlm)
dev.off()

# 3) WLS
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_WLS_summary.txt")
summary(ACE_K_wls)
glance(ACE_K_wls)
model_performance(ACE_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wls) # Performance package summary check for heteroscedasticity
icc(ACE_K_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wls_lev_bar.pdf")
barplot(hatvalues(ACE_K_wls), col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wls_lev.pdf")
leveragePlots(ACE_K_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wls)
dev.off()

# 4) WLS weighted
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_WLS_wt_summary.txt")
summary(ACE_K_wls_wt)
glance(ACE_K_wls_wt)
model_performance(ACE_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_K_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_K_wls_wt), col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wls_wt_lev.pdf")
leveragePlots(ACE_K_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/ACE_K_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wls_wt)
dev.off()

# Linear model summary elemental plots

# K
element_title <- "K"
K_WLS_wt = ACE_K_wt
K_WLS_err_wt = ACE_K_wt_wt
theme_set(theme_classic(10))

ACE_K <- ggplot(ACE_dataset, aes(x = K, y = K_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/K_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = K_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = K_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_K

# Define p value, equation & R2 as a string to add to plots

ACE_K_lm_p <- function(ACE_K_lm) {
  f <- summary(ACE_K_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_lm_p(ACE_K_lm)

ACE_K_lm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_lm_p(ACE_K_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_K_wlm_p <- function(ACE_K_wlm) {
  f <- summary(ACE_K_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wlm_p(ACE_K_wlm)

ACE_K_wlm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wlm_p(ACE_K_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_K_wlm_p <- function(ACE_K_wlm) {
  f <- summary(ACE_K_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wlm_p(ACE_K_wlm)

ACE_K_wlm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/(K_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wlm_p(ACE_K_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_K_wls_p <- function(ACE_K_wls) {
  f <- summary(ACE_K_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wls_p(ACE_K_wls)

ACE_K_wls_eqn <- function(df){
  m <- ACE_K_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wls_p(ACE_K_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_K_wls_wt_p <- function(ACE_K_wls_wt) {
  f <- summary(ACE_K_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wls_wt_p(ACE_K_wls_wt)

ACE_K_wls_wt_eqn <- function(df){
  m <- ACE_K_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wls_wt_p(ACE_K_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_K_final <- ACE_K + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_K_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_K_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_K_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_K_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_K_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
K_x.reg <- ACE_dataset$K
K_y.reg <- ACE_dataset$K_ICP
# Build linear regression model
K_model1 <- lm(K_y.reg ~ K, data = ACE_dataset) # OSL unweighted
K_model2 <- lm(K_y.reg ~ K, data = ACE_dataset, weights=ACE_K_wt) # WLS unweighted
# Model summary statistics
K_model1
summary(K_model1)
confint(K_model1)
K_model2
summary(K_model2)
confint(K_model2)
# Create predicted values and upper/lower CI to check model is working
K_new1.y <- data.frame(K_x.reg = c(1, 2, 3))
predict(K_model1, K_newdata1 = K_new1.y)
predict(K_model1, K_newdata1 = K_new1.y, K_interval1 = "confidence")
K_new2.y <- data.frame(K_x.reg = c(1, 2, 3))
predict(K_model2, K_newdata2 = K_new2.y)
predict(K_model2, K_newdata2 = K_new2.y, K_interval2 = "confidence")
# Add prediction intervals to model data frame
K_pred.int1 <- predict(K_model1, interval = "prediction")
K_data_1_out <- bind_cols(ACE_dataset, K_pred.int1) %>%
  rename(K_fit_OLS = fit, K_lwr_OLS = lwr, K_upr_OLS = upr) %>% 
  select(c(Location:midpoint, K, K_sd, K_ICP, K_ICP_sd, K_fit_OLS, K_lwr_OLS, K_upr_OLS))
K_pred.int2 <- predict(K_model2, interval = "prediction")
K_data_2_out <- bind_cols(K_data_1_out, K_pred.int2) %>%
  rename(K_fit_WLS = fit, K_lwr_WLS = lwr, K_upr_WLS = upr) %>% 
  select(c(Location:midpoint, K, K_sd, K_ICP, K_ICP_sd, K_fit_OLS, K_lwr_OLS, K_upr_OLS, K_fit_WLS, K_lwr_WLS, K_upr_WLS)) %>% 
  print()
write.csv(K_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_K_predict <- ACE_K_final + 
  geom_line(data = K_data_1_out, aes(y = K_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = K_data_1_out, aes(y = K_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = K_data_2_out, aes(y = K_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = K_data_2_out, aes(y = K_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_K_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/K/K_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")



 
# ACE_Co LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Co_lm <- lm(Co_ICP ~ Co, data = ACE_dataset)
summary(ACE_Co_lm)
glance(ACE_Co_lm)
model_performance(ACE_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Co_wlm <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/(Co_ICP_sd)^2)
summary(ACE_Co_wlm)
glance(ACE_Co_wlm)
model_performance(ACE_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Co_model <- lm(Co_ICP ~ Co, data = ACE_dataset) # define model
ACE_Co_wt <- 1 / lm(abs(ACE_Co_model$residuals) ~ ACE_Co_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Co_wls <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=ACE_Co_wt) #perform weighted least squares regression
# Checks
summary(ACE_Co_wls) # summary stats
glance(ACE_Co_wls) # summary stats including AIC
model_performance(ACE_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Co_model_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/Co_ICP_sd^2) # define model
ACE_Co_wt_wt <- 1 / lm(abs(ACE_Co_model_wt$residuals) ~ ACE_Co_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Co_wls_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=ACE_Co_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Co_wls_wt) # summary stats
glance(ACE_Co_wls_wt) # summary stats including AIC
model_performance(ACE_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_lm_hats <- as.data.frame(hatvalues(ACE_Co_lm))
ACE_Co_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_lm_cooksD <- cooks.distance(ACE_Co_lm)
ACE_Co_lm_influential <- ACE_Co_lm_cooksD[(ACE_Co_lm_cooksD > (3 * mean(ACE_Co_lm_cooksD, na.rm = TRUE)))]
ACE_Co_lm_influential
ACE_Co_lm_influential_names <- names(ACE_Co_lm_influential)
ACE_Co_lm_outliers <- ACE_dataset[ACE_Co_lm_influential_names,] # outliers only using of index values
ACE_Co_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Co_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_wlm_hats <- as.data.frame(hatvalues(ACE_Co_wlm))
ACE_Co_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_wlm_cooksD <- cooks.distance(ACE_Co_wlm)
ACE_Co_wlm_influential <- ACE_Co_wlm_cooksD[(ACE_Co_wlm_cooksD > (3 * mean(ACE_Co_wlm_cooksD, na.rm = TRUE)))]
ACE_Co_wlm_influential
ACE_Co_wlm_influential_names <- names(ACE_Co_wlm_influential)
ACE_Co_wlm_outliers <- ACE_dataset[ACE_Co_wlm_influential_names,] # outliers only using of index values
ACE_Co_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Co_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_wls_hats <- as.data.frame(hatvalues(ACE_Co_wls))
ACE_Co_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_wls_cooksD <- cooks.distance(ACE_Co_wls)
ACE_Co_wls_influential <- ACE_Co_wls_cooksD[(ACE_Co_wls_cooksD > (3 * mean(ACE_Co_wls_cooksD, na.rm = TRUE)))]
ACE_Co_wls_influential
ACE_Co_wls_influential_names <- names(ACE_Co_wls_influential)
ACE_Co_wls_outliers <- ACE_dataset[ACE_Co_wls_influential_names,] # outliers only using of index values
ACE_Co_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Co_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_wls_wt_hats <- as.data.frame(hatvalues(ACE_Co_wls_wt))
ACE_Co_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_wls_wt_cooksD <- cooks.distance(ACE_Co_wls_wt)
ACE_Co_wls_wt_influential <- ACE_Co_wls_wt_cooksD[(ACE_Co_wls_wt_cooksD > (3 * mean(ACE_Co_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Co_wls_wt_influential
ACE_Co_wls_wt_influential_names <- names(ACE_Co_wls_wt_influential)
ACE_Co_wls_wt_outliers <- ACE_dataset[ACE_Co_wls_wt_influential_names,] # outliers only using of index values
ACE_Co_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Co_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_OLS_summary.txt")
summary(ACE_Co_lm)
glance(ACE_Co_lm)
model_performance(ACE_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Co_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_lm_lev.pdf")
leveragePlots(ACE_Co_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_OLS_wt_summary.txt")
summary(ACE_Co_wlm)
glance(ACE_Co_wlm)
model_performance(ACE_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Co_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wlm_lev.pdf")
leveragePlots(ACE_Co_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_WLS_summary.txt")
summary(ACE_Co_wls)
glance(ACE_Co_wls)
model_performance(ACE_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Co_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Co_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wls_lev.pdf")
leveragePlots(ACE_Co_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_WLS_wt_summary.txt")
summary(ACE_Co_wls_wt)
glance(ACE_Co_wls_wt)
model_performance(ACE_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Co_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Co_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wls_wt_lev.pdf")
leveragePlots(ACE_Co_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/ACE_Co_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Co"
Co_WLS_wt = ACE_Co_wt
Co_WLS_err_wt = ACE_Co_wt_wt
theme_set(theme_classic(10))

ACE_Co <- ggplot(ACE_dataset, aes(x = Co, y = Co_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Co_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Co_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Co_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Co

# Define p value, equation & R2 as a string to add to plots

ACE_Co_lm_p <- function(ACE_Co_lm) {
  f <- summary(ACE_Co_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_lm_p(ACE_Co_lm)

ACE_Co_lm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_lm_p(ACE_Co_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Co_wlm_p <- function(ACE_Co_wlm) {
  f <- summary(ACE_Co_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wlm_p(ACE_Co_wlm)

ACE_Co_wlm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wlm_p(ACE_Co_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Co_wlm_p <- function(ACE_Co_wlm) {
  f <- summary(ACE_Co_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wlm_p(ACE_Co_wlm)

ACE_Co_wlm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/(Co_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wlm_p(ACE_Co_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Co_wls_p <- function(ACE_Co_wls) {
  f <- summary(ACE_Co_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wls_p(ACE_Co_wls)

ACE_Co_wls_eqn <- function(df){
  m <- ACE_Co_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wls_p(ACE_Co_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Co_wls_wt_p <- function(ACE_Co_wls_wt) {
  f <- summary(ACE_Co_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wls_wt_p(ACE_Co_wls_wt)

ACE_Co_wls_wt_eqn <- function(df){
  m <- ACE_Co_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wls_wt_p(ACE_Co_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Co_final <- ACE_Co + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Co_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Co_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Co_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Co_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Co_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Co_x.reg <- ACE_dataset$Co
Co_y.reg <- ACE_dataset$Co_ICP
# Build linear regression model
Co_model1 <- lm(Co_y.reg ~ Co, data = ACE_dataset) # OSL unweighted
Co_model2 <- lm(Co_y.reg ~ Co, data = ACE_dataset, weights=ACE_Co_wt) # WLS unweighted
# Model summary statistics
Co_model1
summary(Co_model1)
confint(Co_model1)
Co_model2
summary(Co_model2)
confint(Co_model2)
# Create predicted values and upper/lower CI to check model is working
Co_new1.y <- data.frame(Co_x.reg = c(1, 2, 3))
predict(Co_model1, Co_newdata1 = Co_new1.y)
predict(Co_model1, Co_newdata1 = Co_new1.y, Co_interval1 = "confidence")
Co_new2.y <- data.frame(Co_x.reg = c(1, 2, 3))
predict(Co_model2, Co_newdata2 = Co_new2.y)
predict(Co_model2, Co_newdata2 = Co_new2.y, Co_interval2 = "confidence")
# Add prediction intervals to model data frame
Co_pred.int1 <- predict(Co_model1, interval = "prediction")
Co_data_1_out <- bind_cols(ACE_dataset, Co_pred.int1) %>%
  rename(Co_fit_OLS = fit, Co_lwr_OLS = lwr, Co_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Co, Co_sd, Co_ICP, Co_ICP_sd, Co_fit_OLS, Co_lwr_OLS, Co_upr_OLS))
Co_pred.int2 <- predict(Co_model2, interval = "prediction")
Co_data_2_out <- bind_cols(Co_data_1_out, Co_pred.int2) %>%
  rename(Co_fit_WLS = fit, Co_lwr_WLS = lwr, Co_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Co, Co_sd, Co_ICP, Co_ICP_sd, Co_fit_OLS, Co_lwr_OLS, Co_upr_OLS, Co_fit_WLS, Co_lwr_WLS, Co_upr_WLS)) %>% 
  print()
write.csv(Co_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Co_predict <- ACE_Co_final + 
  geom_line(data = Co_data_1_out, aes(y = Co_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Co_data_1_out, aes(y = Co_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Co_data_2_out, aes(y = Co_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Co_data_2_out, aes(y = Co_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Co_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Co/Co_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Ni LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ni_lm <- lm(Ni_ICP ~ Ni, data = ACE_dataset)
summary(ACE_Ni_lm)
glance(ACE_Ni_lm)
model_performance(ACE_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Ni_wlm <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/(Ni_ICP_sd)^2)
summary(ACE_Ni_wlm)
glance(ACE_Ni_wlm)
model_performance(ACE_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Ni_model <- lm(Ni_ICP ~ Ni, data = ACE_dataset) # define model
ACE_Ni_wt <- 1 / lm(abs(ACE_Ni_model$residuals) ~ ACE_Ni_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ni_wls <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ni_wls) # summary stats
glance(ACE_Ni_wls) # summary stats including AIC
model_performance(ACE_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Ni_model_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/Ni_ICP_sd^2) # define model
ACE_Ni_wt_wt <- 1 / lm(abs(ACE_Ni_model_wt$residuals) ~ ACE_Ni_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ni_wls_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ni_wls_wt) # summary stats
glance(ACE_Ni_wls_wt) # summary stats including AIC
model_performance(ACE_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_lm_hats <- as.data.frame(hatvalues(ACE_Ni_lm))
ACE_Ni_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_lm_cooksD <- cooks.distance(ACE_Ni_lm)
ACE_Ni_lm_influential <- ACE_Ni_lm_cooksD[(ACE_Ni_lm_cooksD > (3 * mean(ACE_Ni_lm_cooksD, na.rm = TRUE)))]
ACE_Ni_lm_influential
ACE_Ni_lm_influential_names <- names(ACE_Ni_lm_influential)
ACE_Ni_lm_outliers <- ACE_dataset[ACE_Ni_lm_influential_names,] # outliers only using of index values
ACE_Ni_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ni_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_wlm_hats <- as.data.frame(hatvalues(ACE_Ni_wlm))
ACE_Ni_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_wlm_cooksD <- cooks.distance(ACE_Ni_wlm)
ACE_Ni_wlm_influential <- ACE_Ni_wlm_cooksD[(ACE_Ni_wlm_cooksD > (3 * mean(ACE_Ni_wlm_cooksD, na.rm = TRUE)))]
ACE_Ni_wlm_influential
ACE_Ni_wlm_influential_names <- names(ACE_Ni_wlm_influential)
ACE_Ni_wlm_outliers <- ACE_dataset[ACE_Ni_wlm_influential_names,] # outliers only using of index values
ACE_Ni_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ni_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_wls_hats <- as.data.frame(hatvalues(ACE_Ni_wls))
ACE_Ni_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_wls_cooksD <- cooks.distance(ACE_Ni_wls)
ACE_Ni_wls_influential <- ACE_Ni_wls_cooksD[(ACE_Ni_wls_cooksD > (3 * mean(ACE_Ni_wls_cooksD, na.rm = TRUE)))]
ACE_Ni_wls_influential
ACE_Ni_wls_influential_names <- names(ACE_Ni_wls_influential)
ACE_Ni_wls_outliers <- ACE_dataset[ACE_Ni_wls_influential_names,] # outliers only using of index values
ACE_Ni_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ni_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_wls_wt_hats <- as.data.frame(hatvalues(ACE_Ni_wls_wt))
ACE_Ni_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_wls_wt_cooksD <- cooks.distance(ACE_Ni_wls_wt)
ACE_Ni_wls_wt_influential <- ACE_Ni_wls_wt_cooksD[(ACE_Ni_wls_wt_cooksD > (3 * mean(ACE_Ni_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Ni_wls_wt_influential
ACE_Ni_wls_wt_influential_names <- names(ACE_Ni_wls_wt_influential)
ACE_Ni_wls_wt_outliers <- ACE_dataset[ACE_Ni_wls_wt_influential_names,] # outliers only using of index values
ACE_Ni_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Ni_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_OLS_summary.txt")
summary(ACE_Ni_lm)
glance(ACE_Ni_lm)
model_performance(ACE_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_lm_lev.pdf")
leveragePlots(ACE_Ni_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_OLS_wt_summary.txt")
summary(ACE_Ni_wlm)
glance(ACE_Ni_wlm)
model_performance(ACE_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wlm_lev.pdf")
leveragePlots(ACE_Ni_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_WLS_summary.txt")
summary(ACE_Ni_wls)
glance(ACE_Ni_wls)
model_performance(ACE_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wls_lev.pdf")
leveragePlots(ACE_Ni_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_WLS_wt_summary.txt")
summary(ACE_Ni_wls_wt)
glance(ACE_Ni_wls_wt)
model_performance(ACE_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wls_wt_lev.pdf")
leveragePlots(ACE_Ni_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/ACE_Ni_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Ni"
Ni_WLS_wt = ACE_Ni_wt
Ni_WLS_err_wt = ACE_Ni_wt_wt
theme_set(theme_classic(10))

ACE_Ni <- ggplot(ACE_dataset, aes(x = Ni, y = Ni_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Ni_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Ni_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Ni_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Ni

# Define p value, equation & R2 as a string to add to plots

ACE_Ni_lm_p <- function(ACE_Ni_lm) {
  f <- summary(ACE_Ni_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_lm_p(ACE_Ni_lm)

ACE_Ni_lm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_lm_p(ACE_Ni_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Ni_wlm_p <- function(ACE_Ni_wlm) {
  f <- summary(ACE_Ni_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wlm_p(ACE_Ni_wlm)

ACE_Ni_wlm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wlm_p(ACE_Ni_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Ni_wlm_p <- function(ACE_Ni_wlm) {
  f <- summary(ACE_Ni_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wlm_p(ACE_Ni_wlm)

ACE_Ni_wlm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/(Ni_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wlm_p(ACE_Ni_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Ni_wls_p <- function(ACE_Ni_wls) {
  f <- summary(ACE_Ni_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wls_p(ACE_Ni_wls)

ACE_Ni_wls_eqn <- function(df){
  m <- ACE_Ni_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wls_p(ACE_Ni_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Ni_wls_wt_p <- function(ACE_Ni_wls_wt) {
  f <- summary(ACE_Ni_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wls_wt_p(ACE_Ni_wls_wt)

ACE_Ni_wls_wt_eqn <- function(df){
  m <- ACE_Ni_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wls_wt_p(ACE_Ni_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Ni_final <- ACE_Ni + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ni_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ni_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ni_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ni_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ni_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Ni_x.reg <- ACE_dataset$Ni
Ni_y.reg <- ACE_dataset$Ni_ICP
# Build linear regression model
Ni_model1 <- lm(Ni_y.reg ~ Ni, data = ACE_dataset) # OSL unweighted
Ni_model2 <- lm(Ni_y.reg ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt) # WLS unweighted
# Model summary statistics
Ni_model1
summary(Ni_model1)
confint(Ni_model1)
Ni_model2
summary(Ni_model2)
confint(Ni_model2)
# Create predicted values and upper/lower CI to check model is working
Ni_new1.y <- data.frame(Ni_x.reg = c(1, 2, 3))
predict(Ni_model1, Ni_newdata1 = Ni_new1.y)
predict(Ni_model1, Ni_newdata1 = Ni_new1.y, Ni_interval1 = "confidence")
Ni_new2.y <- data.frame(Ni_x.reg = c(1, 2, 3))
predict(Ni_model2, Ni_newdata2 = Ni_new2.y)
predict(Ni_model2, Ni_newdata2 = Ni_new2.y, Ni_interval2 = "confidence")
# Add prediction intervals to model data frame
Ni_pred.int1 <- predict(Ni_model1, interval = "prediction")
Ni_data_1_out <- bind_cols(ACE_dataset, Ni_pred.int1) %>%
  rename(Ni_fit_OLS = fit, Ni_lwr_OLS = lwr, Ni_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Ni, Ni_sd, Ni_ICP, Ni_ICP_sd, Ni_fit_OLS, Ni_lwr_OLS, Ni_upr_OLS))
Ni_pred.int2 <- predict(Ni_model2, interval = "prediction")
Ni_data_2_out <- bind_cols(Ni_data_1_out, Ni_pred.int2) %>%
  rename(Ni_fit_WLS = fit, Ni_lwr_WLS = lwr, Ni_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Ni, Ni_sd, Ni_ICP, Ni_ICP_sd, Ni_fit_OLS, Ni_lwr_OLS, Ni_upr_OLS, Ni_fit_WLS, Ni_lwr_WLS, Ni_upr_WLS)) %>% 
  print()
write.csv(Ni_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Ni_predict <- ACE_Ni_final + 
  geom_line(data = Ni_data_1_out, aes(y = Ni_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Ni_data_1_out, aes(y = Ni_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Ni_data_2_out, aes(y = Ni_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Ni_data_2_out, aes(y = Ni_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Ni_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Ni/Ni_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Cu LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Cu_lm <- lm(Cu_ICP ~ Cu, data = ACE_dataset)
summary(ACE_Cu_lm)
glance(ACE_Cu_lm)
model_performance(ACE_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Cu_wlm <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/(Cu_ICP_sd)^2)
summary(ACE_Cu_wlm)
glance(ACE_Cu_wlm)
model_performance(ACE_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Cu_model <- lm(Cu_ICP ~ Cu, data = ACE_dataset) # define model
ACE_Cu_wt <- 1 / lm(abs(ACE_Cu_model$residuals) ~ ACE_Cu_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Cu_wls <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt) #perform weighted least squares regression
# Checks
summary(ACE_Cu_wls) # summary stats
glance(ACE_Cu_wls) # summary stats including AIC
model_performance(ACE_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Cu_model_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/Cu_ICP_sd^2) # define model
ACE_Cu_wt_wt <- 1 / lm(abs(ACE_Cu_model_wt$residuals) ~ ACE_Cu_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Cu_wls_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Cu_wls_wt) # summary stats
glance(ACE_Cu_wls_wt) # summary stats including AIC
model_performance(ACE_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_lm_hats <- as.data.frame(hatvalues(ACE_Cu_lm))
ACE_Cu_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_lm_cooksD <- cooks.distance(ACE_Cu_lm)
ACE_Cu_lm_influential <- ACE_Cu_lm_cooksD[(ACE_Cu_lm_cooksD > (3 * mean(ACE_Cu_lm_cooksD, na.rm = TRUE)))]
ACE_Cu_lm_influential
ACE_Cu_lm_influential_names <- names(ACE_Cu_lm_influential)
ACE_Cu_lm_outliers <- ACE_dataset[ACE_Cu_lm_influential_names,] # outliers only using of index values
ACE_Cu_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Cu_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_wlm_hats <- as.data.frame(hatvalues(ACE_Cu_wlm))
ACE_Cu_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_wlm_cooksD <- cooks.distance(ACE_Cu_wlm)
ACE_Cu_wlm_influential <- ACE_Cu_wlm_cooksD[(ACE_Cu_wlm_cooksD > (3 * mean(ACE_Cu_wlm_cooksD, na.rm = TRUE)))]
ACE_Cu_wlm_influential
ACE_Cu_wlm_influential_names <- names(ACE_Cu_wlm_influential)
ACE_Cu_wlm_outliers <- ACE_dataset[ACE_Cu_wlm_influential_names,] # outliers only using of index values
ACE_Cu_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Cu_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_wls_hats <- as.data.frame(hatvalues(ACE_Cu_wls))
ACE_Cu_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_wls_cooksD <- cooks.distance(ACE_Cu_wls)
ACE_Cu_wls_influential <- ACE_Cu_wls_cooksD[(ACE_Cu_wls_cooksD > (3 * mean(ACE_Cu_wls_cooksD, na.rm = TRUE)))]
ACE_Cu_wls_influential
ACE_Cu_wls_influential_names <- names(ACE_Cu_wls_influential)
ACE_Cu_wls_outliers <- ACE_dataset[ACE_Cu_wls_influential_names,] # outliers only using of index values
ACE_Cu_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Cu_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_wls_wt_hats <- as.data.frame(hatvalues(ACE_Cu_wls_wt))
ACE_Cu_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_wls_wt_cooksD <- cooks.distance(ACE_Cu_wls_wt)
ACE_Cu_wls_wt_influential <- ACE_Cu_wls_wt_cooksD[(ACE_Cu_wls_wt_cooksD > (3 * mean(ACE_Cu_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Cu_wls_wt_influential
ACE_Cu_wls_wt_influential_names <- names(ACE_Cu_wls_wt_influential)
ACE_Cu_wls_wt_outliers <- ACE_dataset[ACE_Cu_wls_wt_influential_names,] # outliers only using of index values
ACE_Cu_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Cu_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_OLS_summary.txt")
summary(ACE_Cu_lm)
glance(ACE_Cu_lm)
model_performance(ACE_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_lm_lev.pdf")
leveragePlots(ACE_Cu_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_OLS_wt_summary.txt")
summary(ACE_Cu_wlm)
glance(ACE_Cu_wlm)
model_performance(ACE_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wlm_lev.pdf")
leveragePlots(ACE_Cu_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_WLS_summary.txt")
summary(ACE_Cu_wls)
glance(ACE_Cu_wls)
model_performance(ACE_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wls_lev.pdf")
leveragePlots(ACE_Cu_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_WLS_wt_summary.txt")
summary(ACE_Cu_wls_wt)
glance(ACE_Cu_wls_wt)
model_performance(ACE_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wls_wt_lev.pdf")
leveragePlots(ACE_Cu_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/ACE_Cu_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Cu"
Cu_WLS_wt = ACE_Cu_wt
Cu_WLS_err_wt = ACE_Cu_wt_wt
theme_set(theme_classic(10))

ACE_Cu <- ggplot(ACE_dataset, aes(x = Cu, y = Cu_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Cu_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Cu_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Cu_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Cu

# Define p value, equation & R2 as a string to add to plots

ACE_Cu_lm_p <- function(ACE_Cu_lm) {
  f <- summary(ACE_Cu_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_lm_p(ACE_Cu_lm)

ACE_Cu_lm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_lm_p(ACE_Cu_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Cu_wlm_p <- function(ACE_Cu_wlm) {
  f <- summary(ACE_Cu_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wlm_p(ACE_Cu_wlm)

ACE_Cu_wlm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wlm_p(ACE_Cu_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Cu_wlm_p <- function(ACE_Cu_wlm) {
  f <- summary(ACE_Cu_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wlm_p(ACE_Cu_wlm)

ACE_Cu_wlm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/(Cu_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wlm_p(ACE_Cu_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Cu_wls_p <- function(ACE_Cu_wls) {
  f <- summary(ACE_Cu_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wls_p(ACE_Cu_wls)

ACE_Cu_wls_eqn <- function(df){
  m <- ACE_Cu_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wls_p(ACE_Cu_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Cu_wls_wt_p <- function(ACE_Cu_wls_wt) {
  f <- summary(ACE_Cu_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wls_wt_p(ACE_Cu_wls_wt)

ACE_Cu_wls_wt_eqn <- function(df){
  m <- ACE_Cu_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wls_wt_p(ACE_Cu_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Cu_final <- ACE_Cu + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Cu_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Cu_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Cu_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Cu_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Cu_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Cu_x.reg <- ACE_dataset$Cu
Cu_y.reg <- ACE_dataset$Cu_ICP
# Build linear regression model
Cu_model1 <- lm(Cu_y.reg ~ Cu, data = ACE_dataset) # OSL unweighted
Cu_model2 <- lm(Cu_y.reg ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt) # WLS unweighted
# Model summary statistics
Cu_model1
summary(Cu_model1)
confint(Cu_model1)
Cu_model2
summary(Cu_model2)
confint(Cu_model2)
# Create predicted values and upper/lower CI to check model is working
Cu_new1.y <- data.frame(Cu_x.reg = c(1, 2, 3))
predict(Cu_model1, Cu_newdata1 = Cu_new1.y)
predict(Cu_model1, Cu_newdata1 = Cu_new1.y, Cu_interval1 = "confidence")
Cu_new2.y <- data.frame(Cu_x.reg = c(1, 2, 3))
predict(Cu_model2, Cu_newdata2 = Cu_new2.y)
predict(Cu_model2, Cu_newdata2 = Cu_new2.y, Cu_interval2 = "confidence")
# Add prediction intervals to model data frame
Cu_pred.int1 <- predict(Cu_model1, interval = "prediction")
Cu_data_1_out <- bind_cols(ACE_dataset, Cu_pred.int1) %>%
  rename(Cu_fit_OLS = fit, Cu_lwr_OLS = lwr, Cu_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Cu, Cu_sd, Cu_ICP, Cu_ICP_sd, Cu_fit_OLS, Cu_lwr_OLS, Cu_upr_OLS))
Cu_pred.int2 <- predict(Cu_model2, interval = "prediction")
Cu_data_2_out <- bind_cols(Cu_data_1_out, Cu_pred.int2) %>%
  rename(Cu_fit_WLS = fit, Cu_lwr_WLS = lwr, Cu_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Cu, Cu_sd, Cu_ICP, Cu_ICP_sd, Cu_fit_OLS, Cu_lwr_OLS, Cu_upr_OLS, Cu_fit_WLS, Cu_lwr_WLS, Cu_upr_WLS)) %>% 
  print()
write.csv(Cu_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Cu_predict <- ACE_Cu_final + 
  geom_line(data = Cu_data_1_out, aes(y = Cu_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Cu_data_1_out, aes(y = Cu_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Cu_data_2_out, aes(y = Cu_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Cu_data_2_out, aes(y = Cu_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Cu_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Cu/Cu_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Zn LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Zn_lm <- lm(Zn_ICP ~ Zn, data = ACE_dataset)
summary(ACE_Zn_lm)
glance(ACE_Zn_lm)
model_performance(ACE_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Zn_wlm <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/(Zn_ICP_sd)^2)
summary(ACE_Zn_wlm)
glance(ACE_Zn_wlm)
model_performance(ACE_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Zn_model <- lm(Zn_ICP ~ Zn, data = ACE_dataset) # define model
ACE_Zn_wt <- 1 / lm(abs(ACE_Zn_model$residuals) ~ ACE_Zn_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Zn_wls <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zn_wls) # summary stats
glance(ACE_Zn_wls) # summary stats including AIC
model_performance(ACE_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Zn_model_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/Zn_ICP_sd^2) # define model
ACE_Zn_wt_wt <- 1 / lm(abs(ACE_Zn_model_wt$residuals) ~ ACE_Zn_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Zn_wls_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zn_wls_wt) # summary stats
glance(ACE_Zn_wls_wt) # summary stats including AIC
model_performance(ACE_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_lm_hats <- as.data.frame(hatvalues(ACE_Zn_lm))
ACE_Zn_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_lm_cooksD <- cooks.distance(ACE_Zn_lm)
ACE_Zn_lm_influential <- ACE_Zn_lm_cooksD[(ACE_Zn_lm_cooksD > (3 * mean(ACE_Zn_lm_cooksD, na.rm = TRUE)))]
ACE_Zn_lm_influential
ACE_Zn_lm_influential_names <- names(ACE_Zn_lm_influential)
ACE_Zn_lm_outliers <- ACE_dataset[ACE_Zn_lm_influential_names,] # outliers only using of index values
ACE_Zn_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zn_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_wlm_hats <- as.data.frame(hatvalues(ACE_Zn_wlm))
ACE_Zn_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_wlm_cooksD <- cooks.distance(ACE_Zn_wlm)
ACE_Zn_wlm_influential <- ACE_Zn_wlm_cooksD[(ACE_Zn_wlm_cooksD > (3 * mean(ACE_Zn_wlm_cooksD, na.rm = TRUE)))]
ACE_Zn_wlm_influential
ACE_Zn_wlm_influential_names <- names(ACE_Zn_wlm_influential)
ACE_Zn_wlm_outliers <- ACE_dataset[ACE_Zn_wlm_influential_names,] # outliers only using of index values
ACE_Zn_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zn_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_wls_hats <- as.data.frame(hatvalues(ACE_Zn_wls))
ACE_Zn_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_wls_cooksD <- cooks.distance(ACE_Zn_wls)
ACE_Zn_wls_influential <- ACE_Zn_wls_cooksD[(ACE_Zn_wls_cooksD > (3 * mean(ACE_Zn_wls_cooksD, na.rm = TRUE)))]
ACE_Zn_wls_influential
ACE_Zn_wls_influential_names <- names(ACE_Zn_wls_influential)
ACE_Zn_wls_outliers <- ACE_dataset[ACE_Zn_wls_influential_names,] # outliers only using of index values
ACE_Zn_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zn_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_wls_wt_hats <- as.data.frame(hatvalues(ACE_Zn_wls_wt))
ACE_Zn_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_wls_wt_cooksD <- cooks.distance(ACE_Zn_wls_wt)
ACE_Zn_wls_wt_influential <- ACE_Zn_wls_wt_cooksD[(ACE_Zn_wls_wt_cooksD > (3 * mean(ACE_Zn_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Zn_wls_wt_influential
ACE_Zn_wls_wt_influential_names <- names(ACE_Zn_wls_wt_influential)
ACE_Zn_wls_wt_outliers <- ACE_dataset[ACE_Zn_wls_wt_influential_names,] # outliers only using of index values
ACE_Zn_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Zn_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_OLS_summary.txt")
summary(ACE_Zn_lm)
glance(ACE_Zn_lm)
model_performance(ACE_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_lm_lev.pdf")
leveragePlots(ACE_Zn_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_OLS_wt_summary.txt")
summary(ACE_Zn_wlm)
glance(ACE_Zn_wlm)
model_performance(ACE_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wlm_lev.pdf")
leveragePlots(ACE_Zn_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_WLS_summary.txt")
summary(ACE_Zn_wls)
glance(ACE_Zn_wls)
model_performance(ACE_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wls_lev.pdf")
leveragePlots(ACE_Zn_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_WLS_wt_summary.txt")
summary(ACE_Zn_wls_wt)
glance(ACE_Zn_wls_wt)
model_performance(ACE_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wls_wt_lev.pdf")
leveragePlots(ACE_Zn_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/ACE_Zn_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Zn"
Zn_WLS_wt = ACE_Zn_wt
Zn_WLS_err_wt = ACE_Zn_wt_wt
theme_set(theme_classic(10))

ACE_Zn <- ggplot(ACE_dataset, aes(x = Zn, y = Zn_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Zn_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Zn_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Zn_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Zn

# Define p value, equation & R2 as a string to add to plots

ACE_Zn_lm_p <- function(ACE_Zn_lm) {
  f <- summary(ACE_Zn_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_lm_p(ACE_Zn_lm)

ACE_Zn_lm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_lm_p(ACE_Zn_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Zn_wlm_p <- function(ACE_Zn_wlm) {
  f <- summary(ACE_Zn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wlm_p(ACE_Zn_wlm)

ACE_Zn_wlm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wlm_p(ACE_Zn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Zn_wlm_p <- function(ACE_Zn_wlm) {
  f <- summary(ACE_Zn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wlm_p(ACE_Zn_wlm)

ACE_Zn_wlm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/(Zn_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wlm_p(ACE_Zn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Zn_wls_p <- function(ACE_Zn_wls) {
  f <- summary(ACE_Zn_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wls_p(ACE_Zn_wls)

ACE_Zn_wls_eqn <- function(df){
  m <- ACE_Zn_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wls_p(ACE_Zn_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Zn_wls_wt_p <- function(ACE_Zn_wls_wt) {
  f <- summary(ACE_Zn_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wls_wt_p(ACE_Zn_wls_wt)

ACE_Zn_wls_wt_eqn <- function(df){
  m <- ACE_Zn_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wls_wt_p(ACE_Zn_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Zn_final <- ACE_Zn + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Zn_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Zn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Zn_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Zn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Zn_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Zn_x.reg <- ACE_dataset$Zn
Zn_y.reg <- ACE_dataset$Zn_ICP
# Build linear regression model
Zn_model1 <- lm(Zn_y.reg ~ Zn, data = ACE_dataset) # OSL unweighted
Zn_model2 <- lm(Zn_y.reg ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt) # WLS unweighted
# Model summary statistics
Zn_model1
summary(Zn_model1)
confint(Zn_model1)
Zn_model2
summary(Zn_model2)
confint(Zn_model2)
# Create predicted values and upper/lower CI to check model is working
Zn_new1.y <- data.frame(Zn_x.reg = c(1, 2, 3))
predict(Zn_model1, Zn_newdata1 = Zn_new1.y)
predict(Zn_model1, Zn_newdata1 = Zn_new1.y, Zn_interval1 = "confidence")
Zn_new2.y <- data.frame(Zn_x.reg = c(1, 2, 3))
predict(Zn_model2, Zn_newdata2 = Zn_new2.y)
predict(Zn_model2, Zn_newdata2 = Zn_new2.y, Zn_interval2 = "confidence")
# Add prediction intervals to model data frame
Zn_pred.int1 <- predict(Zn_model1, interval = "prediction")
Zn_data_1_out <- bind_cols(ACE_dataset, Zn_pred.int1) %>%
  rename(Zn_fit_OLS = fit, Zn_lwr_OLS = lwr, Zn_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Zn, Zn_sd, Zn_ICP, Zn_ICP_sd, Zn_fit_OLS, Zn_lwr_OLS, Zn_upr_OLS))
Zn_pred.int2 <- predict(Zn_model2, interval = "prediction")
Zn_data_2_out <- bind_cols(Zn_data_1_out, Zn_pred.int2) %>%
  rename(Zn_fit_WLS = fit, Zn_lwr_WLS = lwr, Zn_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Zn, Zn_sd, Zn_ICP, Zn_ICP_sd, Zn_fit_OLS, Zn_lwr_OLS, Zn_upr_OLS, Zn_fit_WLS, Zn_lwr_WLS, Zn_upr_WLS)) %>% 
  print()
write.csv(Zn_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Zn_predict <- ACE_Zn_final + 
  geom_line(data = Zn_data_1_out, aes(y = Zn_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Zn_data_1_out, aes(y = Zn_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Zn_data_2_out, aes(y = Zn_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Zn_data_2_out, aes(y = Zn_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Zn_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Zn/Zn_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Rb LM  -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Rb_lm <- lm(Rb_ICP ~ Rb, data = ACE_dataset)
summary(ACE_Rb_lm)
glance(ACE_Rb_lm)
model_performance(ACE_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Rb_wlm <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/(Rb_ICP_sd)^2)
summary(ACE_Rb_wlm)
glance(ACE_Rb_wlm)
model_performance(ACE_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Rb_model <- lm(Rb_ICP ~ Rb, data = ACE_dataset) # define model
ACE_Rb_wt <- 1 / lm(abs(ACE_Rb_model$residuals) ~ ACE_Rb_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Rb_wls <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt) #perform weighted least squares regression
# Checks
summary(ACE_Rb_wls) # summary stats
glance(ACE_Rb_wls) # summary stats including AIC
model_performance(ACE_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted - model
ACE_Rb_model_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/Rb_ICP_sd^2) # define model
ACE_Rb_wt_wt <- 1 / lm(abs(ACE_Rb_model_wt$residuals) ~ ACE_Rb_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Rb_wls_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Rb_wls_wt) # summary stats
glance(ACE_Rb_wls_wt) # summary stats including AIC
model_performance(ACE_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_lm_hats <- as.data.frame(hatvalues(ACE_Rb_lm))
ACE_Rb_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_lm_cooksD <- cooks.distance(ACE_Rb_lm)
ACE_Rb_lm_influential <- ACE_Rb_lm_cooksD[(ACE_Rb_lm_cooksD > (3 * mean(ACE_Rb_lm_cooksD, na.rm = TRUE)))]
ACE_Rb_lm_influential
ACE_Rb_lm_influential_names <- names(ACE_Rb_lm_influential)
ACE_Rb_lm_outliers <- ACE_dataset[ACE_Rb_lm_influential_names,] # outliers only using of index values
ACE_Rb_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Rb_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_wlm_hats <- as.data.frame(hatvalues(ACE_Rb_wlm))
ACE_Rb_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_wlm_cooksD <- cooks.distance(ACE_Rb_wlm)
ACE_Rb_wlm_influential <- ACE_Rb_wlm_cooksD[(ACE_Rb_wlm_cooksD > (3 * mean(ACE_Rb_wlm_cooksD, na.rm = TRUE)))]
ACE_Rb_wlm_influential
ACE_Rb_wlm_influential_names <- names(ACE_Rb_wlm_influential)
ACE_Rb_wlm_outliers <- ACE_dataset[ACE_Rb_wlm_influential_names,] # outliers only using of index values
ACE_Rb_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Rb_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_wls_hats <- as.data.frame(hatvalues(ACE_Rb_wls))
ACE_Rb_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_wls_cooksD <- cooks.distance(ACE_Rb_wls)
ACE_Rb_wls_influential <- ACE_Rb_wls_cooksD[(ACE_Rb_wls_cooksD > (3 * mean(ACE_Rb_wls_cooksD, na.rm = TRUE)))]
ACE_Rb_wls_influential
ACE_Rb_wls_influential_names <- names(ACE_Rb_wls_influential)
ACE_Rb_wls_outliers <- ACE_dataset[ACE_Rb_wls_influential_names,] # outliers only using of index values
ACE_Rb_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Rb_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_wls_wt_hats <- as.data.frame(hatvalues(ACE_Rb_wls_wt))
ACE_Rb_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_wls_wt_cooksD <- cooks.distance(ACE_Rb_wls_wt)
ACE_Rb_wls_wt_influential <- ACE_Rb_wls_wt_cooksD[(ACE_Rb_wls_wt_cooksD > (3 * mean(ACE_Rb_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Rb_wls_wt_influential
ACE_Rb_wls_wt_influential_names <- names(ACE_Rb_wls_wt_influential)
ACE_Rb_wls_wt_outliers <- ACE_dataset[ACE_Rb_wls_wt_influential_names,] # outliers only using of index values
ACE_Rb_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_Rb_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_OLS_summary.txt")
summary(ACE_Rb_lm)
glance(ACE_Rb_lm)
model_performance(ACE_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_lm_lev.pdf")
leveragePlots(ACE_Rb_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_OLS_wt_summary.txt")
summary(ACE_Rb_wlm)
glance(ACE_Rb_wlm)
model_performance(ACE_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wlm_lev.pdf")
leveragePlots(ACE_Rb_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_WLS_summary.txt")
summary(ACE_Rb_wls)
glance(ACE_Rb_wls)
model_performance(ACE_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wls_lev.pdf")
leveragePlots(ACE_Rb_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_WLS_wt_summary.txt")
summary(ACE_Rb_wls_wt)
glance(ACE_Rb_wls_wt)
model_performance(ACE_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wls_wt_lev.pdf")
leveragePlots(ACE_Rb_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/ACE_Rb_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Rb"
Rb_WLS_wt = ACE_Rb_wt
Rb_WLS_err_wt = ACE_Rb_wt_wt
theme_set(theme_classic(10))

ACE_Rb <- ggplot(ACE_dataset, aes(x = Rb, y = Rb_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              linetype = "dashed", aes(weight = 1/Rb_ICP_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = Rb_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Rb_WLS_wt), colour="black") + # WLS unweighted
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title,  itrax_dataset, "[XRF-CS]") , y = paste0(element_title, " clr [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title, " clr"))
ACE_Rb

# Define p value, equation & R2 as a string to add to plots

ACE_Rb_lm_p <- function(ACE_Rb_lm) {
  f <- summary(ACE_Rb_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_lm_p(ACE_Rb_lm)

ACE_Rb_lm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_lm_p(ACE_Rb_lm), digits = 2)))
  as.character(as.expression(eq));
}
# Define p value, OLS equation & R2 as a string to add to plot

ACE_Rb_wlm_p <- function(ACE_Rb_wlm) {
  f <- summary(ACE_Rb_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wlm_p(ACE_Rb_wlm)

ACE_Rb_wlm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wlm_p(ACE_Rb_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_Rb_wlm_p <- function(ACE_Rb_wlm) {
  f <- summary(ACE_Rb_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wlm_p(ACE_Rb_wlm)

ACE_Rb_wlm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/(Rb_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wlm_p(ACE_Rb_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_Rb_wls_p <- function(ACE_Rb_wls) {
  f <- summary(ACE_Rb_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wls_p(ACE_Rb_wls)

ACE_Rb_wls_eqn <- function(df){
  m <- ACE_Rb_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wls_p(ACE_Rb_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_Rb_wls_wt_p <- function(ACE_Rb_wls_wt) {
  f <- summary(ACE_Rb_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wls_wt_p(ACE_Rb_wls_wt)

ACE_Rb_wls_wt_eqn <- function(df){
  m <- ACE_Rb_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wls_wt_p(ACE_Rb_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Rb_final <- ACE_Rb + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Rb_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Rb_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Rb_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Rb_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Rb_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Rb_x.reg <- ACE_dataset$Rb
Rb_y.reg <- ACE_dataset$Rb_ICP
# Build linear regression model
Rb_model1 <- lm(Rb_y.reg ~ Rb, data = ACE_dataset) # OSL unweighted
Rb_model2 <- lm(Rb_y.reg ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt) # WLS unweighted
# Model summary statistics
Rb_model1
summary(Rb_model1)
confint(Rb_model1)
Rb_model2
summary(Rb_model2)
confint(Rb_model2)
# Create predicted values and upper/lower CI to check model is working
Rb_new1.y <- data.frame(Rb_x.reg = c(1, 2, 3))
predict(Rb_model1, Rb_newdata1 = Rb_new1.y)
predict(Rb_model1, Rb_newdata1 = Rb_new1.y, Rb_interval1 = "confidence")
Rb_new2.y <- data.frame(Rb_x.reg = c(1, 2, 3))
predict(Rb_model2, Rb_newdata2 = Rb_new2.y)
predict(Rb_model2, Rb_newdata2 = Rb_new2.y, Rb_interval2 = "confidence")
# Add prediction intervals to model data frame
Rb_pred.int1 <- predict(Rb_model1, interval = "prediction")
Rb_data_1_out <- bind_cols(ACE_dataset, Rb_pred.int1) %>%
  rename(Rb_fit_OLS = fit, Rb_lwr_OLS = lwr, Rb_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Rb, Rb_sd, Rb_ICP, Rb_ICP_sd, Rb_fit_OLS, Rb_lwr_OLS, Rb_upr_OLS))
Rb_pred.int2 <- predict(Rb_model2, interval = "prediction")
Rb_data_2_out <- bind_cols(Rb_data_1_out, Rb_pred.int2) %>%
  rename(Rb_fit_WLS = fit, Rb_lwr_WLS = lwr, Rb_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Rb, Rb_sd, Rb_ICP, Rb_ICP_sd, Rb_fit_OLS, Rb_lwr_OLS, Rb_upr_OLS, Rb_fit_WLS, Rb_lwr_WLS, Rb_upr_WLS)) %>% 
  print()
write.csv(Rb_data_2_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Rb_predict <- ACE_Rb_final + 
  geom_line(data = Rb_data_1_out, aes(y = Rb_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Rb_data_1_out, aes(y = Rb_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Rb_data_2_out, aes(y = Rb_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Rb_data_2_out, aes(y = Rb_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: " , element_title, " clr [95% CI & PI]"))
ACE_Rb_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/Rb/Rb_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# ACE_Dry Mass (DM) -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_DM_lm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset)
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_OLS_summary.txt")
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_DM_wlm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/(coh_inc_sd)^2)
summary(ACE_DM_wlm)
glance(ACE_DM_wlm)
model_performance(ACE_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_OLS_wt_summary.txt")
summary(ACE_DM_wlm)
glance(ACE_DM_wlm)
model_performance(ACE_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
ACE_DM_model <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset) # define model
ACE_DM_wt <- 1 / lm(abs(ACE_DM_model$residuals) ~ ACE_DM_model$fitted.values)$fitted.values^2 #define weights to use
ACE_DM_wls <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=ACE_DM_wt) #perform weighted least squares regression
# Checks
summary(ACE_DM_wls) # summary stats
glance(ACE_DM_wls) # summary stats including AIC
model_performance(ACE_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_WLS_summary.txt")
summary(ACE_DM_wls)
glance(ACE_DM_wls)
model_performance(ACE_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_wls) # Performance package summary check for heteroscedasticity
icc(ACE_DM_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_DM_model_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/coh_inc_sd^2) # define model
ACE_DM_wt_wt <- 1 / lm(abs(ACE_DM_model_wt$residuals) ~ ACE_DM_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_DM_wls_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=ACE_DM_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_DM_wls_wt) # summary stats
glance(ACE_DM_wls_wt) # summary stats including AIC
model_performance(ACE_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_WLS_wt_summary.txt")
summary(ACE_DM_wls_wt)
glance(ACE_DM_wls_wt)
model_performance(ACE_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_DM_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance

# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_DM_lm_hats <- as.data.frame(hatvalues(ACE_DM_lm))
ACE_DM_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_DM_lm_cooksD <- cooks.distance(ACE_DM_lm)
ACE_DM_lm_influential <- ACE_DM_lm_cooksD[(ACE_DM_lm_cooksD > (3 * mean(ACE_DM_lm_cooksD, na.rm = TRUE)))]
ACE_DM_lm_influential
ACE_DM_lm_influential_names <- names(ACE_DM_lm_influential)
ACE_DM_lm_outliers <- ACE_dataset[ACE_DM_lm_influential_names,] # outliers only using of index values
ACE_DM_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_DM_lm_lev.pdf")
leveragePlots(ACE_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_DM_wlm_hats <- as.data.frame(hatvalues(ACE_DM_wlm))
ACE_DM_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_DM_wlm_cooksD <- cooks.distance(ACE_DM_wlm)
ACE_DM_wlm_influential <- ACE_DM_wlm_cooksD[(ACE_DM_wlm_cooksD > (3 * mean(ACE_DM_wlm_cooksD, na.rm = TRUE)))]
ACE_DM_wlm_influential
ACE_DM_wlm_influential_names <- names(ACE_DM_wlm_influential)
ACE_DM_wlm_outliers <- ACE_dataset[ACE_DM_wlm_influential_names,] # outliers only using of index values
ACE_DM_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_DM_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/OLS_ACE_DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_DM_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_DM_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_DM_wlm_lev.pdf")
leveragePlots(ACE_DM_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_DM_wls_hats <- as.data.frame(hatvalues(ACE_DM_wls))
ACE_DM_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_DM_wls_cooksD <- cooks.distance(ACE_DM_wls)
ACE_DM_wls_influential <- ACE_DM_wls_cooksD[(ACE_DM_wls_cooksD > (3 * mean(ACE_DM_wls_cooksD, na.rm = TRUE)))]
ACE_DM_wls_influential
ACE_DM_wls_influential_names <- names(ACE_DM_wls_influential)
ACE_DM_wls_outliers <- ACE_dataset[ACE_DM_wls_influential_names,] # outliers only using of index values
ACE_DM_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_DM_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_wls_lev_bar.pdf")
barplot(hatvalues(ACE_DM_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_wls_lev.pdf")
leveragePlots(ACE_DM_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_DM_wls_wt_hats <- as.data.frame(hatvalues(ACE_DM_wls_wt))
ACE_DM_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_DM_wls_wt_cooksD <- cooks.distance(ACE_DM_wls_wt)
ACE_DM_wls_wt_influential <- ACE_DM_wls_wt_cooksD[(ACE_DM_wls_wt_cooksD > (3 * mean(ACE_DM_wls_wt_cooksD, na.rm = TRUE)))]
ACE_DM_wls_wt_influential
ACE_DM_wls_wt_influential_names <- names(ACE_DM_wls_wt_influential)
ACE_DM_wls_wt_outliers <- ACE_dataset[ACE_DM_wls_wt_influential_names,] # outliers only using of index values
ACE_DM_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_DM_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_DM_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_wls_wt_lev.pdf")
leveragePlots(ACE_DM_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/ACE_ACE_DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_wls_wt)
dev.off()

# Linear model summary elemental plots
element_title <- "coh_inc"
theme_set(theme_classic(10))
ACE_DM <- ggplot(ACE_dataset, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/coh_inc_ICP_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = ACE_DM_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = ACE_DM_wt_wt), colour="black") + # weighted WLS regression
  #scale_shape_manual(values = c(21)) +
  #scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, "[XRF-CS]") , y = paste0("Dry mass (%)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", element_title))
ACE_DM
# Define p value, OLS equation & R2 as a string to add to plot

ACE_DM_lm_p <- function(ACE_DM_lm) {
  f <- summary(ACE_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_lm_p(ACE_DM_lm)

ACE_DM_lm_eqn <- function(df){
  m <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_lm_p(ACE_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_DM_wlm_p <- function(ACE_DM_wlm) {
  f <- summary(ACE_DM_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_wlm_p(ACE_DM_wlm)

ACE_DM_wlm_eqn <- function(df){
  m <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/(coh_inc_ICP_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_wlm_p(ACE_DM_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_DM_wls_p <- function(ACE_DM_wls) {
  f <- summary(ACE_DM_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_wls_p(ACE_DM_wls)

ACE_DM_wls_eqn <- function(df){
  m <- ACE_DM_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_wls_p(ACE_DM_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_DM_wls_wt_p <- function(ACE_DM_wls_wt) {
  f <- summary(ACE_DM_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_wls_wt_p(ACE_DM_wls_wt)

ACE_DM_wls_wt_eqn <- function(df){
  m <- ACE_DM_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_wls_wt_p(ACE_DM_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_DM_final <- ACE_DM + 
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_DM_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_DM_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_DM_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_DM_final
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/DM/DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# Supplementary Figure - All Matched Elements -------------------------------------------------------------

# Summary 3x2 matrix plots of ITRAX-acf 2ndry matched elements
# OLS & WLS summary - unweighted stats on plot
ggarrange(ACE_K_predict, ACE_Co_predict, ACE_Ni_predict, 
          ACE_Cu_predict, ACE_Zn_predict, ACE_Rb_predict,
          ncol = 2, nrow = 3, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/ACE_Fig3a_OLS_WLS_log_inc_2ndry_predict.pdf",
       height = c(40), width = c(30), dpi = 600, units = "cm")

# Summary 4x3 matrix plots of ITRAX-acf matched elements - no Co - replaced by DM
# OLS & WLS summary - unweighted stats on plot
ggarrange(ACE_K_predict, ACE_Ca_predict, ACE_Ti_predict,
          ACE_Mn_predict, ACE_Fe_predict, ACE_Ni_predict, 
          ACE_Cu_predict, ACE_Zn_predict, ACE_Rb_predict, 
          ACE_Sr_predict, ACE_Zr_predict, ACE_DM_final,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3b_OLS_WLS_clr/ACE_FigS4b_OLS_WLS_clr_predict.pdf",
       height = c(50), width = c(50), dpi = 600, units = "cm")

# END ---------------------------------------------------------------------



