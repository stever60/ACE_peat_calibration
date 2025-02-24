# ITRAX-ICPMS Matching of ACE data - KER3 element/inc mean +/- SD

# Set up -----------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()

# Load libraries ---------------------------------------------------------------

packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'readr', 'ggpubr', 'patchwork',
              'gridExtra', 'cowplot', 'vegan', 'rioja', 'ellipse', 'factoextra',
              'reshape2', 'GGally', 'ggsci', 'ggdendro', 'dendextend', 'dynamicTreeCut',
              'colorspace', 'cluster', 'magrittr', 'mgcv', 'gtable', 'repr',
              'bestNormalize','sjmisc', 'chemometrics', 'compositions', 
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis',
              'ggrepel', 'itraxR','PeriodicTable','errors','patchwork',
              'forecast','directlabels','broom','performance','lmtest','ggpmisc',
              'cowplot','Hmisc')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

# Set up ------------------------------------------------------------------

# Set working directory - Macbook Pro 2013
#setwd("/Users/Steve/Dropbox/BAS/Data/R/")
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()
# Element lists  ----------------------------------------------------------

# Defined from acfanalysis of ACE09 composite ITRAX dataframe analysis file:
# https://www.dropbox.com/s/rhtlkp6uwp71ryc/ACE_ITRAX_COMPOSITE.R?dl=0

# ICP elements as defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")

# XRF elements defined by Francois & ITRAX acf
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

# ICP elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

xrf_icp_Elements_min <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                          "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                          "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                          "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                          "Mo_inc", "Mo_coh", "coh_inc", "dry_mass_pc")
                        
# Matching ----------------------------------------------------------------

# ITRAX to ICP dataframe usind min-max depths

# ICPMS 

# define element list to use (ITRAX autocorrelation (acf) defined & Francois ICP matched)
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")

icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# Import ICPMS composite dataframe

ACE_ICP <- read_csv("Papers_R/2024_DeVleeschouwer/Data/ACE_SHW_ICPMS_Composite.csv", 
                      col_names = TRUE, skip = 0) %>% 
  rename(sample = Sample_ID)
ACE_ICP

# Import KER3 ICP dataframe - use find/replace for other sites
ACE_ICP_KER3 <- ACE_ICP %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all inf to NA
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "KER3") %>% # select site
  rename(top = field_depth_top) %>% # rename depths for matching
  rename(bottom = field_depth_bottom) %>%
  rename(mp = field_depth_mid) %>% 
  print()
write.csv(ACE_ICP_KER3,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/ACE_ICP_KER3.csv", 
          row.names = FALSE)

# Import ITRAX QC & ACF cps dataframe
ACE_itrax <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE_itrax.csv",
         col_names = TRUE, skip = 0) 
ACE_itrax

# ITRAX - cps -------------------------------------------------------------

# # Select KER3 site data - use find/replace for other sites
# ACE_itrax_KER3_cps <- ACE_itrax %>% 
#   filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
#   filter(Site == "KER3") %>% 
#   filter(qc == TRUE) %>%
#   select(depth:surface, kcps:MSE, all_of(acf_icp_Elements_min), coh_inc) %>% 
#   print()
# write.csv(ACE_itrax_KER3_cps,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/ACE_itrax_KER3_cps.csv", 
#           row.names = FALSE)

# ITRAX - Inc normalised --------------------------------------------------

# Produce E/inc or Ln(E/inc) dataset
# Replace 0 values in each column with half column min for log function to work
ACE_itrax_KER3_inc <- ACE_itrax %>%
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "KER3") %>%
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:MSE, all_of(acf_icp_Elements_min), coh_inc) %>%
  mutate(across(any_of(acf_icp_Elements_min), .fns = ~./ Mo_inc)) %>% # uncomment to run as inc normalised data
  mutate_at(vars(any_of(acf_icp_Elements_min)), ~ (. == 0) * min(.[. != 0])/2 + .) # %>%
# mutate(across(any_of(acf_icp_Elements_min), log))  # natural log of inc normalised values
is.na(ACE_itrax_KER3_inc)<-sapply(ACE_itrax_KER3_inc, is.infinite)
ACE_itrax_KER3_inc
write.csv(ACE_itrax_KER3_inc,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/ACE_itrax_KER3_inc.csv",
          row.names = FALSE)

# # ITRAX - %cps sum  ----------------------------
# # cps as % of cps_sum - recalculated for acf-icp elements
# ACE_itrax_KER3_cps_sum <- ACE_itrax %>%
#   filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
#   filter(Site == "KER3") %>%
#   filter(qc == TRUE) %>%
#   select(depth:surface, kcps:MSE, all_of(acf_icp_Elements_min), coh_inc) %>%
#   mutate(cps_acf_sum = rowSums(across(all_of(acf_icp_Elements_min)))) %>%
#   mutate(across(all_of(acf_icp_Elements_min)) / cps_acf_sum) %>%
#   mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
#   replace(is.na(.), 0) %>%
#   mutate(across(all_of(acf_icp_Elements_min)) *100) %>%
#   mutate(sum = rowSums(across(all_of(acf_icp_Elements_min)))) %>%
#   mutate(scatter_sum = Mo_inc+Mo_coh) %>%
#   relocate(cps_acf_sum, .after = cps) %>% 
#   print()
# # check sum <100% & write file
# head(ACE_itrax_KER3_cps_sum$scatter_sum)
# write.csv(ACE_itrax_KER3_cps_sum,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/pc_cps/ACE_itrax_KER3_cps_sum.csv", row.names = FALSE)
# 
# # ITRAX - clr -------------------------------------------------------------
# 
# library (compositions)
# library(tidyverse)
# 
# # Produce clr dataframe - replace 0 values in each column with half column min for clr function to work
# ACE_itrax_KER3_clr1 <- ACE_itrax %>%
#   filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
#   filter(Site == "KER3") %>% # select site KER3 data only
#   filter(qc == TRUE) %>% # select only spectra that pass itrax.R based quality control (qc)
#   select(all_of(acf_icp_Elements_min)) %>% # select elements to use - all acf elements use:
#   mutate_at(vars(any_of(acf_icp_Elements_min)), ~ (. == 0) * min(.[. != 0])/2 + .) %>% # replace 0 values in each column with half column min
#   clr() %>% # make centred log ratio dataframe
#   print() 
# is.na(ACE_itrax_KER3_clr1)<-sapply(ACE_itrax_KER3_clr1, is.infinite) # replace any infinite values with NA
# 
# ACE_itrax_KER3_text <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE_itrax.csv",
#                                   col_names = TRUE, skip = 0) %>%
#   filter(Site == "KER3") %>%
#   filter(qc == TRUE) %>%
#   select(depth:surface, kcps:MSE, inc_coh, coh_inc) %>% # add Location:Sitehere if not being used in matching
# 
# ACE_itrax_KER3_clr <- bind_cols(ACE_itrax_KER3_text, ACE_itrax_KER3_clr1) %>%
#   relocate(inc_coh, .after = Mo_coh) %>%
#   relocate(coh_inc, .after = inc_coh) %>% 
#   print()
# write.csv(ACE_itrax_KER3_clr,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/clr/ACE_itrax_KER3_clr.csv",
#           row.names = FALSE)

# MATCHING ----------------------------------------------------------------

# define adjacent matched elements
xrf_icp_Elements_min <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                          "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                          "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                          "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                          "Mo_inc", "Mo_coh", "coh_inc", "dry_mass_pc")

# Choose dataset to use for matching
#ACE_itrax_KER3 <- ACE_itrax_KER3_cps
ACE_itrax_KER3 <- ACE_itrax_KER3_inc # use find/replace cps to run inc normalised data
#ACE_itrax_KER3 <- ACE_itrax_KER3_cps_sum
#ACE_itrax_KER3 <- ACE_itrax_KER3_clr

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
KER3_itrax_matched <- itrax_reduce(ACE_itrax_KER3, names = ACE_ICP_KER3$sample,
                                        breaks_lower = ACE_ICP_KER3$top,
                                        breaks_upper = ACE_ICP_KER3$bottom,
                                        fun = mean,
                                        edges = c(">=", "<=")) %>%
  select(resample_names, SH20_age, kcps:MSE, all_of(acf_icp_Elements_min), coh_inc) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_ICP_KER3, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  select(sample, min_depth, max_depth, midpoint, SH20_age, kcps:MSE, all_of(acf_icp_Elements_min), coh_inc)
KER3_itrax_matched

# Calculate stdev for ITRAX matching data
KER3_itrax_matched_sd <- itrax_reduce(ACE_itrax_KER3, names = ACE_ICP_KER3$sample,
                                     breaks_lower = ACE_ICP_KER3$top,
                                     breaks_upper = ACE_ICP_KER3$bottom,
                                     fun = sd,
                                     edges = c(">=", "<=")) %>%
  select(resample_names, SH20_age, kcps:MSE, all_of(acf_icp_Elements_min), coh_inc) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_ICP_KER3, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  select(sample, all_of(acf_icp_Elements_min), coh_inc) %>% 
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd)
KER3_itrax_matched_sd

# Merge mean & sd data into the same dataframe
KER3_xrf_matched <-  KER3_itrax_matched %>% 
  inner_join(., KER3_itrax_matched_sd, by = "sample")
KER3_xrf_matched

# Create matched ICPMS and ITRAX dataframe for calibration & plotting analysis
KER3_xrf_icp_matched <-  KER3_xrf_matched %>% 
  select(sample, min_depth:MSE, all_of(acf_icp_Elements_min), coh_inc, all_of(acf_icp_Elements_min_sd), coh_inc_sd) %>% 
  inner_join(., ACE_ICP_KER3, by = "sample") %>%
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(Location:mp, .before = sample) %>% 
  select(-c(top:mp)) #%>% 
  #mutate(across(SH20_age:cps, round, 0)) %>% 
  #mutate(across(c(MSE:Mo_coh, K_sd:Mo_coh_sd, water_content_pc:density_gcm3_err), round, 2))
KER3_xrf_icp_matched
write.csv(KER3_xrf_icp_matched,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/KER3_xrf_icp_matched.csv", row.names = FALSE)

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
KER3_xrf_icp_matched <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/KER3_xrf_icp_matched.csv")
is.na(KER3_xrf_icp_matched)<-sapply(KER3_xrf_icp_matched, is.infinite) # replace any infinite values with NA
KER3_xrf_icp_matched <- KER3_xrf_icp_matched %>% 
  mutate_at(vars(all_of(c(icp_Elements_min, acf_icp_Elements_min))), 
            ~ (. == 0) * min(.[. != 0])/2 + .) %>%
  mutate_at(vars(all_of(c(icp_Elements_min_sd, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
KER3_xrf_icp_matched$Zn

# Correlation matrices --------------------------------------------------

# ITRAX
theme_set(theme_bw(base_size=2))
ggcorr(KER3_xrf_icp_matched[,acf_icp_Elements_min], method = c("everything", "pearson"), 
       size = 5, label = TRUE, label_alpha = TRUE, label_round=2, label_size= 5)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Correlation/KER3_itrax_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ICP
theme_set(theme_bw(base_size=2))
ggcorr(KER3_xrf_icp_matched[,icp_Elements_min], method = c("everything", "pearson"), 
       size = 5, label = TRUE, label_alpha = TRUE, label_round=2, label_size= 5)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Correlation/KER3_ICP_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP
theme_set(theme_bw(base_size=2))
ggcorr(KER3_xrf_icp_matched[,xrf_icp_Elements_min], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size= 3)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Correlation/KER3_itrax_ICP_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation & density matrices ----------------------------------------

# ITRAX
theme_set(theme_bw(base_size=8))
ggpairs(KER3_xrf_icp_matched, columns = acf_icp_Elements_min, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Correlation/KER3_itrax_ICP_Corr-den_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ICP
theme_set(theme_bw(base_size=8))
ggpairs(KER3_xrf_icp_matched, columns = icp_Elements_min, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Correlation/KER3_itrax_ICP_Corr-den_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP
theme_set(theme_bw(base_size=8))
ggpairs(KER3_xrf_icp_matched, columns = xrf_icp_Elements_min, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Correlation/KER3_itrax_ICP_Corr-den_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Linear models for ITRAX-ACF matched elements -----------------------------------------------------------

library(performance) # linear model testing and graphical outputs
library(ggpmisc) # this package is used for simple/quick for labelling (stat_poly_eq)
library(lmtest) # linear model testing 

# Linear models & performance assessment /stats -----------------------------------------------

# Set up labels
ACE_dataset <- "KER3_xrf_icp_matched"
site_title <- "KER3"
#dataset <- " (cps)"
dataset <- " (E/inc)"
#dataset <- " [Ln(E/inc)]"
#dataset <- " (%cps sum)"
#dataset <- " (clr)"

# KER3_K -----------------------------------------------------------------------

# Unweighted OLS (Ordinary Least Squares) - linear model & checks
KER3_K_lm <- lm(K_ICP ~ K, data = ACE_dataset)
summary(KER3_K_lm)
glance(KER3_K_lm)
model_performance(KER3_K_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_K_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_OLS_summary.txt")
summary(KER3_K_lm)
glance(KER3_K_lm)
model_performance(KER3_K_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_K_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_K_lm) # Performance package summary check for heteroscedasticity
icc(KER3_K_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_K_wlm <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/(K_sd)^2)
summary(KER3_K_wlm)
glance(KER3_K_wlm)
model_performance(KER3_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_K_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_OLS_wt_summary.txt")
summary(KER3_K_wlm)
glance(KER3_K_wlm)
model_performance(KER3_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_K_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_K_lm) # Performance package summary check for heteroscedasticity
icc(KER3_K_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_K_model <- lm(K_ICP ~ K, data = ACE_dataset) # define model
KER3_K_wt <- 1 / lm(abs(KER3_K_model$residuals) ~ KER3_K_model$fitted.values)$fitted.values^2 #define weights to use
KER3_K_wls <- lm(K_ICP ~ K, data = ACE_dataset, weights=KER3_K_wt) #perform weighted least squares regression
# Checks
summary(KER3_K_wls) # summary stats
glance(KER3_K_wls) # summary stats including AIC
model_performance(KER3_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_K_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_WLS_summary.txt")
summary(KER3_K_wls)
glance(KER3_K_wls)
model_performance(KER3_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_K_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_K_wls) # Performance package summary check for heteroscedasticity
icc(KER3_K_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_K_model_wt <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/K_sd^2) # define model
KER3_K_wt_wt <- 1 / lm(abs(KER3_K_model_wt$residuals) ~ KER3_K_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_K_wls_wt <- lm(K_ICP ~ K, data = ACE_dataset, weights=KER3_K_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_K_wls_wt) # summary stats
glance(KER3_K_wls_wt) # summary stats including AIC
model_performance(KER3_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_K_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_WLS_wt_summary.txt")
summary(KER3_K_wls_wt)
glance(KER3_K_wls_wt)
model_performance(KER3_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_K_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_K_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_K_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "K"
theme_set(theme_bw(10))
KER3_K <- ggplot(ACE_dataset, aes(x = K, y = K_ICP)) +
  geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/K_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_K_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_K_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_K_lm_p <- function(KER3_K_lm) {
  f <- summary(KER3_K_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_K_lm_p(KER3_K_lm)

KER3_K_lm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_K_lm_p(KER3_K_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_K_wlm_p <- function(KER3_K_wlm) {
  f <- summary(KER3_K_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_K_wlm_p(KER3_K_wlm)

KER3_K_wlm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/(K_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_K_wlm_p(KER3_K_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_K_wls_p <- function(KER3_K_wls) {
  f <- summary(KER3_K_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_K_wls_p(KER3_K_wls)

KER3_K_wls_eqn <- function(df){
  m <- KER3_K_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_K_wls_p(KER3_K_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_K_wls_wt_p <- function(KER3_K_wls_wt) {
  f <- summary(KER3_K_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_K_wls_wt_p(KER3_K_wls_wt)

KER3_K_wls_wt_eqn <- function(df){
  m <- KER3_K_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_K_wls_wt_p(KER3_K_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_K_final <- KER3_K + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_K_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_K_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_K_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_K_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_K_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/K/K_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# KER3_Ca -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Ca_lm <- lm(Ca_ICP ~ Ca, data = ACE_dataset)
summary(KER3_Ca_lm)
glance(KER3_Ca_lm)
model_performance(KER3_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ca_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_OLS_summary.txt")
summary(KER3_Ca_lm)
glance(KER3_Ca_lm)
model_performance(KER3_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ca_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ca_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Ca_wlm <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/(Ca_sd)^2)
summary(KER3_Ca_wlm)
glance(KER3_Ca_wlm)
model_performance(KER3_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ca_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_OLS_wt_summary.txt")
summary(KER3_Ca_wlm)
glance(KER3_Ca_wlm)
model_performance(KER3_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ca_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ca_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Ca_model <- lm(Ca_ICP ~ Ca, data = ACE_dataset) # define model
KER3_Ca_wt <- 1 / lm(abs(KER3_Ca_model$residuals) ~ KER3_Ca_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Ca_wls <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=KER3_Ca_wt) #perform weighted least squares regression
# Checks
summary(KER3_Ca_wls) # summary stats
glance(KER3_Ca_wls) # summary stats including AIC
model_performance(KER3_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ca_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_WLS_summary.txt")
summary(KER3_Ca_wls)
glance(KER3_Ca_wls)
model_performance(KER3_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ca_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ca_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Ca_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Ca_model_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/Ca_sd^2) # define model
KER3_Ca_wt_wt <- 1 / lm(abs(KER3_Ca_model_wt$residuals) ~ KER3_Ca_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Ca_wls_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=KER3_Ca_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Ca_wls_wt) # summary stats
glance(KER3_Ca_wls_wt) # summary stats including AIC
model_performance(KER3_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ca_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_WLS_wt_summary.txt")
summary(KER3_Ca_wls_wt)
glance(KER3_Ca_wls_wt)
model_performance(KER3_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ca_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ca_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Ca_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Ca"
theme_set(theme_bw(10))
KER3_Ca <- ggplot(ACE_dataset, aes(x = Ca, y = Ca_ICP)) +
  geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Ca_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Ca_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Ca_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ca_lm_p <- function(KER3_Ca_lm) {
  f <- summary(KER3_Ca_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ca_lm_p(KER3_Ca_lm)

KER3_Ca_lm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ca_lm_p(KER3_Ca_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ca_wlm_p <- function(KER3_Ca_wlm) {
  f <- summary(KER3_Ca_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ca_wlm_p(KER3_Ca_wlm)

KER3_Ca_wlm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/(Ca_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ca_wlm_p(KER3_Ca_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ca_wls_p <- function(KER3_Ca_wls) {
  f <- summary(KER3_Ca_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ca_wls_p(KER3_Ca_wls)

KER3_Ca_wls_eqn <- function(df){
  m <- KER3_Ca_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ca_wls_p(KER3_Ca_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ca_wls_wt_p <- function(KER3_Ca_wls_wt) {
  f <- summary(KER3_Ca_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ca_wls_wt_p(KER3_Ca_wls_wt)

KER3_Ca_wls_wt_eqn <- function(df){
  m <- KER3_Ca_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ca_wls_wt_p(KER3_Ca_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Ca_final <- KER3_Ca + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Ca_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Ca_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Ca_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Ca_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Ca_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ca/Ca_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# KER3_Ti -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Ti_lm <- lm(Ti_ICP ~ Ti, data = ACE_dataset)
summary(KER3_Ti_lm)
glance(KER3_Ti_lm)
model_performance(KER3_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ti_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_OLS_summary.txt")
summary(KER3_Ti_lm)
glance(KER3_Ti_lm)
model_performance(KER3_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ti_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ti_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Ti_wlm <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/(Ti_sd)^2)
summary(KER3_Ti_wlm)
glance(KER3_Ti_wlm)
model_performance(KER3_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ti_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_OLS_wt_summary.txt")
summary(KER3_Ti_wlm)
glance(KER3_Ti_wlm)
model_performance(KER3_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ti_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ti_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Ti_model <- lm(Ti_ICP ~ Ti, data = ACE_dataset) # define model
KER3_Ti_wt <- 1 / lm(abs(KER3_Ti_model$residuals) ~ KER3_Ti_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Ti_wls <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=KER3_Ti_wt) #perform weighted least squares regression
# Checks
summary(KER3_Ti_wls) # summary stats
glance(KER3_Ti_wls) # summary stats including AIC
model_performance(KER3_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ti_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_WLS_summary.txt")
summary(KER3_Ti_wls)
glance(KER3_Ti_wls)
model_performance(KER3_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ti_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ti_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Ti_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Ti_model_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/Ti_sd^2) # define model
KER3_Ti_wt_wt <- 1 / lm(abs(KER3_Ti_model_wt$residuals) ~ KER3_Ti_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Ti_wls_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=KER3_Ti_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Ti_wls_wt) # summary stats
glance(KER3_Ti_wls_wt) # summary stats including AIC
model_performance(KER3_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ti_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_WLS_wt_summary.txt")
summary(KER3_Ti_wls_wt)
glance(KER3_Ti_wls_wt)
model_performance(KER3_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ti_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ti_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Ti_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Ti"
theme_set(theme_bw(10))
KER3_Ti <- ggplot(ACE_dataset, aes(x = Ti, y = Ti_ICP)) +
  geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Ti_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Ti_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Ti_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ti_lm_p <- function(KER3_Ti_lm) {
  f <- summary(KER3_Ti_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ti_lm_p(KER3_Ti_lm)

KER3_Ti_lm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ti_lm_p(KER3_Ti_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ti_wlm_p <- function(KER3_Ti_wlm) {
  f <- summary(KER3_Ti_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ti_wlm_p(KER3_Ti_wlm)

KER3_Ti_wlm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/(Ti_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ti_wlm_p(KER3_Ti_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ti_wls_p <- function(KER3_Ti_wls) {
  f <- summary(KER3_Ti_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ti_wls_p(KER3_Ti_wls)

KER3_Ti_wls_eqn <- function(df){
  m <- KER3_Ti_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ti_wls_p(KER3_Ti_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ti_wls_wt_p <- function(KER3_Ti_wls_wt) {
  f <- summary(KER3_Ti_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ti_wls_wt_p(KER3_Ti_wls_wt)

KER3_Ti_wls_wt_eqn <- function(df){
  m <- KER3_Ti_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ti_wls_wt_p(KER3_Ti_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Ti_final <- KER3_Ti + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Ti_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Ti_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Ti_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Ti_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Ti_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ti/Ti_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# KER3_Mn -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Mn_lm <- lm(Mn_ICP ~ Mn, data = ACE_dataset)
summary(KER3_Mn_lm)
glance(KER3_Mn_lm)
model_performance(KER3_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Mn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_OLS_summary.txt")
summary(KER3_Mn_lm)
glance(KER3_Mn_lm)
model_performance(KER3_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Mn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Mn_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Mn_wlm <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/(Mn_sd)^2)
summary(KER3_Mn_wlm)
glance(KER3_Mn_wlm)
model_performance(KER3_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Mn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_OLS_wt_summary.txt")
summary(KER3_Mn_wlm)
glance(KER3_Mn_wlm)
model_performance(KER3_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Mn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Mn_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Mn_model <- lm(Mn_ICP ~ Mn, data = ACE_dataset) # define model
KER3_Mn_wt <- 1 / lm(abs(KER3_Mn_model$residuals) ~ KER3_Mn_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Mn_wls <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=KER3_Mn_wt) #perform weighted least squares regression
# Checks
summary(KER3_Mn_wls) # summary stats
glance(KER3_Mn_wls) # summary stats including AIC
model_performance(KER3_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Mn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_WLS_summary.txt")
summary(KER3_Mn_wls)
glance(KER3_Mn_wls)
model_performance(KER3_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Mn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Mn_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Mn_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Mn_model_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/Mn_sd^2) # define model
KER3_Mn_wt_wt <- 1 / lm(abs(KER3_Mn_model_wt$residuals) ~ KER3_Mn_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Mn_wls_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=KER3_Mn_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Mn_wls_wt) # summary stats
glance(KER3_Mn_wls_wt) # summary stats including AIC
model_performance(KER3_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Mn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_WLS_wt_summary.txt")
summary(KER3_Mn_wls_wt)
glance(KER3_Mn_wls_wt)
model_performance(KER3_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Mn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Mn_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Mn_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Mn"
theme_set(theme_bw(10))
KER3_Mn <- ggplot(ACE_dataset, aes(x = Mn, y = Mn_ICP)) +
  geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Mn_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Mn_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Mn_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Mn_lm_p <- function(KER3_Mn_lm) {
  f <- summary(KER3_Mn_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Mn_lm_p(KER3_Mn_lm)

KER3_Mn_lm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Mn_lm_p(KER3_Mn_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Mn_wlm_p <- function(KER3_Mn_wlm) {
  f <- summary(KER3_Mn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Mn_wlm_p(KER3_Mn_wlm)

KER3_Mn_wlm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/(Mn_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Mn_wlm_p(KER3_Mn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Mn_wls_p <- function(KER3_Mn_wls) {
  f <- summary(KER3_Mn_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Mn_wls_p(KER3_Mn_wls)

KER3_Mn_wls_eqn <- function(df){
  m <- KER3_Mn_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Mn_wls_p(KER3_Mn_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Mn_wls_wt_p <- function(KER3_Mn_wls_wt) {
  f <- summary(KER3_Mn_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Mn_wls_wt_p(KER3_Mn_wls_wt)

KER3_Mn_wls_wt_eqn <- function(df){
  m <- KER3_Mn_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Mn_wls_wt_p(KER3_Mn_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Mn_final <- KER3_Mn + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Mn_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Mn_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Mn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Mn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Mn_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Mn/Mn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# KER3_Fe -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Fe_lm <- lm(Fe_ICP ~ Fe, data = ACE_dataset)
summary(KER3_Fe_lm)
glance(KER3_Fe_lm)
model_performance(KER3_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Fe_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_OLS_summary.txt")
summary(KER3_Fe_lm)
glance(KER3_Fe_lm)
model_performance(KER3_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Fe_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Fe_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Fe_wlm <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/(Fe_sd)^2)
summary(KER3_Fe_wlm)
glance(KER3_Fe_wlm)
model_performance(KER3_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Fe_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_OLS_wt_summary.txt")
summary(KER3_Fe_wlm)
glance(KER3_Fe_wlm)
model_performance(KER3_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Fe_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Fe_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Fe_model <- lm(Fe_ICP ~ Fe, data = ACE_dataset) # define model
KER3_Fe_wt <- 1 / lm(abs(KER3_Fe_model$residuals) ~ KER3_Fe_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Fe_wls <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=KER3_Fe_wt) #perform weighted least squares regression
# Checks
summary(KER3_Fe_wls) # summary stats
glance(KER3_Fe_wls) # summary stats including AIC
model_performance(KER3_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Fe_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_WLS_summary.txt")
summary(KER3_Fe_wls)
glance(KER3_Fe_wls)
model_performance(KER3_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Fe_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Fe_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Fe_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Fe_model_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/Fe_sd^2) # define model
KER3_Fe_wt_wt <- 1 / lm(abs(KER3_Fe_model_wt$residuals) ~ KER3_Fe_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Fe_wls_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=KER3_Fe_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Fe_wls_wt) # summary stats
glance(KER3_Fe_wls_wt) # summary stats including AIC
model_performance(KER3_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Fe_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_WLS_wt_summary.txt")
summary(KER3_Fe_wls_wt)
glance(KER3_Fe_wls_wt)
model_performance(KER3_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Fe_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Fe_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Fe_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Fe"
theme_set(theme_bw(10))
KER3_Fe <- ggplot(ACE_dataset, aes(x = Fe, y = Fe_ICP)) +
  geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Fe_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Fe_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Fe_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Fe_lm_p <- function(KER3_Fe_lm) {
  f <- summary(KER3_Fe_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Fe_lm_p(KER3_Fe_lm)

KER3_Fe_lm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Fe_lm_p(KER3_Fe_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Fe_wlm_p <- function(KER3_Fe_wlm) {
  f <- summary(KER3_Fe_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Fe_wlm_p(KER3_Fe_wlm)

KER3_Fe_wlm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/(Fe_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Fe_wlm_p(KER3_Fe_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Fe_wls_p <- function(KER3_Fe_wls) {
  f <- summary(KER3_Fe_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Fe_wls_p(KER3_Fe_wls)

KER3_Fe_wls_eqn <- function(df){
  m <- KER3_Fe_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Fe_wls_p(KER3_Fe_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Fe_wls_wt_p <- function(KER3_Fe_wls_wt) {
  f <- summary(KER3_Fe_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Fe_wls_wt_p(KER3_Fe_wls_wt)

KER3_Fe_wls_wt_eqn <- function(df){
  m <- KER3_Fe_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Fe_wls_wt_p(KER3_Fe_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Fe_final <- KER3_Fe + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Fe_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Fe_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Fe_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Fe_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3) %>% 
  print()
KER3_Fe_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Fe/Fe_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# KER3_Co -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Co_lm <- lm(Co_ICP ~ Co, data = ACE_dataset)
summary(KER3_Co_lm)
glance(KER3_Co_lm)
model_performance(KER3_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Co_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_OLS_summary.txt")
summary(KER3_Co_lm)
glance(KER3_Co_lm)
model_performance(KER3_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Co_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Co_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Co_wlm <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/(Co_sd)^2)
summary(KER3_Co_wlm)
glance(KER3_Co_wlm)
model_performance(KER3_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Co_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_OLS_wt_summary.txt")
summary(KER3_Co_wlm)
glance(KER3_Co_wlm)
model_performance(KER3_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Co_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Co_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Co_model <- lm(Co_ICP ~ Co, data = ACE_dataset) # define model
KER3_Co_wt <- 1 / lm(abs(KER3_Co_model$residuals) ~ KER3_Co_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Co_wls <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=KER3_Co_wt) #perform weighted least squares regression
# Checks
summary(KER3_Co_wls) # summary stats
glance(KER3_Co_wls) # summary stats including AIC
model_performance(KER3_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Co_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_WLS_summary.txt")
summary(KER3_Co_wls)
glance(KER3_Co_wls)
model_performance(KER3_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Co_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Co_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Co_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Co_model_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/Co_sd^2) # define model
KER3_Co_wt_wt <- 1 / lm(abs(KER3_Co_model_wt$residuals) ~ KER3_Co_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Co_wls_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=KER3_Co_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Co_wls_wt) # summary stats
glance(KER3_Co_wls_wt) # summary stats including AIC
model_performance(KER3_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Co_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_WLS_wt_summary.txt")
summary(KER3_Co_wls_wt)
glance(KER3_Co_wls_wt)
model_performance(KER3_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Co_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Co_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Co_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Co"
theme_set(theme_bw(10))
KER3_Co <- ggplot(ACE_dataset, aes(x = Co, y = Co_ICP)) +
  geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Co_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Co_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Co_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Co_lm_p <- function(KER3_Co_lm) {
  f <- summary(KER3_Co_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Co_lm_p(KER3_Co_lm)

KER3_Co_lm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Co_lm_p(KER3_Co_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Co_wlm_p <- function(KER3_Co_wlm) {
  f <- summary(KER3_Co_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Co_wlm_p(KER3_Co_wlm)

KER3_Co_wlm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/(Co_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Co_wlm_p(KER3_Co_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Co_wls_p <- function(KER3_Co_wls) {
  f <- summary(KER3_Co_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Co_wls_p(KER3_Co_wls)

KER3_Co_wls_eqn <- function(df){
  m <- KER3_Co_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Co_wls_p(KER3_Co_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Co_wls_wt_p <- function(KER3_Co_wls_wt) {
  f <- summary(KER3_Co_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Co_wls_wt_p(KER3_Co_wls_wt)

KER3_Co_wls_wt_eqn <- function(df){
  m <- KER3_Co_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Co_wls_wt_p(KER3_Co_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Co_final <- KER3_Co + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Co_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Co_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Co_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Co_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Co_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Co/Co_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# KER3_Ni -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Ni_lm <- lm(Ni_ICP ~ Ni, data = ACE_dataset)
summary(KER3_Ni_lm)
glance(KER3_Ni_lm)
model_performance(KER3_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ni_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_OLS_summary.txt")
summary(KER3_Ni_lm)
glance(KER3_Ni_lm)
model_performance(KER3_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ni_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ni_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Ni_wlm <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/(Ni_sd)^2)
summary(KER3_Ni_wlm)
glance(KER3_Ni_wlm)
model_performance(KER3_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ni_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_OLS_wt_summary.txt")
summary(KER3_Ni_wlm)
glance(KER3_Ni_wlm)
model_performance(KER3_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ni_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ni_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Ni_model <- lm(Ni_ICP ~ Ni, data = ACE_dataset) # define model
KER3_Ni_wt <- 1 / lm(abs(KER3_Ni_model$residuals) ~ KER3_Ni_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Ni_wls <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=KER3_Ni_wt) #perform weighted least squares regression
# Checks
summary(KER3_Ni_wls) # summary stats
glance(KER3_Ni_wls) # summary stats including AIC
model_performance(KER3_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ni_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_WLS_summary.txt")
summary(KER3_Ni_wls)
glance(KER3_Ni_wls)
model_performance(KER3_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ni_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ni_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Ni_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Ni_model_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/Ni_sd^2) # define model
KER3_Ni_wt_wt <- 1 / lm(abs(KER3_Ni_model_wt$residuals) ~ KER3_Ni_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Ni_wls_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=KER3_Ni_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Ni_wls_wt) # summary stats
glance(KER3_Ni_wls_wt) # summary stats including AIC
model_performance(KER3_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Ni_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_WLS_wt_summary.txt")
summary(KER3_Ni_wls_wt)
glance(KER3_Ni_wls_wt)
model_performance(KER3_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Ni_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ni_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Ni_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Ni"
theme_set(theme_bw(10))
KER3_Ni <- ggplot(ACE_dataset, aes(x = Ni, y = Ni_ICP)) +
  geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Ni_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Ni_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Ni_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ni_lm_p <- function(KER3_Ni_lm) {
  f <- summary(KER3_Ni_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ni_lm_p(KER3_Ni_lm)

KER3_Ni_lm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ni_lm_p(KER3_Ni_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ni_wlm_p <- function(KER3_Ni_wlm) {
  f <- summary(KER3_Ni_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ni_wlm_p(KER3_Ni_wlm)

KER3_Ni_wlm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/(Ni_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ni_wlm_p(KER3_Ni_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ni_wls_p <- function(KER3_Ni_wls) {
  f <- summary(KER3_Ni_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ni_wls_p(KER3_Ni_wls)

KER3_Ni_wls_eqn <- function(df){
  m <- KER3_Ni_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ni_wls_p(KER3_Ni_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Ni_wls_wt_p <- function(KER3_Ni_wls_wt) {
  f <- summary(KER3_Ni_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Ni_wls_wt_p(KER3_Ni_wls_wt)

KER3_Ni_wls_wt_eqn <- function(df){
  m <- KER3_Ni_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Ni_wls_wt_p(KER3_Ni_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Ni_final <- KER3_Ni + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Ni_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Ni_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Ni_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Ni_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Ni_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Ni/Ni_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# KER3_Cu -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Cu_lm <- lm(Cu_ICP ~ Cu, data = ACE_dataset)
summary(KER3_Cu_lm)
glance(KER3_Cu_lm)
model_performance(KER3_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Cu_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_OLS_summary.txt")
summary(KER3_Cu_lm)
glance(KER3_Cu_lm)
model_performance(KER3_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Cu_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Cu_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Cu_lm) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Cu_wlm <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/(Cu_sd)^2)
summary(KER3_Cu_wlm)
glance(KER3_Cu_wlm)
model_performance(KER3_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Cu_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_OLS_wt_summary.txt")
summary(KER3_Cu_wlm)
glance(KER3_Cu_wlm)
model_performance(KER3_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Cu_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Cu_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Cu_lm) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Cu_model <- lm(Cu_ICP ~ Cu, data = ACE_dataset) # define model
KER3_Cu_wt <- 1 / lm(abs(KER3_Cu_model$residuals) ~ KER3_Cu_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Cu_wls <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=KER3_Cu_wt) #perform weighted least squares regression
# Checks
summary(KER3_Cu_wls) # summary stats
glance(KER3_Cu_wls) # summary stats including AIC
model_performance(KER3_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Cu_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_WLS_summary.txt")
summary(KER3_Cu_wls)
glance(KER3_Cu_wls)
model_performance(KER3_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Cu_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Cu_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Cu_wls) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Cu_model_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/Cu_sd^2) # define model
KER3_Cu_wt_wt <- 1 / lm(abs(KER3_Cu_model_wt$residuals) ~ KER3_Cu_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Cu_wls_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=KER3_Cu_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Cu_wls_wt) # summary stats
glance(KER3_Cu_wls_wt) # summary stats including AIC
model_performance(KER3_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Cu_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_WLS_wt_summary.txt")
summary(KER3_Cu_wls_wt)
glance(KER3_Cu_wls_wt)
model_performance(KER3_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Cu_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Cu_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Cu_wls_wt) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Cu"
theme_set(theme_bw(10))
KER3_Cu <- ggplot(ACE_dataset, aes(x = Cu, y = Cu_ICP)) +
  geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Cu_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Cu_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Cu_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Cu_lm_p <- function(KER3_Cu_lm) {
  f <- summary(KER3_Cu_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Cu_lm_p(KER3_Cu_lm)

KER3_Cu_lm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Cu_lm_p(KER3_Cu_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Cu_wlm_p <- function(KER3_Cu_wlm) {
  f <- summary(KER3_Cu_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Cu_wlm_p(KER3_Cu_wlm)

KER3_Cu_wlm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/(Cu_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Cu_wlm_p(KER3_Cu_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Cu_wls_p <- function(KER3_Cu_wls) {
  f <- summary(KER3_Cu_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Cu_wls_p(KER3_Cu_wls)

KER3_Cu_wls_eqn <- function(df){
  m <- KER3_Cu_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Cu_wls_p(KER3_Cu_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Cu_wls_wt_p <- function(KER3_Cu_wls_wt) {
  f <- summary(KER3_Cu_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Cu_wls_wt_p(KER3_Cu_wls_wt)

KER3_Cu_wls_wt_eqn <- function(df){
  m <- KER3_Cu_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Cu_wls_wt_p(KER3_Cu_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Cu_final <- KER3_Cu + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Cu_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Cu_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Cu_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Cu_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Cu_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Cu/Cu_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# KER3_Zn -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Zn_lm <- lm(Zn_ICP ~ Zn, data = ACE_dataset)
summary(KER3_Zn_lm)
glance(KER3_Zn_lm)
model_performance(KER3_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_OLS_summary.txt")
summary(KER3_Zn_lm)
glance(KER3_Zn_lm)
model_performance(KER3_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zn_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Zn_lm) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Zn_wlm <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/(Zn_sd)^2)
summary(KER3_Zn_wlm)
glance(KER3_Zn_wlm)
model_performance(KER3_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_OLS_wt_summary.txt")
summary(KER3_Zn_wlm)
glance(KER3_Zn_wlm)
model_performance(KER3_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zn_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Zn_lm) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Zn_model <- lm(Zn_ICP ~ Zn, data = ACE_dataset) # define model
KER3_Zn_wt <- 1 / lm(abs(KER3_Zn_model$residuals) ~ KER3_Zn_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Zn_wls <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=KER3_Zn_wt) #perform weighted least squares regression
# Checks
summary(KER3_Zn_wls) # summary stats
glance(KER3_Zn_wls) # summary stats including AIC
model_performance(KER3_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_WLS_summary.txt")
summary(KER3_Zn_wls)
glance(KER3_Zn_wls)
model_performance(KER3_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zn_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Zn_wls) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Zn_model_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/Zn_sd^2) # define model
KER3_Zn_wt_wt <- 1 / lm(abs(KER3_Zn_model_wt$residuals) ~ KER3_Zn_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Zn_wls_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=KER3_Zn_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Zn_wls_wt) # summary stats
glance(KER3_Zn_wls_wt) # summary stats including AIC
model_performance(KER3_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_WLS_wt_summary.txt")
summary(KER3_Zn_wls_wt)
glance(KER3_Zn_wls_wt)
model_performance(KER3_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zn_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Zn_wls_wt) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Zn"
theme_set(theme_bw(10))
KER3_Zn <- ggplot(ACE_dataset, aes(x = Zn, y = Zn_ICP)) +
  geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Zn_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Zn_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Zn_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zn_lm_p <- function(KER3_Zn_lm) {
  f <- summary(KER3_Zn_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zn_lm_p(KER3_Zn_lm)

KER3_Zn_lm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zn_lm_p(KER3_Zn_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zn_wlm_p <- function(KER3_Zn_wlm) {
  f <- summary(KER3_Zn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zn_wlm_p(KER3_Zn_wlm)

KER3_Zn_wlm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/(Zn_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zn_wlm_p(KER3_Zn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zn_wls_p <- function(KER3_Zn_wls) {
  f <- summary(KER3_Zn_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zn_wls_p(KER3_Zn_wls)

KER3_Zn_wls_eqn <- function(df){
  m <- KER3_Zn_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zn_wls_p(KER3_Zn_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zn_wls_wt_p <- function(KER3_Zn_wls_wt) {
  f <- summary(KER3_Zn_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zn_wls_wt_p(KER3_Zn_wls_wt)

KER3_Zn_wls_wt_eqn <- function(df){
  m <- KER3_Zn_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zn_wls_wt_p(KER3_Zn_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Zn_final <- KER3_Zn + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Zn_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Zn_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Zn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Zn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Zn_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zn/Zn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# KER3_Rb -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Rb_lm <- lm(Rb_ICP ~ Rb, data = ACE_dataset)
summary(KER3_Rb_lm)
glance(KER3_Rb_lm)
model_performance(KER3_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Rb_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_OLS_summary.txt")
summary(KER3_Rb_lm)
glance(KER3_Rb_lm)
model_performance(KER3_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Rb_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Rb_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Rb_lm) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Rb_wlm <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/(Rb_sd)^2)
summary(KER3_Rb_wlm)
glance(KER3_Rb_wlm)
model_performance(KER3_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Rb_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_OLS_wt_summary.txt")
summary(KER3_Rb_wlm)
glance(KER3_Rb_wlm)
model_performance(KER3_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Rb_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Rb_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Rb_lm) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Rb_model <- lm(Rb_ICP ~ Rb, data = ACE_dataset) # define model
KER3_Rb_wt <- 1 / lm(abs(KER3_Rb_model$residuals) ~ KER3_Rb_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Rb_wls <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=KER3_Rb_wt) #perform weighted least squares regression
# Checks
summary(KER3_Rb_wls) # summary stats
glance(KER3_Rb_wls) # summary stats including AIC
model_performance(KER3_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Rb_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_WLS_summary.txt")
summary(KER3_Rb_wls)
glance(KER3_Rb_wls)
model_performance(KER3_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Rb_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Rb_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Rb_wls) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Rb_model_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/Rb_sd^2) # define model
KER3_Rb_wt_wt <- 1 / lm(abs(KER3_Rb_model_wt$residuals) ~ KER3_Rb_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Rb_wls_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=KER3_Rb_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Rb_wls_wt) # summary stats
glance(KER3_Rb_wls_wt) # summary stats including AIC
model_performance(KER3_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Rb_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_WLS_wt_summary.txt")
summary(KER3_Rb_wls_wt)
glance(KER3_Rb_wls_wt)
model_performance(KER3_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Rb_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Rb_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Rb_wls_wt) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Rb"
theme_set(theme_bw(10))
KER3_Rb <- ggplot(ACE_dataset, aes(x = Rb, y = Rb_ICP)) +
  geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Rb_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Rb_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Rb_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Rb_lm_p <- function(KER3_Rb_lm) {
  f <- summary(KER3_Rb_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Rb_lm_p(KER3_Rb_lm)

KER3_Rb_lm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Rb_lm_p(KER3_Rb_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Rb_wlm_p <- function(KER3_Rb_wlm) {
  f <- summary(KER3_Rb_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Rb_wlm_p(KER3_Rb_wlm)

KER3_Rb_wlm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/(Rb_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Rb_wlm_p(KER3_Rb_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Rb_wls_p <- function(KER3_Rb_wls) {
  f <- summary(KER3_Rb_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Rb_wls_p(KER3_Rb_wls)

KER3_Rb_wls_eqn <- function(df){
  m <- KER3_Rb_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Rb_wls_p(KER3_Rb_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Rb_wls_wt_p <- function(KER3_Rb_wls_wt) {
  f <- summary(KER3_Rb_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Rb_wls_wt_p(KER3_Rb_wls_wt)

KER3_Rb_wls_wt_eqn <- function(df){
  m <- KER3_Rb_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Rb_wls_wt_p(KER3_Rb_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Rb_final <- KER3_Rb + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Rb_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Rb_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Rb_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Rb_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Rb_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Rb/Rb_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# KER3_Sr -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Sr_lm <- lm(Sr_ICP ~ Sr, data = ACE_dataset)
summary(KER3_Sr_lm)
glance(KER3_Sr_lm)
model_performance(KER3_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Sr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_OLS_summary.txt")
summary(KER3_Sr_lm)
glance(KER3_Sr_lm)
model_performance(KER3_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Sr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Sr_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Sr_lm) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Sr_wlm <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/(Sr_sd)^2)
summary(KER3_Sr_wlm)
glance(KER3_Sr_wlm)
model_performance(KER3_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Sr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_OLS_wt_summary.txt")
summary(KER3_Sr_wlm)
glance(KER3_Sr_wlm)
model_performance(KER3_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Sr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Sr_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Sr_lm) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Sr_model <- lm(Sr_ICP ~ Sr, data = ACE_dataset) # define model
KER3_Sr_wt <- 1 / lm(abs(KER3_Sr_model$residuals) ~ KER3_Sr_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Sr_wls <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=KER3_Sr_wt) #perform weighted least squares regression
# Checks
summary(KER3_Sr_wls) # summary stats
glance(KER3_Sr_wls) # summary stats including AIC
model_performance(KER3_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Sr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_WLS_summary.txt")
summary(KER3_Sr_wls)
glance(KER3_Sr_wls)
model_performance(KER3_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Sr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Sr_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Sr_wls) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Sr_model_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/Sr_sd^2) # define model
KER3_Sr_wt_wt <- 1 / lm(abs(KER3_Sr_model_wt$residuals) ~ KER3_Sr_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Sr_wls_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=KER3_Sr_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Sr_wls_wt) # summary stats
glance(KER3_Sr_wls_wt) # summary stats including AIC
model_performance(KER3_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Sr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_WLS_wt_summary.txt")
summary(KER3_Sr_wls_wt)
glance(KER3_Sr_wls_wt)
model_performance(KER3_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Sr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Sr_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Sr_wls_wt) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Sr"
theme_set(theme_bw(10))
KER3_Sr <- ggplot(ACE_dataset, aes(x = Sr, y = Sr_ICP)) +
  geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Sr_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Sr_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Sr_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Sr_lm_p <- function(KER3_Sr_lm) {
  f <- summary(KER3_Sr_lm_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Sr_lm_p(KER3_Sr_lm_wt)

KER3_Sr_lm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Sr_lm_p(KER3_Sr_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Sr_wlm_p <- function(KER3_Sr_wlm) {
  f <- summary(KER3_Sr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Sr_wlm_p(KER3_Sr_wlm)

KER3_Sr_wlm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/(Sr_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Sr_wlm_p(KER3_Sr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Sr_wls_p <- function(KER3_Sr_wls) {
  f <- summary(KER3_Sr_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Sr_wls_p(KER3_Sr_wls)

KER3_Sr_wls_eqn <- function(df){
  m <- KER3_Sr_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Sr_wls_p(KER3_Sr_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Sr_wls_wt_p <- function(KER3_Sr_wls_wt) {
  f <- summary(KER3_Sr_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Sr_wls_wt_p(KER3_Sr_wls)

KER3_Sr_wls_wt_eqn <- function(df){
  m <- KER3_Sr_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Sr_wls_wt_p(KER3_Sr_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Sr_final <- KER3_Sr + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Sr_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Sr_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Sr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Sr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Sr_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Sr/Sr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# KER3_Zr -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_Zr_lm <- lm(Zr_ICP ~ Zr, data = ACE_dataset)
summary(KER3_Zr_lm)
glance(KER3_Zr_lm)
model_performance(KER3_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_OLS_summary.txt")
summary(KER3_Zr_lm)
glance(KER3_Zr_lm)
model_performance(KER3_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zr_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Zr_lm) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_Zr_wlm <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/(Zr_sd)^2)
summary(KER3_Zr_wlm)
glance(KER3_Zr_wlm)
model_performance(KER3_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_OLS_wt_summary.txt")
summary(KER3_Zr_wlm)
glance(KER3_Zr_wlm)
model_performance(KER3_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zr_lm) # Performance package summary check for heteroscedasticity
icc(KER3_Zr_lm) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_Zr_model <- lm(Zr_ICP ~ Zr, data = ACE_dataset) # define model
KER3_Zr_wt <- 1 / lm(abs(KER3_Zr_model$residuals) ~ KER3_Zr_model$fitted.values)$fitted.values^2 #define weights to use
KER3_Zr_wls <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=KER3_Zr_wt) #perform weighted least squares regression
# Checks
summary(KER3_Zr_wls) # summary stats
glance(KER3_Zr_wls) # summary stats including AIC
model_performance(KER3_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_WLS_summary.txt")
summary(KER3_Zr_wls)
glance(KER3_Zr_wls)
model_performance(KER3_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zr_wls) # Performance package summary check for heteroscedasticity
icc(KER3_Zr_wls) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_Zr_model_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/Zr_sd^2) # define model
KER3_Zr_wt_wt <- 1 / lm(abs(KER3_Zr_model_wt$residuals) ~ KER3_Zr_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_Zr_wls_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=KER3_Zr_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_Zr_wls_wt) # summary stats
glance(KER3_Zr_wls_wt) # summary stats including AIC
model_performance(KER3_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_Zr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_WLS_wt_summary.txt")
summary(KER3_Zr_wls_wt)
glance(KER3_Zr_wls_wt)
model_performance(KER3_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_Zr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zr_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_Zr_wls_wt) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Zr"
theme_set(theme_bw(10))
KER3_Zr <- ggplot(ACE_dataset, aes(x = Zr, y = Zr_ICP)) +
  geom_errorbar(aes(ymin=Zr_ICP-Zr_ICP_sd, ymax=Zr_ICP+Zr_ICP_sd), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=Zr-Zr_sd, xmax=Zr+Zr_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Zr_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_Zr_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_Zr_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste(element_title, " [XRF-CS]", dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, dataset))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zr_lm_p <- function(KER3_Zr_lm) {
  f <- summary(KER3_Zr_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zr_lm_p(KER3_Zr_lm)

KER3_Zr_lm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zr_lm_p(KER3_Zr_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zr_wlm_p <- function(KER3_Zr_wlm) {
  f <- summary(KER3_Zr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zr_wlm_p(KER3_Zr_wlm)

KER3_Zr_wlm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/(Zr_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zr_wlm_p(KER3_Zr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zr_wls_p <- function(KER3_Zr_wls) {
  f <- summary(KER3_Zr_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zr_wls_p(KER3_Zr_wls)

KER3_Zr_wls_eqn <- function(df){
  m <- KER3_Zr_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zr_wls_p(KER3_Zr_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_Zr_wls_wt_p <- function(KER3_Zr_wls_wt) {
  f <- summary(KER3_Zr_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_Zr_wls_wt_p(KER3_Zr_wls_wt)

KER3_Zr_wls_wt_eqn <- function(df){
  m <- KER3_Zr_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_Zr_wls_wt_p(KER3_Zr_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_Zr_final <- KER3_Zr + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_Zr_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_Zr_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_Zr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_Zr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_Zr_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/Zr/Zr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# KER3_DM -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
KER3_DM_lm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset)
summary(KER3_DM_lm)
glance(KER3_DM_lm)
model_performance(KER3_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_OLS_summary.txt")
summary(KER3_DM_lm)
glance(KER3_DM_lm)
model_performance(KER3_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_lm) # Performance package summary check for heteroscedasticity
icc(KER3_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
KER3_DM_wlm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/(coh_inc_sd)^2)
summary(KER3_DM_wlm)
glance(KER3_DM_wlm)
model_performance(KER3_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_DM_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_OLS_wt_summary.txt")
summary(KER3_DM_wlm)
glance(KER3_DM_wlm)
model_performance(KER3_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_DM_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_lm) # Performance package summary check for heteroscedasticity
icc(KER3_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
KER3_DM_model <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset) # define model
KER3_DM_wt <- 1 / lm(abs(KER3_DM_model$residuals) ~ KER3_DM_model$fitted.values)$fitted.values^2 #define weights to use
KER3_DM_wls <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=KER3_DM_wt) #perform weighted least squares regression
# Checks
summary(KER3_DM_wls) # summary stats
glance(KER3_DM_wls) # summary stats including AIC
model_performance(KER3_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_DM_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_WLS_summary.txt")
summary(KER3_DM_wls)
glance(KER3_DM_wls)
model_performance(KER3_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_DM_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_wls) # Performance package summary check for heteroscedasticity
icc(KER3_DM_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
KER3_DM_model_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/coh_inc_sd^2) # define model
KER3_DM_wt_wt <- 1 / lm(abs(KER3_DM_model_wt$residuals) ~ KER3_DM_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_DM_wls_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=KER3_DM_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_DM_wls_wt) # summary stats
glance(KER3_DM_wls_wt) # summary stats including AIC
model_performance(KER3_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_DM_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_WLS_wt_summary.txt")
summary(KER3_DM_wls_wt)
glance(KER3_DM_wls_wt)
model_performance(KER3_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_DM_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_DM_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots
element_title <- "Dry mass, coh/inc"
theme_set(theme_bw(10))
KER3_DM <- ggplot(ACE_dataset, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=0, colour = "grey") +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(colour = "OLS")) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/coh_inc_sd^2, colour = "OLS_wt")) +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = 'top', label.x = 'right') + 
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = KER3_DM_wt, colour="WLS")) + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = KER3_DM_wt_wt, colour="WLS_wt")) + # weighted WLS regression
  geom_point(size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_colour_manual(name="legend", values=c("#FF9999", "red", "#6699FF", "blue")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.position = "bottom", 
        #legend.justification = c("top", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = paste0("coh/inc (XRF-CS) / ", dataset), y = paste0("Dry mass (%)")) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))

# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_DM_lm_p <- function(KER3_DM_lm) {
  f <- summary(KER3_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_DM_lm_p(KER3_DM_lm)

KER3_DM_lm_eqn <- function(df){
  m <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_DM_lm_p(KER3_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_DM_wlm_p <- function(KER3_DM_wlm) {
  f <- summary(KER3_DM_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_DM_wlm_p(KER3_DM_wlm)

KER3_DM_wlm_eqn <- function(df){
  m <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/(coh_inc_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_DM_wlm_p(KER3_DM_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_DM_wls_p <- function(KER3_DM_wls) {
  f <- summary(KER3_DM_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_DM_wls_p(KER3_DM_wls)

KER3_DM_wls_eqn <- function(df){
  m <- KER3_DM_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_DM_wls_p(KER3_DM_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
KER3_DM_wls_wt_p <- function(KER3_DM_wls_wt) {
  f <- summary(KER3_DM_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_DM_wls_wt_p(KER3_DM_wls_wt)

KER3_DM_wls_wt_eqn <- function(df){
  m <- KER3_DM_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_DM_wls_wt_p(KER3_DM_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_DM_final <- KER3_DM + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_DM_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_DM_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#6699FF", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_DM_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "red", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "#FF9999", hjust = -0.1, vjust = 5, size = 3)
KER3_DM_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/DM/DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Summary 4x3 matrix plots of ITRAX-ACF matched elements --------

# OLS & WLS summary - unweighted stats on plot
ggarrange(KER3_K_final, KER3_Ca_final, KER3_Ti_final,
          KER3_Mn_final, KER3_Fe_final, KER3_Co_final,
          KER3_Ni_final, KER3_Cu_final, KER3_Zn_final,
          KER3_Rb_final, KER3_Sr_final, KER3_Zr_final,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/KER3_OLS_WLS_Summary.pdf",
       height = c(50), width = c(50), dpi = 600, units = "cm")














# -------------------------------------------------------------------------
# OLS Linear models - forced through origin
# K unweighted (black) forced through origin ---------------------------
KER3_K_origin <- ggplot(KER3_xrf_icp_matched, aes(x = K, y = K_ICP)) +
  geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=.1, colour = "grey") +
  geom_smooth(method=glm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/K_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "K (XRF-CS) / cps", y = "K (ICPMS) / ppm" ) +
  ggtitle("KER3: K (ICP SD weighted = blue; unweighted = black)")
KER3_K_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/K_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# Ca unweighted (black) forced through origin ---------------------------
KER3_Ca_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Ca, y = Ca_ICP)) +
  geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Ca_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Ca (XRF-CS) / cps", y = "Ca (ICPMS) / ppm" ) +
  ggtitle("KER3: Ca (ICP SD weighted = blue; unweighted = black)")
KER3_Ca_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ca_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Ti 

# Ti unweighted (black) forced through origin ---------------------------
KER3_Ti_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Ti, y = Ti_ICP)) +
  geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Ti_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Ti (XRF-CS) / cps", y = "Ti (ICPMS) / ppm" ) +
  ggtitle("KER3: Ti (ICP SD weighted = blue; unweighted = black)")
KER3_Ti_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ti_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Mn

# Mn unweighted (black) forced through origin ---------------------------
KER3_Mn_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Mn, y = Mn_ICP)) +
  geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Mn_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Mn (XRF-CS) / cps", y = "Mn (ICPMS) / ppm" ) +
  ggtitle("KER3: Mn (ICP SD weighted = blue; unweighted = black)")
KER3_Mn_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/incorigin//Mn_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Fe unweighted (black) forced through origin ---------------------------
KER3_Fe_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Fe, y = Fe_ICP)) +
  geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Fe_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Fe (XRF-CS) / cps", y = "Fe (ICPMS) / ppm" ) +
  ggtitle("KER3: Fe (ICP SD weighted = blue; unweighted = black)")
KER3_Fe_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Fe_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Co

# Co unweighted (black) forced through origin ---------------------------
KER3_Co_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Co, y = Co_ICP)) +
  geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Co_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Co (XRF-CS) / cps", y = "Co (ICPMS) / ppm" ) +
  ggtitle("KER3: Co (ICP SD weighted = blue; unweighted = black)")
KER3_Co_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Co_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Ni

# Ni unweighted (black) forced through origin ---------------------------
KER3_Ni_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Ni, y = Ni_ICP)) +
  geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Ni_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Ni (XRF-CS) / cps", y = "Ni (ICPMS) / ppm" ) +
  ggtitle("KER3: Ni (ICP SD weighted = blue; unweighted = black)")
KER3_Ni_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ni_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Cu

# Cu unweighted (black) forced through origin ---------------------------
KER3_Cu_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Cu, y = Cu_ICP)) +
  geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Cu_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Cu (XRF-CS) / cps", y = "Cu (ICPMS) / ppm" ) +
  ggtitle("KER3: Cu (ICP SD weighted = blue; unweighted = black)")
KER3_Cu_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Cu_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Zn

# Zn unweighted (black) forced through origin ---------------------------
KER3_Zn_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Zn, y = Zn_ICP)) +
  geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Zn_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Zn (XRF-CS) / cps", y = "Zn (ICPMS) / ppm" ) +
  ggtitle("KER3: Zn (ICP SD weighted = blue; unweighted = black)")
KER3_Zn_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Zn_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Rb

# Rb unweighted (black) forced through origin ---------------------------
KER3_Rb_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Rb, y = Rb_ICP)) +
  geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Rb_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Rb (XRF-CS) / cps", y = "Rb (ICPMS) / ppm" ) +
  ggtitle("KER3: Rb (ICP SD weighted = blue; unweighted = black)")
KER3_Rb_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Rb_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Sr

# Sr unweighted (black) forced through origin ---------------------------
KER3_Sr_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Sr, y = Sr_ICP)) +
  geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Sr_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Sr (XRF-CS) / cps", y = "Sr (ICPMS) / ppm" ) +
  ggtitle("KER3: Sr (ICP SD weighted = blue; unweighted = black)")
KER3_Sr_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Sr_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# DM

# DM unweighted (black) forced through origin ---------------------------
KER3_DM_origin <- ggplot(KER3_xrf_icp_matched, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x-1, aes(weight = 1/dry_mass_err)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "coh/inc (XRF-CS) / cps", y = "DM (ICPMS) / ppm" ) +
  ggtitle("KER3: DM (ICP SD weighted = blue; unweighted = black)")
KER3_DM_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/DM_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")






# Performance tests & stats - forced through origin - lm and wlm ------------------------

# K

#Performance: K unweighted lm stats forced through origin
KER3_K_lm_origin <- lm(K_ICP ~ K-1, data = KER3_xrf_icp_matched)
summary(KER3_K_lm_origin) # summary stats
check_model(KER3_K_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/K_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_K_lm_origin_mp <- model_performance(KER3_K_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_K_lm_origin_mp)
bptest(KER3_K_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_K_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_K_lm)

#Performance: K weighted lm forced through origin
KER3_K_wlm_origin <- lm(K_ICP ~ K-1, weight = 1/K_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_K_wlm_origin) # summary stats

check_model(KER3_K_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/K_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_K_wlm_origin_mp <- model_performance(KER3_K_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_K_wlm_origin_mp)
bptest(KER3_K_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_K_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_K_wlm)

# Ca

#Performance: Ca unweighted lm stats forced through origin
KER3_Ca_lm_origin <- lm(Ca_ICP ~ Ca-1, data = KER3_xrf_icp_matched)
summary(KER3_Ca_lm_origin) # summary stats
check_model(KER3_Ca_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ca_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Ca_lm_origin_mp <- model_performance(KER3_Ca_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Ca_lm_origin_mp)
bptest(KER3_Ca_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ca_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Ca_lm)

#Performance: Ca weighted lm forced through origin
KER3_Ca_wlm_origin <- lm(Ca_ICP ~ Ca-1, weight = 1/Ca_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Ca_wlm_origin) # summary stats
check_model(KER3_Ca_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ca_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Ca_wlm_origin_mp <- model_performance(KER3_Ca_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Ca_wlm_origin_mp)
bptest(KER3_Ca_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ca_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Ca_wlm)

# Ti

#Performance: Ti unweighted lm stats forced through origin
KER3_Ti_lm_origin <- lm(Ti_ICP ~ Ti-1, data = KER3_xrf_icp_matched)
summary(KER3_Ti_lm_origin) # summary stats
check_model(KER3_Ti_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ti_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Ti_lm_origin_mp <- model_performance(KER3_Ti_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Ti_lm_origin_mp)
bptest(KER3_Ti_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ti_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Ti_lm)

#Performance: Ti weighted lm forced through origin
KER3_Ti_wlm_origin <- lm(Ti_ICP ~ Ti-1, weight = 1/Ti_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Ti_wlm_origin) # summary stats
check_model(KER3_Ti_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ti_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Ti_wlm_origin_mp <- model_performance(KER3_Ti_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Ti_wlm_origin_mp)
bptest(KER3_Ti_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ti_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Ti_wlm)

# Mn

#Performance: Mn unweighted lm stats forced through origin
KER3_Mn_lm_origin <- lm(Mn_ICP ~ Mn-1, data = KER3_xrf_icp_matched)
summary(KER3_Mn_lm_origin) # summary stats
check_model(KER3_Mn_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Mn_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Mn_lm_origin_mp <- model_performance(KER3_Mn_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Mn_lm_origin_mp)
bptest(KER3_Mn_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Mn_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Mn_lm)

#Performance: Mn weighted lm forced through origin
KER3_Mn_wlm_origin <- lm(Mn_ICP ~ Mn-1, weight = 1/Mn_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Mn_wlm_origin) # summary stats
check_model(KER3_Mn_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Mn_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Mn_wlm_origin_mp <- model_performance(KER3_Mn_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Mn_wlm_origin_mp)
bptest(KER3_Mn_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Mn_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Mn_wlm)

# Fe

#Performance: Fe unweighted lm stats forced through origin
KER3_Fe_lm_origin <- lm(Fe_ICP ~ Fe-1, data = KER3_xrf_icp_matched)
summary(KER3_Fe_lm_origin) # summary stats
check_model(KER3_Fe_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Fe_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Fe_lm_origin_mp <- model_performance(KER3_Fe_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Fe_lm_origin_mp)
bptest(KER3_Fe_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Fe_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Fe_lm)

#Performance: Fe weighted lm forced through origin
KER3_Fe_wlm_origin <- lm(Fe_ICP ~ Fe-1, weight = 1/Fe_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Fe_wlm_origin) # summary stats
check_model(KER3_Fe_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Fe_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Fe_wlm_origin_mp <- model_performance(KER3_Fe_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Fe_wlm_origin_mp)
bptest(KER3_Fe_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Fe_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Fe_wlm)

# Co 

#Performance: Co unweighted lm stats forced through origin
KER3_Co_lm_origin <- lm(Co_ICP ~ Co-1, data = KER3_xrf_icp_matched)
summary(KER3_Co_lm_origin) # summary stats
check_model(KER3_Co_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Co_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Co_lm_origin_mp <- model_performance(KER3_Co_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Co_lm_origin_mp)
bptest(KER3_Co_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Co_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Co_lm)

#Performance: Co weighted lm forced through origin
KER3_Co_wlm_origin <- lm(Co_ICP ~ Co-1, weight = 1/Co_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Co_wlm_origin) # summary stats
check_model(KER3_Co_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Co_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Co_wlm_origin_mp <- model_performance(KER3_Co_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Co_wlm_origin_mp)
bptest(KER3_Co_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Co_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Co_wlm)

# Ni

#Performance: Ni unweighted lm stats forced through origin
KER3_Ni_lm_origin <- lm(Ni_ICP ~ Ni-1, data = KER3_xrf_icp_matched)
summary(KER3_Ni_lm_origin) # summary stats
check_model(KER3_Ni_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ni_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Ni_lm_origin_mp <- model_performance(KER3_Ni_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Ni_lm_origin_mp)
bptest(KER3_Ni_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ni_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Ni_lm)

#Performance: Ni weighted lm forced through origin
KER3_Ni_wlm_origin <- lm(Ni_ICP ~ Ni-1, weight = 1/Ni_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Ni_wlm_origin) # summary stats
check_model(KER3_Ni_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Ni_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Ni_wlm_origin_mp <- model_performance(KER3_Ni_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Ni_wlm_origin_mp)
bptest(KER3_Ni_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Ni_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Ni_wlm)

# Cu

#Performance: Cu unweighted lm stats forced through origin
KER3_Cu_lm_origin <- lm(Cu_ICP ~ Cu-1, data = KER3_xrf_icp_matched)
summary(KER3_Cu_lm_origin) # summary stats
check_model(KER3_Cu_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Cu_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Cu_lm_origin_mp <- model_performance(KER3_Cu_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Cu_lm_origin_mp)
bptest(KER3_Cu_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Cu_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Cu_lm)

#Performance: Cu weighted lm forced through origin
KER3_Cu_wlm_origin <- lm(Cu_ICP ~ Cu-1, weight = 1/Cu_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Cu_wlm_origin) # summary stats
check_model(KER3_Cu_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Cu_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Cu_wlm_origin_mp <- model_performance(KER3_Cu_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Cu_wlm_origin_mp)
bptest(KER3_Cu_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Cu_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Cu_wlm)

# Zn

#Performance: Zn unweighted lm stats forced through origin
KER3_Zn_lm_origin <- lm(Zn_ICP ~ Zn-1, data = KER3_xrf_icp_matched)
summary(KER3_Zn_lm_origin) # summary stats
check_model(KER3_Zn_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Zn_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Zn_lm_origin_mp <- model_performance(KER3_Zn_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Zn_lm_origin_mp)
bptest(KER3_Zn_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zn_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Zn_lm)

#Performance: Zn weighted lm forced through origin
KER3_Zn_wlm_origin <- lm(Zn_ICP ~ Zn-1, weight = 1/Zn_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Zn_wlm_origin) # summary stats
check_model(KER3_Zn_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Zn_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Zn_wlm_origin_mp <- model_performance(KER3_Zn_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Zn_wlm_origin_mp)
bptest(KER3_Zn_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Zn_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Zn_wlm)

# Rb

#Performance: Rb unweighted lm stats forced through origin
KER3_Rb_lm_origin <- lm(Rb_ICP ~ Rb-1, data = KER3_xrf_icp_matched)
summary(KER3_Rb_lm_origin) # summary stats
check_model(KER3_Rb_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Rb_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Rb_lm_origin_mp <- model_performance(KER3_Rb_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Rb_lm_origin_mp)
bptest(KER3_Rb_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Rb_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Rb_lm)

#Performance: Rb weighted lm forced through origin
KER3_Rb_wlm_origin <- lm(Rb_ICP ~ Rb-1, weight = 1/Rb_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Rb_wlm_origin) # summary stats
check_model(KER3_Rb_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Rb_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Rb_wlm_origin_mp <- model_performance(KER3_Rb_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Rb_wlm_origin_mp)
bptest(KER3_Rb_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Rb_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Rb_wlm)

# Sr

#Performance: Sr unweighted lm stats forced through origin
KER3_Sr_lm_origin <- lm(Sr_ICP ~ Sr-1, data = KER3_xrf_icp_matched)
summary(KER3_Sr_lm_origin) # summary stats
check_model(KER3_Sr_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Sr_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Sr_lm_origin_mp <- model_performance(KER3_Sr_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Sr_lm_origin_mp)
bptest(KER3_Sr_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Sr_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Sr_lm)

#Performance: Sr weighted lm forced through origin
KER3_Sr_wlm_origin <- lm(Sr_ICP ~ Sr-1, weight = 1/Sr_ICP_sd, data = KER3_xrf_icp_matched)
summary(KER3_Sr_wlm_origin) # summary stats
check_model(KER3_Sr_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/Sr_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_Sr_wlm_origin_mp <- model_performance(KER3_Sr_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_Sr_wlm_origin_mp)
bptest(KER3_Sr_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_Sr_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_Sr_wlm)

# DM

#Performance: DM unweighted lm stats forced through origin
KER3_DM_lm_origin <- lm(dry_mass_pc ~ coh_inc-1, data = KER3_xrf_icp_matched)
summary(KER3_DM_lm_origin) # summary stats
check_model(KER3_DM_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/DM_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_DM_lm_origin_mp <- model_performance(KER3_DM_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_DM_lm_origin_mp)
bptest(KER3_DM_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_lm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_DM_lm)

#Performance: DM weighted lm forced through origin
KER3_DM_wlm_origin <- lm(dry_mass_pc ~ coh_inc-1, weight = 1/dry_mass_err, data = KER3_xrf_icp_matched)
summary(KER3_DM_wlm_origin) # summary stats
check_model(KER3_DM_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/origin/DM_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
KER3_DM_wlm_origin_mp <- model_performance(KER3_DM_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(KER3_DM_wlm_origin_mp)
bptest(KER3_DM_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_wlm_origin) # Performance package summary check for heteroscedasticity
icc(KER3_DM_wlm)





# Summary 4x3 matrix plots of ITRAX-ACF matched variables --------
# OLS summary - forced through origin - unweighted stats on plot
ggarrange(KER3_K_origin, KER3_Ca_origin, KER3_Ti_origin, 
          KER3_Mn_origin, KER3_Fe_origin, KER3_Co_origin, 
          KER3_Ni_origin, KER3_Cu_origin, KER3_Zn_origin,
          KER3_Rb_origin, KER3_Sr_origin, KER3_DM_origin,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/KER3_OLS_Summary_origin_lm_wlm.pdf", 
       height = c(50), width = c(50), dpi = 600, units = "cm")







# -------------------------------------------------------------------------
# Old code ----------------------------------------------------------------

# Simple OLS regression plotting - unweighted stats on plot - no performance assessment -----------
theme_set(theme_bw(10))
KER3_K <- ggplot(KER3_xrf_icp_matched, aes(x = K, y = K_ICP)) +
  geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/K_ICP_sd)) + 
  stat_cor(aes(label =  after_stat(rr.label)), label.y.npc="bottom", label.x.npc = "right", hjust = 2) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "K (XRF-CS) / cps", y = "K (ICPMS) / ppm" ) +
  ggtitle("KER3: K (unweighted = black; ICP error-weighted = blue)")

KER3_Ca <- ggplot(KER3_xrf_icp_matched, aes(x = Ca, y = Ca_ICP)) +
  geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Ca_ICP_sd^2)) + 
  stat_cor(aes(label =  after_stat(rr.label)), label.y.npc="bottom", label.x.npc = "right", hjust = 2) +
  stat_poly_line(formula=y~x, colour = "black") +
  geom_point() +
  stat_poly_eq(formula=y~x, use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Ca (XRF-CS) / cps", y = "Ca (ICPMS) / ppm" ) +
  ggtitle("KER3: Ca (unweighted = black; ICP error-weighted = blue)")

KER3_Ti <- ggplot(KER3_xrf_icp_matched, aes(x = Ti, y = Ti_ICP)) +
  geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Ti_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Ti (XRF-CS) / cps", y = "Ti-49 (ICPMS) / ppm" ) +
  ggtitle("KER3: Ti (unweighted = black; ICP error-weighted = blue)")

KER3_Mn <- ggplot(KER3_xrf_icp_matched, aes(x = Mn, y = Mn_ICP)) +
  geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Mn_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Mn (XRF-CS) / cps", y = "Mn (ICPMS) / ppm" ) +
  ggtitle("KER3: Mn (unweighted = black; ICP error-weighted = blue)")

KER3_Fe <- ggplot(KER3_xrf_icp_matched, aes(x = Fe, y = Fe_ICP)) +
  geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Fe_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Fe (XRF-CS) / cps", y = "Fe (ICPMS) / ppm" ) +
  ggtitle("KER3: Fe (unweighted = black; ICP error-weighted = blue)")

KER3_Co <- ggplot(KER3_xrf_icp_matched, aes(x = Co, y = Co_ICP)) +
  geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Co_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Co (XRF-CS) / cps", y = "Co (ICPMS) / ppm" ) +
  ggtitle("KER3: Co (unweighted = black; ICP error-weighted = blue)")

KER3_Ni <- ggplot(KER3_xrf_icp_matched, aes(x = Ni, y = Ni_ICP)) +
  geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Ni_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Ni (XRF-CS) / cps", y = "Ni (ICPMS) / ppm" ) +
  ggtitle("KER3: Ni (unweighted = black; ICP error-weighted = blue)")

KER3_Cu <- ggplot(KER3_xrf_icp_matched, aes(x = Cu, y = Cu_ICP)) +
  geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Cu_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Cu (XRF-CS) / cps", y = "Cu (ICPMS) / ppm" ) +
  ggtitle("KER3: Cu (unweighted = black; ICP error-weighted = blue)")

KER3_Zn <- ggplot(KER3_xrf_icp_matched, aes(x = Zn, y = Zn_ICP)) +
  geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Zn_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Zn (XRF-CS) / cps", y = "Zn (ICPMS) / ppm" ) +
  ggtitle("KER3: Zn (unweighted = black; ICP error-weighted = blue)")

KER3_Rb <- ggplot(KER3_xrf_icp_matched, aes(x = Rb, y = Rb_ICP)) +
  geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Rb_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Rb (XRF-CS) / cps", y = "Rb (ICPMS) / ppm" ) +
  ggtitle("KER3: Rb (unweighted = black; ICP error-weighted = blue)")

KER3_Sr <- ggplot(KER3_xrf_icp_matched, aes(x = Sr, y = Sr_ICP)) +
  geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE,formula=y~x, aes(weight = 1/Sr_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Sr (XRF-CS) / cps", y = "Sr (ICPMS) / ppm" ) +
  ggtitle("KER3: Sr (unweighted = black; ICP error-weighted = blue)")

KER3_DM <- ggplot(KER3_xrf_icp_matched, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/dry_mass_err)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "coh/inc (XRF-CS)", y = "Dry mass / %" ) +
  ggtitle("KER3: Dry mass vs coh/inc (unweighted = black; ICP error-weighted = blue)")


# Simple version of OLS regression - forced through origin - unweighted stats on plot --------
theme_set(theme_bw(10))
KER3_K_origin <- ggplot(KER3_xrf_icp_matched, aes(x = K, y = K_ICP)) +
  geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/K_ICP_sd)) + 
  stat_cor(aes(label =  after_stat(rr.label)), label.y.npc="bottom", label.x.npc = "right", hjust = 2) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "K (XRF-CS) / cps", y = "K (ICPMS) / ppm" ) +
  ggtitle("KER3: K (unweighted = black; ICP error-weighted = blue)")

KER3_Ca_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Ca, y = Ca_ICP)) +
  geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Ca_ICP_sd^2)) + 
  stat_cor(aes(label =  after_stat(rr.label)), label.y.npc="bottom", label.x.npc = "right", hjust = 2) +
  stat_poly_line(formula=y~x, colour = "black") +
  geom_point() +
  stat_poly_eq(formula=y~x, use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Ca (XRF-CS) / cps", y = "Ca (ICPMS) / ppm" ) +
  ggtitle("KER3: Ca (unweighted = black; ICP error-weighted = blue)")

KER3_Ti_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Ti, y = Ti_ICP)) +
  geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Ti_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Ti (XRF-CS) / cps", y = "Ti-49 (ICPMS) / ppm" ) +
  ggtitle("KER3: Ti (unweighted = black; ICP error-weighted = blue)")

KER3_Mn_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Mn, y = Mn_ICP)) +
  geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Mn_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Mn (XRF-CS) / cps", y = "Mn (ICPMS) / ppm" ) +
  ggtitle("KER3: Mn (unweighted = black; ICP error-weighted = blue)")

KER3_Fe_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Fe, y = Fe_ICP)) +
  geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Fe_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Fe (XRF-CS) / cps", y = "Fe (ICPMS) / ppm" ) +
  ggtitle("KER3: Fe (unweighted = black; ICP error-weighted = blue)")

KER3_Co_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Co, y = Co_ICP)) +
  geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Co_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Co (XRF-CS) / cps", y = "Co (ICPMS) / ppm" ) +
  ggtitle("KER3: Co (unweighted = black; ICP error-weighted = blue)")

KER3_Ni_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Ni, y = Ni_ICP)) +
  geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Ni_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Ni (XRF-CS) / cps", y = "Ni (ICPMS) / ppm" ) +
  ggtitle("KER3: Ni (unweighted = black; ICP error-weighted = blue)")

KER3_Cu_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Cu, y = Cu_ICP)) +
  geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Cu_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Cu (XRF-CS) / cps", y = "Cu (ICPMS) / ppm" ) +
  ggtitle("KER3: Cu (unweighted = black; ICP error-weighted = blue)")

KER3_Zn_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Zn, y = Zn_ICP)) +
  geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Zn_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Zn (XRF-CS) / cps", y = "Zn (ICPMS) / ppm" ) +
  ggtitle("KER3: Zn (unweighted = black; ICP error-weighted = blue)")

KER3_Rb_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Rb, y = Rb_ICP)) +
  geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/Rb_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Rb (XRF-CS) / cps", y = "Rb (ICPMS) / ppm" ) +
  ggtitle("KER3: Rb (unweighted = black; ICP error-weighted = blue)")

KER3_Sr_origin <- ggplot(KER3_xrf_icp_matched, aes(x = Sr, y = Sr_ICP)) +
  geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE,formula=y~x, aes(weight = 1/Sr_ICP_sd)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "Sr (XRF-CS) / cps", y = "Sr (ICPMS) / ppm" ) +
  ggtitle("KER3: Sr (unweighted = black; ICP error-weighted = blue)")

KER3_DM_origin <- ggplot(KER3_xrf_icp_matched, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x, aes(weight = 1/dry_mass_err)) +
  stat_poly_line(colour = "black") +
  geom_point() +
  stat_poly_eq(use_label(c("R2", "adj.R2", "p", "eq", "f"))) +
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top") +
  scale_shape_manual(values = c(21)) +
  labs(x = "coh/inc (XRF-CS)", y = "Dry mass / %" ) +
  ggtitle("KER3: Dry mass vs coh/inc (unweighted = black; ICP error-weighted = blue)")

# Summary 4x3 matrix plots of ITRAX-ACF matched variables --------

# Simple OLS summary - unweighted stats on plot
ggarrange(KER3_K, KER3_Ca, KER3_Ti, 
          KER3_Mn, KER3_Fe, KER3_Co, 
          KER3_Ni, KER3_Cu, KER3_Zn,
          KER3_Rb, KER3_Sr, KER3_DM,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/KER3_OLS_Summary_lm_wlm.pdf", 
       height = c(50), width = c(50), dpi = 600, units = "cm")

# Simple OLS summary - forced through origin - unweighted stats on plot
ggarrange(KER3_K_origin, KER3_Ca_origin, KER3_Ti_origin, 
          KER3_Mn_origin, KER3_Fe_origin, KER3_Co_origin, 
          KER3_Ni_origin, KER3_Cu_origin, KER3_Zn_origin,
          KER3_Rb_origin, KER3_Sr_origin, KER3_DM_origin,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/KER3_OLS_Summary_origin_lm_wlm.pdf", 
       height = c(50), width = c(50), dpi = 600, units = "cm")





# -------------------------------------------------------------------------


