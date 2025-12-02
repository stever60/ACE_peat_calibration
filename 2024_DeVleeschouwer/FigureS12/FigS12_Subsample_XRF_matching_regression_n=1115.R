# Figure S12: Subsample - XRF dataset matching

# n = 1114 for DM dataset - 8/9/25 - new KER3 basal core alignment = original field depth 160-210 cm
# no PB1 LOI data

# Set up & clear previous ------------------------------------------------------

#clear previous console
remove (list = ls())
# Set working directory - Macbook Pro 2013
#setwd("/Users/Steve/Dropbox/BAS/Data/R/")
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
              'ggrepel', 'itraxR', 'PeriodicTable', 'errors', 'forecast', 'broom',
              'directlabels', 'performance', 'lmtest', 'ggpmisc', 'cowplot', 'Hmisc','car')
lapply(packages, library, character.only=TRUE)


# 1) Import Subsample data  ----------------------------------------------------
#  Define variable to select ---------------------------------------------------
subsample_param <- c("Water_Content", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density_g_cm3","Dry_density_g_cm3",	"Dry_density_err_g_cm3", 
                     "DMAR_g_cm2_yr", "DMAR_err_g_cm2_yr"	)

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density_g_cm3")

ACE_Subsample <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/ACE_Subsample.csv") %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all inf to NA
  #filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  select(c(Location:Section, Sample, strat_depth_top: accrate_err, Water_Content: DMAR_err_pc)) %>% 
  select(-Ash_Mass) %>% 
  rename(sample = Sample, top = strat_depth_top, bottom = strat_depth_bottom, 
         mp = strat_depth_mid) #rename columns for matching 
ACE_Subsample
head(ACE_Subsample)
tail(ACE_Subsample)

# 2) ITRAX XRF-CS cps qc dataset -----------------------------------------------
# Define lists of elements to use ----------------------------------------------

# XRF-CS acf elements matched to Francois ICPMS elements
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

# ICPMS elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# Import ITRAX XRF-CS cps qc dataset -----------------------------
ACE_itrax <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/ACE_xrf_qc.csv") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  # filter to retain data that passed quality control (qc) (should be all data)
  filter(qc == TRUE) %>%
  # filter for XRF-CS acf elements matched to Francois ICPMS elements
  select(!`Fe a*2`) %>% 
  select(Site, depth:surface, kcps:MSE, all_of(acf_icp_Elements_min), Mo_inc:coh_inc) %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% # create log scatter ratio columns
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc) %>%
  relocate(Total_scatter, .after = Ln_coh_inc)
ACE_itrax

# Create depth-matching XRF-CS files for each site ------------------

# BI10
ACE_itrax_BI10 <- ACE_itrax %>% 
  filter(Site == "BI10") %>% 
  select(depth:surface, kcps:Total_scatter)
# HER42PB
ACE_itrax_HER42PB <- ACE_itrax %>% 
  filter(Site == "HER42PB") %>% 
  select(depth:surface, kcps:Total_scatter)
# KER1
ACE_itrax_KER1 <- ACE_itrax %>% 
  filter(Site == "KER1") %>% 
  select(depth:surface, kcps:Total_scatter)
# KER3
ACE_itrax_KER3 <- ACE_itrax %>% 
  filter(Site == "KER3") %>% 
  select(depth:surface, kcps:Total_scatter)
# PB1
ACE_itrax_PB1 <- ACE_itrax %>%
  filter(Site == "PB1") %>% 
  select(depth:surface, kcps:Total_scatter)

tail(ACE_itrax)
tail(ACE_Subsample)


# 3) Depth-match Subsample-XRF datasets ----------------------------------------
# BI10 -------------------------------------------------------------------------

ACE_Subsample_BI10 <- ACE_Subsample %>% 
  filter(Site == "BI10") %>% 
  print()
write.csv(ACE_Subsample_BI10,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/BI10/ACE_Subsample_BI10.csv", 
          row.names = FALSE)

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
BI10_ss_xrf_matched <- itrax_reduce(ACE_itrax_BI10, names = ACE_Subsample_BI10$sample,
                                    breaks_lower = ACE_Subsample_BI10$top,
                                    breaks_upper = ACE_Subsample_BI10$bottom,
                                    fun = mean,
                                    edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_BI10, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
BI10_ss_xrf_matched

# Calculate stdev for ITRAX matching data
BI10_ss_xrf_matched_sd <- itrax_reduce(ACE_itrax_BI10, names = ACE_Subsample_BI10$sample,
                                       breaks_lower = ACE_Subsample_BI10$top,
                                       breaks_upper = ACE_Subsample_BI10$bottom,
                                       fun = sd,
                                       edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_BI10, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
BI10_ss_xrf_matched_sd

# Merge mean & sd data into the same dataframe
BI10_Subsample_xrf_matched0 <-  BI10_ss_xrf_matched_sd %>% 
  inner_join(., BI10_ss_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
BI10_Subsample_xrf_matched0

BI10_Subsample_xrf_matched <- BI10_Subsample_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
BI10_Subsample_xrf_matched
BI10_Subsample_xrf_matched$Zn # check no 0 values present

write.csv(BI10_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/BI10/BI10_Subsample_xrf_matched.csv", 
          row.names = FALSE)

# HER42PB ----------------------------------------------------------------------

ACE_Subsample_HER42PB <- ACE_Subsample %>% 
  filter(Site == "HER42PB") %>% 
  print()
write.csv(ACE_Subsample_HER42PB,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/HER42PB/ACE_Subsample_HER42PB.csv", 
          row.names = FALSE)

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
HER42PB_ss_xrf_matched <- itrax_reduce(ACE_itrax_HER42PB, names = ACE_Subsample_HER42PB$sample,
                                      breaks_lower = ACE_Subsample_HER42PB$top,
                                      breaks_upper = ACE_Subsample_HER42PB$bottom,
                                      fun = mean,
                                      edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_HER42PB, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
HER42PB_ss_xrf_matched

# Calculate stdev for ITRAX matching data
HER42PB_ss_xrf_matched_sd <- itrax_reduce(ACE_itrax_HER42PB, names = ACE_Subsample_HER42PB$sample,
                                         breaks_lower = ACE_Subsample_HER42PB$top,
                                         breaks_upper = ACE_Subsample_HER42PB$bottom,
                                         fun = sd,
                                         edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_HER42PB, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
HER42PB_ss_xrf_matched_sd

# Merge mean & sd data into the same dataframe
HER42PB_Subsample_xrf_matched0 <-  HER42PB_ss_xrf_matched_sd %>% 
  inner_join(., HER42PB_ss_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
HER42PB_Subsample_xrf_matched0

HER42PB_Subsample_xrf_matched <- HER42PB_Subsample_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
HER42PB_Subsample_xrf_matched
HER42PB_Subsample_xrf_matched$Zn # check no 0 values present

write.csv(HER42PB_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/HER42PB/HER42PB_Subsample_xrf_matched.csv", 
          row.names = FALSE)

# KER1 -------------------------------------------------------------------------

ACE_Subsample_KER1 <- ACE_Subsample %>% 
  filter(Site == "KER1") %>% 
  print()
write.csv(ACE_Subsample_KER1,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/KER1/ACE_Subsample_KER1.csv", 
          row.names = FALSE)

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
KER1_ss_xrf_matched <- itrax_reduce(ACE_itrax_KER1, names = ACE_Subsample_KER1$sample,
                                       breaks_lower = ACE_Subsample_KER1$top,
                                       breaks_upper = ACE_Subsample_KER1$bottom,
                                       fun = mean,
                                       edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_KER1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
KER1_ss_xrf_matched

# Calculate stdev for ITRAX matching data
KER1_ss_xrf_matched_sd <- itrax_reduce(ACE_itrax_KER1, names = ACE_Subsample_KER1$sample,
                                          breaks_lower = ACE_Subsample_KER1$top,
                                          breaks_upper = ACE_Subsample_KER1$bottom,
                                          fun = sd,
                                          edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_KER1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
KER1_ss_xrf_matched_sd

# Merge mean & sd data into the same dataframe
KER1_Subsample_xrf_matched0 <-  KER1_ss_xrf_matched_sd %>% 
  inner_join(., KER1_ss_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
KER1_Subsample_xrf_matched0

KER1_Subsample_xrf_matched <- KER1_Subsample_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
KER1_Subsample_xrf_matched
KER1_Subsample_xrf_matched$Dry_mass # check no 0 values present

write.csv(KER1_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/KER1/KER1_Subsample_xrf_matched.csv", 
          row.names = FALSE)

# KER3 -------------------------------------------------------------------------

ACE_Subsample_KER3 <- ACE_Subsample %>% 
  filter(Site == "KER3") %>% 
  print()
write.csv(ACE_Subsample_KER3,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/KER3/ACE_Subsample_KER3.csv", 
          row.names = FALSE)

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
KER3_ss_xrf_matched <- itrax_reduce(ACE_itrax_KER3, names = ACE_Subsample_KER3$sample,
                                       breaks_lower = ACE_Subsample_KER3$top,
                                       breaks_upper = ACE_Subsample_KER3$bottom,
                                       fun = mean,
                                       edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_KER3, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
KER3_ss_xrf_matched

# Calculate stdev for ITRAX matching data
KER3_ss_xrf_matched_sd <- itrax_reduce(ACE_itrax_KER3, names = ACE_Subsample_KER3$sample,
                                          breaks_lower = ACE_Subsample_KER3$top,
                                          breaks_upper = ACE_Subsample_KER3$bottom,
                                          fun = sd,
                                          edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_KER3, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
KER3_ss_xrf_matched_sd

# Merge mean & sd data into the same dataframe
KER3_Subsample_xrf_matched0 <-  KER3_ss_xrf_matched_sd %>% 
  inner_join(., KER3_ss_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
KER3_Subsample_xrf_matched0

KER3_Subsample_xrf_matched <- KER3_Subsample_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
KER3_Subsample_xrf_matched
KER3_Subsample_xrf_matched$Zn # check no 0 values present

write.csv(KER3_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/KER3/KER3_Subsample_xrf_matched.csv", 
          row.names = FALSE)

# PB1 ----------------------------------------------------------------

ACE_Subsample_PB1 <- ACE_Subsample %>% 
  filter(Site == "PB1") %>% 
  print()
write.csv(ACE_Subsample_PB1,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/PB1/ACE_Subsample_PB1.csv", 
          row.names = FALSE)

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
PB1_ss_xrf_matched <- itrax_reduce(ACE_itrax_PB1, names = ACE_Subsample_PB1$sample,
                                    breaks_lower = ACE_Subsample_PB1$top,
                                    breaks_upper = ACE_Subsample_PB1$bottom,
                                    fun = mean,
                                    edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_PB1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
PB1_ss_xrf_matched
PB1_ss_xrf_matched$Zn

# Calculate stdev for ITRAX matching data
PB1_ss_xrf_matched_sd <- itrax_reduce(ACE_itrax_PB1, names = ACE_Subsample_PB1$sample,
                                       breaks_lower = ACE_Subsample_PB1$top,
                                       breaks_upper = ACE_Subsample_PB1$bottom,
                                       fun = sd,
                                       edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_PB1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
PB1_ss_xrf_matched_sd
PB1_ss_xrf_matched_sd$Zn_sd

# Merge mean & sd data into the same dataframe
PB1_Subsample_xrf_matched0 <-  PB1_ss_xrf_matched_sd %>% 
  inner_join(., PB1_ss_xrf_matched, by = "sample") %>% 
  #filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
PB1_Subsample_xrf_matched0
PB1_Subsample_xrf_matched0$Zn

PB1_Subsample_xrf_matched <- PB1_Subsample_xrf_matched0 %>% # replace 0 values with half min value
  mutate(across(
    all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd)),
    ~ if_else(. == 0,
              min(.[. != 0], na.rm = TRUE) / 2,.)))
PB1_Subsample_xrf_matched
PB1_Subsample_xrf_matched$Zn # check no 0 values present

write.csv(PB1_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/PB1/PB1_Subsample_xrf_matched.csv", 
          row.names = FALSE)


# Create ACE Subsample-XRF matched dataset ----------------------------

ACE_Subsample_xrf_matched<- bind_rows(BI10_Subsample_xrf_matched,
                                HER42PB_Subsample_xrf_matched,
                                KER1_Subsample_xrf_matched,
                                KER3_Subsample_xrf_matched,
                                PB1_Subsample_xrf_matched)
ACE_Subsample_xrf_matched

write.csv(ACE_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/ACE/ACE_Subsample_xrf_matched.csv", 
          row.names = FALSE)



# 4) Dry Mass vs coh/inc regression --------------------------------------------
# ACE OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
ACE_DM_lm <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched)
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_summary.txt")
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_wlm <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2)
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_wt_summary.txt")
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model

# This doesnt work because for Subsampled dataframe as it contains NA - means its a different length
#ACE_model <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched) # define model
#ACE_wt <- 1 / lm(abs(ACE_model$residuals) ~ ACE_model$fitted.values)$fitted.values^2 #define weights to use
#ACE_wls <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched, weights=ACE_wt) #perform weighted least squares regression

# 1. Fit initial model
ACE_model <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched)

# 2. Extract the model frame actually used by lm()
#    This drops rows with NA exactly like the model did
mf <- model.frame(ACE_model)

# 3. Build variance model using this model frame
var_model <- lm(abs(residuals(ACE_model)) ~ fitted(ACE_model))

# 4. Compute weights (1115 rows)
wt <- 1 / (fitted(var_model)^2)

# 5. Fit WLS model using the same model frame
ACE_wls <- lm(Dry_mass ~ coh_inc,
              data = mf,
              weights = wt)

# Checks
summary(ACE_wls) # summary stats
glance(ACE_wls) # summary stats including AIC
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_summary.txt")
summary(ACE_wls)
glance(ACE_wls)
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls) # Performance package summary check for heteroscedasticity
icc(ACE_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
# This doesnt work because for Subsampled dataframe as it contains NA - means its a different length
# ACE_model_wt <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched, weight = 1/Dry_mass_err^2) # define model
# ACE_wt_wt <- 1 / lm(abs(ACE_model_wt$residuals) ~ ACE_model_wt$fitted.values)$fitted.values^2 #define weights to use
# ACE_wls_wt <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched, weights=ACE_wt_wt) #perform weighted least squares regression
# 1. Fit initial model (OLS)
ACE_model <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched)

# 2. Extract the rows actually used by lm()
mf <- model.frame(ACE_model)

# 3. Fit variance model on those rows
var_model <- lm(
  abs(residuals(ACE_model)) ~ fitted(ACE_model)
)

# 4. Compute model-appropriate weights (length = number of rows used)
ACE_wt_internal <- 1 / (fitted(var_model)^2)

# 5. Create full-length vector (length = nrow(ACE_Subsample_xrf_matched))
ACE_Subsample_xrf_matched$ACE_wt <- NA_real_

# 6. Insert weights into the correct row indices
ACE_Subsample_xrf_matched$ACE_wt[as.numeric(rownames(mf))] <- ACE_wt_internal

length(ACE_Subsample_xrf_matched$ACE_wt) == 1135

# Checks
summary(ACE_wls_wt) # summary stats
glance(ACE_wls_wt) # summary stats including AIC
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_wt_summary.txt")
summary(ACE_wls_wt)
glance(ACE_wls_wt)
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_wls_wt) # check for random effects - returns NULL if none present
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
ACE_DM_lm_outliers <- ACE_Subsample_xrf_matched[ACE_DM_lm_influential_names,] # outliers only using of index values
ACE_DM_lm_no_outliers <- ACE_Subsample_xrf_matched %>% anti_join(ACE_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_DM_lm_lev.pdf")
leveragePlots(ACE_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wlm_hats <- as.data.frame(hatvalues(ACE_wlm))
ACE_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wlm_cooksD <- cooks.distance(ACE_wlm)
ACE_wlm_influential <- ACE_wlm_cooksD[(ACE_wlm_cooksD > (3 * mean(ACE_wlm_cooksD, na.rm = TRUE)))]
ACE_wlm_influential
ACE_wlm_influential_names <- names(ACE_wlm_influential)
ACE_wlm_outliers <- ACE_Subsample_xrf_matched[ACE_wlm_influential_names,] # outliers only using of index values
ACE_wlm_no_outliers <- ACE_Subsample_xrf_matched %>% anti_join(ACE_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wlm_lev.pdf")
leveragePlots(ACE_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_hats <- as.data.frame(hatvalues(ACE_wls))
ACE_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_cooksD <- cooks.distance(ACE_wls)
ACE_wls_influential <- ACE_wls_cooksD[(ACE_wls_cooksD > (3 * mean(ACE_wls_cooksD, na.rm = TRUE)))]
ACE_wls_influential
ACE_wls_influential_names <- names(ACE_wls_influential)
ACE_wls_outliers <- ACE_Subsample_xrf_matched[ACE_wls_influential_names,] # outliers only using of index values
ACE_wls_no_outliers <- ACE_Subsample_xrf_matched %>% anti_join(ACE_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_lev_bar.pdf")
barplot(hatvalues(ACE_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_lev.pdf")
leveragePlots(ACE_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_wt_hats <- as.data.frame(hatvalues(ACE_wls_wt))
ACE_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_wt_cooksD <- cooks.distance(ACE_wls_wt)
ACE_wls_wt_influential <- ACE_wls_wt_cooksD[(ACE_wls_wt_cooksD > (3 * mean(ACE_wls_wt_cooksD, na.rm = TRUE)))]
ACE_wls_wt_influential
ACE_wls_wt_influential_names <- names(ACE_wls_wt_influential)
ACE_wls_wt_outliers <- ACE_Subsample_xrf_matched[ACE_wls_wt_influential_names,] # outliers only using of index values
ACE_wls_wt_no_outliers <- ACE_Subsample_xrf_matched %>% anti_join(ACE_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_wt_lev.pdf")
leveragePlots(ACE_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "ACE"
element_title <- "coh_inc"
theme_set(theme_classic(10))
ACE <- ggplot(ACE_Subsample_xrf_matched, aes(x = coh_inc, y = Dry_mass)) +
  geom_errorbar(aes(ymin=Dry_mass-Dry_mass_err, ymax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_point(fill = "Site", color = "Site", size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/Dry_mass_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = ACE_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = ACE_wt_wt), colour="black") + # weighted WLS regression
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
  #xlim(0.15, 0.3) +
  #ylim(0, 80) +
  ggtitle(paste("Site: ", site_title, ": Dry mass vs ", element_title))
ACE
# Define p value, OLS equation & R2 as a string to add to plot

ACE_DM_lm_p <- function(ACE_DM_lm) {
  f <- summary(ACE_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_lm_p(ACE_DM_lm)

ACE_DM_lm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_lm_p(ACE_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_wlm_p <- function(ACE_wlm) {
  f <- summary(ACE_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wlm_p(ACE_wlm)

ACE_wlm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wlm_p(ACE_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_wls_p <- function(ACE_wls) {
  f <- summary(ACE_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_p(ACE_wls)

ACE_wls_eqn <- function(df){
  m <- ACE_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_p(ACE_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_wls_wt_p <- function(ACE_wls_wt) {
  f <- summary(ACE_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_wt_p(ACE_wls_wt)

ACE_wls_wt_eqn <- function(df){
  m <- ACE_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_wt_p(ACE_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_final <- ACE + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/FigS12_ACE_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")



# BI10 OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
BI10_DM_lm <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched)
summary(BI10_DM_lm)
glance(BI10_DM_lm)
model_performance(BI10_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(BI10_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_OLS_summary.txt")
summary(BI10_DM_lm)
glance(BI10_DM_lm)
model_performance(BI10_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(BI10_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(BI10_DM_lm) # Performance package summary check for heteroscedasticity
icc(BI10_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
BI10_wlm <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2)
summary(BI10_wlm)
glance(BI10_wlm)
model_performance(BI10_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(BI10_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_OLS_wt_summary.txt")
summary(BI10_wlm)
glance(BI10_wlm)
model_performance(BI10_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(BI10_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(BI10_DM_lm) # Performance package summary check for heteroscedasticity
icc(BI10_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
BI10_model <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched) # define model
BI10_wt <- 1 / lm(abs(BI10_model$residuals) ~ BI10_model$fitted.values)$fitted.values^2 #define weights to use
BI10_wls <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched, weights=BI10_wt) #perform weighted least squares regression
# Checks
summary(BI10_wls) # summary stats
glance(BI10_wls) # summary stats including AIC
model_performance(BI10_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(BI10_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_WLS_summary.txt")
summary(BI10_wls)
glance(BI10_wls)
model_performance(BI10_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(BI10_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(BI10_wls) # Performance package summary check for heteroscedasticity
icc(BI10_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
BI10_model_wt <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched, weight = 1/Dry_mass_err^2) # define model
BI10_wt_wt <- 1 / lm(abs(BI10_model_wt$residuals) ~ BI10_model_wt$fitted.values)$fitted.values^2 #define weights to use
BI10_wls_wt <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched, weights=BI10_wt_wt) #perform weighted least squares regression
# Checks
summary(BI10_wls_wt) # summary stats
glance(BI10_wls_wt) # summary stats including AIC
model_performance(BI10_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(BI10_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_WLS_wt_summary.txt")
summary(BI10_wls_wt)
glance(BI10_wls_wt)
model_performance(BI10_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(BI10_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(BI10_wls_wt) # Performance package summary check for heteroscedasticity
icc(BI10_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
BI10_DM_lm_hats <- as.data.frame(hatvalues(BI10_DM_lm))
BI10_DM_lm_hats
# Cooks distance - 2-3 x difference from mean 
BI10_DM_lm_cooksD <- cooks.distance(BI10_DM_lm)
BI10_DM_lm_influential <- BI10_DM_lm_cooksD[(BI10_DM_lm_cooksD > (3 * mean(BI10_DM_lm_cooksD, na.rm = TRUE)))]
BI10_DM_lm_influential
BI10_DM_lm_influential_names <- names(BI10_DM_lm_influential)
BI10_DM_lm_outliers <- BI10_Subsample_xrf_matched[BI10_DM_lm_influential_names,] # outliers only using of index values
BI10_DM_lm_no_outliers <- BI10_Subsample_xrf_matched %>% anti_join(BI10_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(BI10_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(BI10_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_DM_lm_lev.pdf")
leveragePlots(BI10_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(BI10_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
BI10_wlm_hats <- as.data.frame(hatvalues(BI10_wlm))
BI10_wlm_hats
# Cooks distance - 2-3 x difference from mean 
BI10_wlm_cooksD <- cooks.distance(BI10_wlm)
BI10_wlm_influential <- BI10_wlm_cooksD[(BI10_wlm_cooksD > (3 * mean(BI10_wlm_cooksD, na.rm = TRUE)))]
BI10_wlm_influential
BI10_wlm_influential_names <- names(BI10_wlm_influential)
BI10_wlm_outliers <- BI10_Subsample_xrf_matched[BI10_wlm_influential_names,] # outliers only using of index values
BI10_wlm_no_outliers <- BI10_Subsample_xrf_matched %>% anti_join(BI10_wlm_outliers) # generates a new dataset with outliers removed
write.csv(BI10_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wlm_lev_bar.pdf")
barplot(hatvalues(BI10_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wlm_lev.pdf")
leveragePlots(BI10_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(BI10_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
BI10_wls_hats <- as.data.frame(hatvalues(BI10_wls))
BI10_wls_hats
# Cooks distance - 2-3 x difference from mean 
BI10_wls_cooksD <- cooks.distance(BI10_wls)
BI10_wls_influential <- BI10_wls_cooksD[(BI10_wls_cooksD > (3 * mean(BI10_wls_cooksD, na.rm = TRUE)))]
BI10_wls_influential
BI10_wls_influential_names <- names(BI10_wls_influential)
BI10_wls_outliers <- BI10_Subsample_xrf_matched[BI10_wls_influential_names,] # outliers only using of index values
BI10_wls_no_outliers <- BI10_Subsample_xrf_matched %>% anti_join(BI10_wls_outliers) # generates a new dataset with outliers removed
write.csv(BI10_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wls_lev_bar.pdf")
barplot(hatvalues(BI10_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wls_lev.pdf")
leveragePlots(BI10_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(BI10_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
BI10_wls_wt_hats <- as.data.frame(hatvalues(BI10_wls_wt))
BI10_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
BI10_wls_wt_cooksD <- cooks.distance(BI10_wls_wt)
BI10_wls_wt_influential <- BI10_wls_wt_cooksD[(BI10_wls_wt_cooksD > (3 * mean(BI10_wls_wt_cooksD, na.rm = TRUE)))]
BI10_wls_wt_influential
BI10_wls_wt_influential_names <- names(BI10_wls_wt_influential)
BI10_wls_wt_outliers <- BI10_Subsample_xrf_matched[BI10_wls_wt_influential_names,] # outliers only using of index values
BI10_wls_wt_no_outliers <- BI10_Subsample_xrf_matched %>% anti_join(BI10_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(BI10_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(BI10_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wls_wt_lev.pdf")
leveragePlots(BI10_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(BI10_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "BI10"
element_title <- "coh_inc"
theme_set(theme_classic(10))
BI10 <- ggplot(BI10_Subsample_xrf_matched, aes(x = coh_inc, y = Dry_mass)) +
  geom_errorbar(aes(ymin=Dry_mass-Dry_mass_err, ymax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/Dry_mass_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = BI10_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = BI10_wt_wt), colour="black") + # weighted WLS regression
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
  #xlim(0.15, 0.25) +
  #ylim(10, 80) +
  ggtitle(paste("Site: ", site_title, ": Dry mass vs ", element_title))
BI10
# Define p value, OLS equation & R2 as a string to add to plot

BI10_DM_lm_p <- function(BI10_DM_lm) {
  f <- summary(BI10_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
BI10_DM_lm_p(BI10_DM_lm)

BI10_DM_lm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(BI10_DM_lm_p(BI10_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

BI10_wlm_p <- function(BI10_wlm) {
  f <- summary(BI10_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
BI10_wlm_p(BI10_wlm)

BI10_wlm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = BI10_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(BI10_wlm_p(BI10_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

BI10_wls_p <- function(BI10_wls) {
  f <- summary(BI10_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
BI10_wls_p(BI10_wls)

BI10_wls_eqn <- function(df){
  m <- BI10_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(BI10_wls_p(BI10_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

BI10_wls_wt_p <- function(BI10_wls_wt) {
  f <- summary(BI10_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
BI10_wls_wt_p(BI10_wls_wt)

BI10_wls_wt_eqn <- function(df){
  m <- BI10_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(BI10_wls_wt_p(BI10_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
BI10_final <- BI10 + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", BI10_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", BI10_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", BI10_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", BI10_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
BI10_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/BI10/FigS12_BI10_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# HER42PB OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
HER42PB_DM_lm <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched)
summary(HER42PB_DM_lm)
glance(HER42PB_DM_lm)
model_performance(HER42PB_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(HER42PB_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_OLS_summary.txt")
summary(HER42PB_DM_lm)
glance(HER42PB_DM_lm)
model_performance(HER42PB_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(HER42PB_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(HER42PB_DM_lm) # Performance package summary check for heteroscedasticity
icc(HER42PB_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
HER42PB_wlm <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2)
summary(HER42PB_wlm)
glance(HER42PB_wlm)
model_performance(HER42PB_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(HER42PB_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_OLS_wt_summary.txt")
summary(HER42PB_wlm)
glance(HER42PB_wlm)
model_performance(HER42PB_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(HER42PB_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(HER42PB_DM_lm) # Performance package summary check for heteroscedasticity
icc(HER42PB_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
HER42PB_model <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched) # define model
HER42PB_wt <- 1 / lm(abs(HER42PB_model$residuals) ~ HER42PB_model$fitted.values)$fitted.values^2 #define weights to use
HER42PB_wls <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched, weights=HER42PB_wt) #perform weighted least squares regression
# Checks
summary(HER42PB_wls) # summary stats
glance(HER42PB_wls) # summary stats including AIC
model_performance(HER42PB_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(HER42PB_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_WLS_summary.txt")
summary(HER42PB_wls)
glance(HER42PB_wls)
model_performance(HER42PB_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(HER42PB_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(HER42PB_wls) # Performance package summary check for heteroscedasticity
icc(HER42PB_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
HER42PB_model_wt <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched, weight = 1/Dry_mass_err^2) # define model
HER42PB_wt_wt <- 1 / lm(abs(HER42PB_model_wt$residuals) ~ HER42PB_model_wt$fitted.values)$fitted.values^2 #define weights to use
HER42PB_wls_wt <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched, weights=HER42PB_wt_wt) #perform weighted least squares regression
# Checks
summary(HER42PB_wls_wt) # summary stats
glance(HER42PB_wls_wt) # summary stats including AIC
model_performance(HER42PB_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(HER42PB_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_WLS_wt_summary.txt")
summary(HER42PB_wls_wt)
glance(HER42PB_wls_wt)
model_performance(HER42PB_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(HER42PB_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(HER42PB_wls_wt) # Performance package summary check for heteroscedasticity
icc(HER42PB_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
HER42PB_DM_lm_hats <- as.data.frame(hatvalues(HER42PB_DM_lm))
HER42PB_DM_lm_hats
# Cooks distance - 2-3 x difference from mean 
HER42PB_DM_lm_cooksD <- cooks.distance(HER42PB_DM_lm)
HER42PB_DM_lm_influential <- HER42PB_DM_lm_cooksD[(HER42PB_DM_lm_cooksD > (3 * mean(HER42PB_DM_lm_cooksD, na.rm = TRUE)))]
HER42PB_DM_lm_influential
HER42PB_DM_lm_influential_names <- names(HER42PB_DM_lm_influential)
HER42PB_DM_lm_outliers <- HER42PB_Subsample_xrf_matched[HER42PB_DM_lm_influential_names,] # outliers only using of index values
HER42PB_DM_lm_no_outliers <- HER42PB_Subsample_xrf_matched %>% anti_join(HER42PB_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(HER42PB_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(HER42PB_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_DM_lm_lev.pdf")
leveragePlots(HER42PB_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(HER42PB_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
HER42PB_wlm_hats <- as.data.frame(hatvalues(HER42PB_wlm))
HER42PB_wlm_hats
# Cooks distance - 2-3 x difference from mean 
HER42PB_wlm_cooksD <- cooks.distance(HER42PB_wlm)
HER42PB_wlm_influential <- HER42PB_wlm_cooksD[(HER42PB_wlm_cooksD > (3 * mean(HER42PB_wlm_cooksD, na.rm = TRUE)))]
HER42PB_wlm_influential
HER42PB_wlm_influential_names <- names(HER42PB_wlm_influential)
HER42PB_wlm_outliers <- HER42PB_Subsample_xrf_matched[HER42PB_wlm_influential_names,] # outliers only using of index values
HER42PB_wlm_no_outliers <- HER42PB_Subsample_xrf_matched %>% anti_join(HER42PB_wlm_outliers) # generates a new dataset with outliers removed
write.csv(HER42PB_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wlm_lev_bar.pdf")
barplot(hatvalues(HER42PB_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wlm_lev.pdf")
leveragePlots(HER42PB_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(HER42PB_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
HER42PB_wls_hats <- as.data.frame(hatvalues(HER42PB_wls))
HER42PB_wls_hats
# Cooks distance - 2-3 x difference from mean 
HER42PB_wls_cooksD <- cooks.distance(HER42PB_wls)
HER42PB_wls_influential <- HER42PB_wls_cooksD[(HER42PB_wls_cooksD > (3 * mean(HER42PB_wls_cooksD, na.rm = TRUE)))]
HER42PB_wls_influential
HER42PB_wls_influential_names <- names(HER42PB_wls_influential)
HER42PB_wls_outliers <- HER42PB_Subsample_xrf_matched[HER42PB_wls_influential_names,] # outliers only using of index values
HER42PB_wls_no_outliers <- HER42PB_Subsample_xrf_matched %>% anti_join(HER42PB_wls_outliers) # generates a new dataset with outliers removed
write.csv(HER42PB_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wls_lev_bar.pdf")
barplot(hatvalues(HER42PB_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wls_lev.pdf")
leveragePlots(HER42PB_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(HER42PB_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
HER42PB_wls_wt_hats <- as.data.frame(hatvalues(HER42PB_wls_wt))
HER42PB_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
HER42PB_wls_wt_cooksD <- cooks.distance(HER42PB_wls_wt)
HER42PB_wls_wt_influential <- HER42PB_wls_wt_cooksD[(HER42PB_wls_wt_cooksD > (3 * mean(HER42PB_wls_wt_cooksD, na.rm = TRUE)))]
HER42PB_wls_wt_influential
HER42PB_wls_wt_influential_names <- names(HER42PB_wls_wt_influential)
HER42PB_wls_wt_outliers <- HER42PB_Subsample_xrf_matched[HER42PB_wls_wt_influential_names,] # outliers only using of index values
HER42PB_wls_wt_no_outliers <- HER42PB_Subsample_xrf_matched %>% anti_join(HER42PB_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(HER42PB_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(HER42PB_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wls_wt_lev.pdf")
leveragePlots(HER42PB_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(HER42PB_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "HER42PB"
element_title <- "coh_inc"
theme_set(theme_classic(10))
HER42PB <- ggplot(HER42PB_Subsample_xrf_matched, aes(x = coh_inc, y = Dry_mass)) +
  geom_errorbar(aes(ymin=Dry_mass-Dry_mass_err, ymax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_point(fill = "gold", color = 'gold', size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/Dry_mass_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = HER42PB_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = HER42PB_wt_wt), colour="black") + # weighted WLS regression
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
  xlim(0.15, 0.25) +
  ylim(10, 80) +
  ggtitle(paste("Site: ", site_title, ": Dry mass vs ", element_title))
HER42PB
# Define p value, OLS equation & R2 as a string to add to plot

HER42PB_DM_lm_p <- function(HER42PB_DM_lm) {
  f <- summary(HER42PB_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
HER42PB_DM_lm_p(HER42PB_DM_lm)

HER42PB_DM_lm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(HER42PB_DM_lm_p(HER42PB_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

HER42PB_wlm_p <- function(HER42PB_wlm) {
  f <- summary(HER42PB_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
HER42PB_wlm_p(HER42PB_wlm)

HER42PB_wlm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = HER42PB_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(HER42PB_wlm_p(HER42PB_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

HER42PB_wls_p <- function(HER42PB_wls) {
  f <- summary(HER42PB_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
HER42PB_wls_p(HER42PB_wls)

HER42PB_wls_eqn <- function(df){
  m <- HER42PB_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(HER42PB_wls_p(HER42PB_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

HER42PB_wls_wt_p <- function(HER42PB_wls_wt) {
  f <- summary(HER42PB_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
HER42PB_wls_wt_p(HER42PB_wls_wt)

HER42PB_wls_wt_eqn <- function(df){
  m <- HER42PB_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(HER42PB_wls_wt_p(HER42PB_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
HER42PB_final <- HER42PB + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", HER42PB_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", HER42PB_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", HER42PB_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", HER42PB_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
HER42PB_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/HER42PB/FigS12_HER42PB_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# KER1 OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
KER1_DM_lm <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched)
summary(KER1_DM_lm)
glance(KER1_DM_lm)
model_performance(KER1_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER1_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_OLS_summary.txt")
summary(KER1_DM_lm)
glance(KER1_DM_lm)
model_performance(KER1_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER1_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER1_DM_lm) # Performance package summary check for heteroscedasticity
icc(KER1_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
KER1_wlm <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2)
summary(KER1_wlm)
glance(KER1_wlm)
model_performance(KER1_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER1_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_OLS_wt_summary.txt")
summary(KER1_wlm)
glance(KER1_wlm)
model_performance(KER1_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER1_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER1_DM_lm) # Performance package summary check for heteroscedasticity
icc(KER1_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
KER1_model <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched) # define model
KER1_wt <- 1 / lm(abs(KER1_model$residuals) ~ KER1_model$fitted.values)$fitted.values^2 #define weights to use
KER1_wls <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched, weights=KER1_wt) #perform weighted least squares regression
# Checks
summary(KER1_wls) # summary stats
glance(KER1_wls) # summary stats including AIC
model_performance(KER1_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER1_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_WLS_summary.txt")
summary(KER1_wls)
glance(KER1_wls)
model_performance(KER1_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER1_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER1_wls) # Performance package summary check for heteroscedasticity
icc(KER1_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
KER1_model_wt <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched, weight = 1/Dry_mass_err^2) # define model
KER1_wt_wt <- 1 / lm(abs(KER1_model_wt$residuals) ~ KER1_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER1_wls_wt <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched, weights=KER1_wt_wt) #perform weighted least squares regression
# Checks
summary(KER1_wls_wt) # summary stats
glance(KER1_wls_wt) # summary stats including AIC
model_performance(KER1_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER1_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_WLS_wt_summary.txt")
summary(KER1_wls_wt)
glance(KER1_wls_wt)
model_performance(KER1_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER1_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER1_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER1_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER1_DM_lm_hats <- as.data.frame(hatvalues(KER1_DM_lm))
KER1_DM_lm_hats
# Cooks distance - 2-3 x difference from mean 
KER1_DM_lm_cooksD <- cooks.distance(KER1_DM_lm)
KER1_DM_lm_influential <- KER1_DM_lm_cooksD[(KER1_DM_lm_cooksD > (3 * mean(KER1_DM_lm_cooksD, na.rm = TRUE)))]
KER1_DM_lm_influential
KER1_DM_lm_influential_names <- names(KER1_DM_lm_influential)
KER1_DM_lm_outliers <- KER1_Subsample_xrf_matched[KER1_DM_lm_influential_names,] # outliers only using of index values
KER1_DM_lm_no_outliers <- KER1_Subsample_xrf_matched %>% anti_join(KER1_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(KER1_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(KER1_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_DM_lm_lev.pdf")
leveragePlots(KER1_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER1_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER1_wlm_hats <- as.data.frame(hatvalues(KER1_wlm))
KER1_wlm_hats
# Cooks distance - 2-3 x difference from mean 
KER1_wlm_cooksD <- cooks.distance(KER1_wlm)
KER1_wlm_influential <- KER1_wlm_cooksD[(KER1_wlm_cooksD > (3 * mean(KER1_wlm_cooksD, na.rm = TRUE)))]
KER1_wlm_influential
KER1_wlm_influential_names <- names(KER1_wlm_influential)
KER1_wlm_outliers <- KER1_Subsample_xrf_matched[KER1_wlm_influential_names,] # outliers only using of index values
KER1_wlm_no_outliers <- KER1_Subsample_xrf_matched %>% anti_join(KER1_wlm_outliers) # generates a new dataset with outliers removed
write.csv(KER1_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wlm_lev_bar.pdf")
barplot(hatvalues(KER1_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wlm_lev.pdf")
leveragePlots(KER1_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER1_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER1_wls_hats <- as.data.frame(hatvalues(KER1_wls))
KER1_wls_hats
# Cooks distance - 2-3 x difference from mean 
KER1_wls_cooksD <- cooks.distance(KER1_wls)
KER1_wls_influential <- KER1_wls_cooksD[(KER1_wls_cooksD > (3 * mean(KER1_wls_cooksD, na.rm = TRUE)))]
KER1_wls_influential
KER1_wls_influential_names <- names(KER1_wls_influential)
KER1_wls_outliers <- KER1_Subsample_xrf_matched[KER1_wls_influential_names,] # outliers only using of index values
KER1_wls_no_outliers <- KER1_Subsample_xrf_matched %>% anti_join(KER1_wls_outliers) # generates a new dataset with outliers removed
write.csv(KER1_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wls_lev_bar.pdf")
barplot(hatvalues(KER1_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wls_lev.pdf")
leveragePlots(KER1_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER1_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER1_wls_wt_hats <- as.data.frame(hatvalues(KER1_wls_wt))
KER1_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
KER1_wls_wt_cooksD <- cooks.distance(KER1_wls_wt)
KER1_wls_wt_influential <- KER1_wls_wt_cooksD[(KER1_wls_wt_cooksD > (3 * mean(KER1_wls_wt_cooksD, na.rm = TRUE)))]
KER1_wls_wt_influential
KER1_wls_wt_influential_names <- names(KER1_wls_wt_influential)
KER1_wls_wt_outliers <- KER1_Subsample_xrf_matched[KER1_wls_wt_influential_names,] # outliers only using of index values
KER1_wls_wt_no_outliers <- KER1_Subsample_xrf_matched %>% anti_join(KER1_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(KER1_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(KER1_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wls_wt_lev.pdf")
leveragePlots(KER1_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER1_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "KER1"
element_title <- "coh_inc"
theme_set(theme_classic(10))
KER1 <- ggplot(KER1_Subsample_xrf_matched, aes(x = coh_inc, y = Dry_mass)) +
  geom_errorbar(aes(ymin=Dry_mass-Dry_mass_err, ymax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_point(fill = "gray", color = 'gray', size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/Dry_mass_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = KER1_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = KER1_wt_wt), colour="black") + # weighted WLS regression
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
  #xlim(0.15, 0.25) +
  #ylim(10, 80) +
  ggtitle(paste("Site: ", site_title, ": Dry mass vs ", element_title))
KER1
# Define p value, OLS equation & R2 as a string to add to plot

KER1_DM_lm_p <- function(KER1_DM_lm) {
  f <- summary(KER1_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER1_DM_lm_p(KER1_DM_lm)

KER1_DM_lm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER1_DM_lm_p(KER1_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

KER1_wlm_p <- function(KER1_wlm) {
  f <- summary(KER1_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER1_wlm_p(KER1_wlm)

KER1_wlm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = KER1_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER1_wlm_p(KER1_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

KER1_wls_p <- function(KER1_wls) {
  f <- summary(KER1_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER1_wls_p(KER1_wls)

KER1_wls_eqn <- function(df){
  m <- KER1_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER1_wls_p(KER1_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

KER1_wls_wt_p <- function(KER1_wls_wt) {
  f <- summary(KER1_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER1_wls_wt_p(KER1_wls_wt)

KER1_wls_wt_eqn <- function(df){
  m <- KER1_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER1_wls_wt_p(KER1_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER1_final <- KER1 + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER1_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER1_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER1_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER1_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
KER1_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER1/FigS12_KER1_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# KER3 OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
KER3_DM_lm <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched)
summary(KER3_DM_lm)
glance(KER3_DM_lm)
model_performance(KER3_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_OLS_summary.txt")
summary(KER3_DM_lm)
glance(KER3_DM_lm)
model_performance(KER3_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_lm) # Performance package summary check for heteroscedasticity
icc(KER3_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
KER3_wlm <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2)
summary(KER3_wlm)
glance(KER3_wlm)
model_performance(KER3_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_OLS_wt_summary.txt")
summary(KER3_wlm)
glance(KER3_wlm)
model_performance(KER3_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_DM_lm) # Performance package summary check for heteroscedasticity
icc(KER3_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
KER3_model <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched) # define model
KER3_wt <- 1 / lm(abs(KER3_model$residuals) ~ KER3_model$fitted.values)$fitted.values^2 #define weights to use
KER3_wls <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched, weights=KER3_wt) #perform weighted least squares regression
# Checks
summary(KER3_wls) # summary stats
glance(KER3_wls) # summary stats including AIC
model_performance(KER3_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_WLS_summary.txt")
summary(KER3_wls)
glance(KER3_wls)
model_performance(KER3_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_wls) # Performance package summary check for heteroscedasticity
icc(KER3_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
KER3_model_wt <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched, weight = 1/Dry_mass_err^2) # define model
KER3_wt_wt <- 1 / lm(abs(KER3_model_wt$residuals) ~ KER3_model_wt$fitted.values)$fitted.values^2 #define weights to use
KER3_wls_wt <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched, weights=KER3_wt_wt) #perform weighted least squares regression
# Checks
summary(KER3_wls_wt) # summary stats
glance(KER3_wls_wt) # summary stats including AIC
model_performance(KER3_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(KER3_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_WLS_wt_summary.txt")
summary(KER3_wls_wt)
glance(KER3_wls_wt)
model_performance(KER3_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(KER3_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(KER3_wls_wt) # Performance package summary check for heteroscedasticity
icc(KER3_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER3_DM_lm_hats <- as.data.frame(hatvalues(KER3_DM_lm))
KER3_DM_lm_hats
# Cooks distance - 2-3 x difference from mean 
KER3_DM_lm_cooksD <- cooks.distance(KER3_DM_lm)
KER3_DM_lm_influential <- KER3_DM_lm_cooksD[(KER3_DM_lm_cooksD > (3 * mean(KER3_DM_lm_cooksD, na.rm = TRUE)))]
KER3_DM_lm_influential
KER3_DM_lm_influential_names <- names(KER3_DM_lm_influential)
KER3_DM_lm_outliers <- KER3_Subsample_xrf_matched[KER3_DM_lm_influential_names,] # outliers only using of index values
KER3_DM_lm_no_outliers <- KER3_Subsample_xrf_matched %>% anti_join(KER3_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(KER3_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(KER3_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_DM_lm_lev.pdf")
leveragePlots(KER3_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER3_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER3_wlm_hats <- as.data.frame(hatvalues(KER3_wlm))
KER3_wlm_hats
# Cooks distance - 2-3 x difference from mean 
KER3_wlm_cooksD <- cooks.distance(KER3_wlm)
KER3_wlm_influential <- KER3_wlm_cooksD[(KER3_wlm_cooksD > (3 * mean(KER3_wlm_cooksD, na.rm = TRUE)))]
KER3_wlm_influential
KER3_wlm_influential_names <- names(KER3_wlm_influential)
KER3_wlm_outliers <- KER3_Subsample_xrf_matched[KER3_wlm_influential_names,] # outliers only using of index values
KER3_wlm_no_outliers <- KER3_Subsample_xrf_matched %>% anti_join(KER3_wlm_outliers) # generates a new dataset with outliers removed
write.csv(KER3_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wlm_lev_bar.pdf")
barplot(hatvalues(KER3_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wlm_lev.pdf")
leveragePlots(KER3_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER3_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER3_wls_hats <- as.data.frame(hatvalues(KER3_wls))
KER3_wls_hats
# Cooks distance - 2-3 x difference from mean 
KER3_wls_cooksD <- cooks.distance(KER3_wls)
KER3_wls_influential <- KER3_wls_cooksD[(KER3_wls_cooksD > (3 * mean(KER3_wls_cooksD, na.rm = TRUE)))]
KER3_wls_influential
KER3_wls_influential_names <- names(KER3_wls_influential)
KER3_wls_outliers <- KER3_Subsample_xrf_matched[KER3_wls_influential_names,] # outliers only using of index values
KER3_wls_no_outliers <- KER3_Subsample_xrf_matched %>% anti_join(KER3_wls_outliers) # generates a new dataset with outliers removed
write.csv(KER3_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wls_lev_bar.pdf")
barplot(hatvalues(KER3_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wls_lev.pdf")
leveragePlots(KER3_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER3_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
KER3_wls_wt_hats <- as.data.frame(hatvalues(KER3_wls_wt))
KER3_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
KER3_wls_wt_cooksD <- cooks.distance(KER3_wls_wt)
KER3_wls_wt_influential <- KER3_wls_wt_cooksD[(KER3_wls_wt_cooksD > (3 * mean(KER3_wls_wt_cooksD, na.rm = TRUE)))]
KER3_wls_wt_influential
KER3_wls_wt_influential_names <- names(KER3_wls_wt_influential)
KER3_wls_wt_outliers <- KER3_Subsample_xrf_matched[KER3_wls_wt_influential_names,] # outliers only using of index values
KER3_wls_wt_no_outliers <- KER3_Subsample_xrf_matched %>% anti_join(KER3_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(KER3_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(KER3_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wls_wt_lev.pdf")
leveragePlots(KER3_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(KER3_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "KER3"
element_title <- "coh_inc"
theme_set(theme_classic(10))
KER3 <- ggplot(KER3_Subsample_xrf_matched, aes(x = coh_inc, y = Dry_mass)) +
  geom_errorbar(aes(ymin=Dry_mass-Dry_mass_err, ymax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_point(fill = "coral3", color = 'coral3', size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/Dry_mass_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = KER3_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = KER3_wt_wt), colour="black") + # weighted WLS regression
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
  #xlim(0.15, 0.25) +
  #ylim(10, 80) +
  ggtitle(paste("Site: ", site_title, ": Dry mass vs ", element_title))
KER3
# Define p value, OLS equation & R2 as a string to add to plot

KER3_DM_lm_p <- function(KER3_DM_lm) {
  f <- summary(KER3_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_DM_lm_p(KER3_DM_lm)

KER3_DM_lm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_DM_lm_p(KER3_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

KER3_wlm_p <- function(KER3_wlm) {
  f <- summary(KER3_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_wlm_p(KER3_wlm)

KER3_wlm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = KER3_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_wlm_p(KER3_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

KER3_wls_p <- function(KER3_wls) {
  f <- summary(KER3_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_wls_p(KER3_wls)

KER3_wls_eqn <- function(df){
  m <- KER3_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_wls_p(KER3_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

KER3_wls_wt_p <- function(KER3_wls_wt) {
  f <- summary(KER3_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
KER3_wls_wt_p(KER3_wls_wt)

KER3_wls_wt_eqn <- function(df){
  m <- KER3_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(KER3_wls_wt_p(KER3_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
KER3_final <- KER3 + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", KER3_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", KER3_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", KER3_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", KER3_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
KER3_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/KER3/FigS12_KER3_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# PB1 OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
PB1_DM_lm <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched)
summary(PB1_DM_lm)
glance(PB1_DM_lm)
model_performance(PB1_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(PB1_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_OLS_summary.txt")
summary(PB1_DM_lm)
glance(PB1_DM_lm)
model_performance(PB1_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(PB1_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(PB1_DM_lm) # Performance package summary check for heteroscedasticity
icc(PB1_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
PB1_wlm <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2)
summary(PB1_wlm)
glance(PB1_wlm)
model_performance(PB1_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(PB1_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_OLS_wt_summary.txt")
summary(PB1_wlm)
glance(PB1_wlm)
model_performance(PB1_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(PB1_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(PB1_DM_lm) # Performance package summary check for heteroscedasticity
icc(PB1_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
# This doesnt work because of NAs 
# PB1_model <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched) # define model
# PB1_wt <- 1 / lm(abs(PB1_model$residuals) ~ PB1_model$fitted.values)$fitted.values^2 #define weights to use
# PB1_wls <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched, weights=PB1_wt) #perform weighted least squares regression

# As for ACE dataset 
PB1_model <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched)
mf_PB1 <- model.frame(PB1_model)
var_model_PB1 <- lm(abs(residuals(PB1_model)) ~ fitted(PB1_model))
PB1_wt_internal <- 1 / (fitted(var_model_PB1)^2)
PB1_Subsample_xrf_matched$PB1_wt <- NA_real_
PB1_Subsample_xrf_matched$PB1_wt[as.numeric(rownames(mf_PB1))] <- PB1_wt_internal
length(PB1_Subsample_xrf_matched$PB1_wt) # check TRUE = same as nrow(PB1_Subsample_xrf_matched)

PB1_wls <- lm(
  Dry_mass ~ coh_inc,
  data = PB1_Subsample_xrf_matched,
  weights = PB1_wt
)

# Checks
summary(PB1_wls) # summary stats
glance(PB1_wls) # summary stats including AIC
model_performance(PB1_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(PB1_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_WLS_summary.txt")
summary(PB1_wls)
glance(PB1_wls)
model_performance(PB1_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(PB1_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(PB1_wls) # Performance package summary check for heteroscedasticity
icc(PB1_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model

# This doesnt work - use below - due to NAs 
# PB1_model_wt <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched, weight = 1/Dry_mass_err^2) # define model
# PB1_wt_wt <- 1 / lm(abs(PB1_model_wt$residuals) ~ PB1_model_wt$fitted.values)$fitted.values^2 #define weights to use
# PB1_wls_wt <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched, weights=PB1_wt_wt) #perform weighted least squares regression

# As for ACE correction - do this instead
PB1_model_wt <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched, weights = 1 / Dry_mass_err^2)
mf_PB1_wt <- model.frame(PB1_model_wt)
var_model_PB1_wt <- lm(abs(residuals(PB1_model_wt)) ~ fitted(PB1_model_wt))
PB1_wt_wt_internal <- 1 / (fitted(var_model_PB1_wt)^2)
PB1_Subsample_xrf_matched$PB1_wt_wt <- NA_real_
#Insert the internal weights into the original rows
PB1_Subsample_xrf_matched$PB1_wt_wt[
  as.numeric(rownames(mf_PB1_wt))
] <- PB1_wt_wt_internal
# check equals nrow(PB1_Subsample_xrf_matched)
length(PB1_Subsample_xrf_matched$PB1_wt_wt) # check equals nrow(PB1_Subsample_xrf_matched)

PB1_wls_wt <- lm(
  Dry_mass ~ coh_inc,
  data = PB1_Subsample_xrf_matched,
  weights = PB1_wt_wt
)

# Checks
summary(PB1_wls_wt) # summary stats
glance(PB1_wls_wt) # summary stats including AIC
model_performance(PB1_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(PB1_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_WLS_wt_summary.txt")
summary(PB1_wls_wt)
glance(PB1_wls_wt)
model_performance(PB1_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(PB1_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(PB1_wls_wt) # Performance package summary check for heteroscedasticity
icc(PB1_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
PB1_DM_lm_hats <- as.data.frame(hatvalues(PB1_DM_lm))
PB1_DM_lm_hats
# Cooks distance - 2-3 x difference from mean 
PB1_DM_lm_cooksD <- cooks.distance(PB1_DM_lm)
PB1_DM_lm_influential <- PB1_DM_lm_cooksD[(PB1_DM_lm_cooksD > (3 * mean(PB1_DM_lm_cooksD, na.rm = TRUE)))]
PB1_DM_lm_influential
PB1_DM_lm_influential_names <- names(PB1_DM_lm_influential)
PB1_DM_lm_outliers <- PB1_Subsample_xrf_matched[PB1_DM_lm_influential_names,] # outliers only using of index values
PB1_DM_lm_no_outliers <- PB1_Subsample_xrf_matched %>% anti_join(PB1_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(PB1_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(PB1_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_DM_lm_lev.pdf")
leveragePlots(PB1_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(PB1_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
PB1_wlm_hats <- as.data.frame(hatvalues(PB1_wlm))
PB1_wlm_hats
# Cooks distance - 2-3 x difference from mean 
PB1_wlm_cooksD <- cooks.distance(PB1_wlm)
PB1_wlm_influential <- PB1_wlm_cooksD[(PB1_wlm_cooksD > (3 * mean(PB1_wlm_cooksD, na.rm = TRUE)))]
PB1_wlm_influential
PB1_wlm_influential_names <- names(PB1_wlm_influential)
PB1_wlm_outliers <- PB1_Subsample_xrf_matched[PB1_wlm_influential_names,] # outliers only using of index values
PB1_wlm_no_outliers <- PB1_Subsample_xrf_matched %>% anti_join(PB1_wlm_outliers) # generates a new dataset with outliers removed
write.csv(PB1_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wlm_lev_bar.pdf")
barplot(hatvalues(PB1_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wlm_lev.pdf")
leveragePlots(PB1_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(PB1_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
PB1_wls_hats <- as.data.frame(hatvalues(PB1_wls))
PB1_wls_hats
# Cooks distance - 2-3 x difference from mean 
PB1_wls_cooksD <- cooks.distance(PB1_wls)
PB1_wls_influential <- PB1_wls_cooksD[(PB1_wls_cooksD > (3 * mean(PB1_wls_cooksD, na.rm = TRUE)))]
PB1_wls_influential
PB1_wls_influential_names <- names(PB1_wls_influential)
PB1_wls_outliers <- PB1_Subsample_xrf_matched[PB1_wls_influential_names,] # outliers only using of index values
PB1_wls_no_outliers <- PB1_Subsample_xrf_matched %>% anti_join(PB1_wls_outliers) # generates a new dataset with outliers removed
write.csv(PB1_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wls_lev_bar.pdf")
barplot(hatvalues(PB1_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wls_lev.pdf")
leveragePlots(PB1_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(PB1_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
PB1_wls_wt_hats <- as.data.frame(hatvalues(PB1_wls_wt))
PB1_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
PB1_wls_wt_cooksD <- cooks.distance(PB1_wls_wt)
PB1_wls_wt_influential <- PB1_wls_wt_cooksD[(PB1_wls_wt_cooksD > (3 * mean(PB1_wls_wt_cooksD, na.rm = TRUE)))]
PB1_wls_wt_influential
PB1_wls_wt_influential_names <- names(PB1_wls_wt_influential)
PB1_wls_wt_outliers <- PB1_Subsample_xrf_matched[PB1_wls_wt_influential_names,] # outliers only using of index values
PB1_wls_wt_no_outliers <- PB1_Subsample_xrf_matched %>% anti_join(PB1_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(PB1_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(PB1_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wls_wt_lev.pdf")
leveragePlots(PB1_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(PB1_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "PB1"
element_title <- "coh_inc"
theme_set(theme_classic(10))
PB1 <- ggplot(PB1_Subsample_xrf_matched, aes(x = coh_inc, y = Dry_mass)) +
  geom_errorbar(aes(ymin=Dry_mass-Dry_mass_err, ymax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_point(fill = "lightblue", color = 'lightblue', size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/Dry_mass_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = PB1_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = PB1_wt_wt), colour="black") + # weighted WLS regression
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
  #xlim(0.15, 0.25) +
  #ylim(10, 80) +
  ggtitle(paste("Site: ", site_title, ": Dry mass vs ", element_title))
PB1
# Define p value, OLS equation & R2 as a string to add to plot

PB1_DM_lm_p <- function(PB1_DM_lm) {
  f <- summary(PB1_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
PB1_DM_lm_p(PB1_DM_lm)

PB1_DM_lm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(PB1_DM_lm_p(PB1_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

PB1_wlm_p <- function(PB1_wlm) {
  f <- summary(PB1_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
PB1_wlm_p(PB1_wlm)

PB1_wlm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = PB1_Subsample_xrf_matched, weight = 1/(Dry_mass_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(PB1_wlm_p(PB1_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

PB1_wls_p <- function(PB1_wls) {
  f <- summary(PB1_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
PB1_wls_p(PB1_wls)

PB1_wls_eqn <- function(df){
  m <- PB1_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(PB1_wls_p(PB1_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

PB1_wls_wt_p <- function(PB1_wls_wt) {
  f <- summary(PB1_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
PB1_wls_wt_p(PB1_wls_wt)

PB1_wls_wt_eqn <- function(df){
  m <- PB1_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(PB1_wls_wt_p(PB1_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
PB1_final <- PB1 + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", PB1_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", PB1_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", PB1_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", PB1_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
PB1_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/PB1/FigS12_PB1_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")



# 5) LOI & Corg vs inc/coh regression ------------------------------------------
# Remove PB1 from ACE matched dataset - no PB1 LOI550 data available -----------


ACE_Subsample_xrf_matched_no_PB1 <- bind_rows(BI10_Subsample_xrf_matched,
                                      HER42PB_Subsample_xrf_matched,
                                      KER1_Subsample_xrf_matched,
                                      KER3_Subsample_xrf_matched)
ACE_Subsample_xrf_matched_no_PB1

write.csv(ACE_Subsample_xrf_matched_no_PB1,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/ACE/ACE_Subsample_xrf_matched_no_PB1.csv", 
          row.names = FALSE)

# ACE LOI550 OLS & WLS ---------------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
ACE_LOI_lm <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1)
summary(ACE_LOI_lm)
glance(ACE_LOI_lm)
model_performance(ACE_LOI_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_LOI_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_OLS_summary.txt")
summary(ACE_LOI_lm)
glance(ACE_LOI_lm)
model_performance(ACE_LOI_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_LOI_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_LOI_lm) # Performance package summary check for heteroscedasticity
icc(ACE_LOI_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_wlm <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weight = 1/(LOI550_err)^2)
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_OLS_wt_summary.txt")
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_LOI_lm) # Performance package summary check for heteroscedasticity
icc(ACE_LOI_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
ACE_model <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1) # define model
ACE_wt <- 1 / lm(abs(ACE_model$residuals) ~ ACE_model$fitted.values)$fitted.values^2 #define weights to use
ACE_wls <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weights=ACE_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls) # summary stats
glance(ACE_wls) # summary stats including AIC
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_WLS_summary.txt")
summary(ACE_wls)
glance(ACE_wls)
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls) # Performance package summary check for heteroscedasticity
icc(ACE_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
ACE_model_wt <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weight = 1/LOI550_err^2) # define model
ACE_wt_wt <- 1 / lm(abs(ACE_model_wt$residuals) ~ ACE_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_wls_wt <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weights=ACE_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls_wt) # summary stats
glance(ACE_wls_wt) # summary stats including AIC
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_WLS_wt_summary.txt")
summary(ACE_wls_wt)
glance(ACE_wls_wt)
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_LOI_lm_hats <- as.data.frame(hatvalues(ACE_LOI_lm))
ACE_LOI_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_LOI_lm_cooksD <- cooks.distance(ACE_LOI_lm)
ACE_LOI_lm_influential <- ACE_LOI_lm_cooksD[(ACE_LOI_lm_cooksD > (3 * mean(ACE_LOI_lm_cooksD, na.rm = TRUE)))]
ACE_LOI_lm_influential
ACE_LOI_lm_influential_names <- names(ACE_LOI_lm_influential)
ACE_LOI_lm_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_LOI_lm_influential_names,] # outliers only using of index values
ACE_LOI_lm_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_LOI_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_LOI_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_LOI_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_DM_lm_lev.pdf")
leveragePlots(ACE_LOI_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_LOI_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wlm_hats <- as.data.frame(hatvalues(ACE_wlm))
ACE_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wlm_cooksD <- cooks.distance(ACE_wlm)
ACE_wlm_influential <- ACE_wlm_cooksD[(ACE_wlm_cooksD > (3 * mean(ACE_wlm_cooksD, na.rm = TRUE)))]
ACE_wlm_influential
ACE_wlm_influential_names <- names(ACE_wlm_influential)
ACE_wlm_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_wlm_influential_names,] # outliers only using of index values
ACE_wlm_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wlm_lev.pdf")
leveragePlots(ACE_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_hats <- as.data.frame(hatvalues(ACE_wls))
ACE_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_cooksD <- cooks.distance(ACE_wls)
ACE_wls_influential <- ACE_wls_cooksD[(ACE_wls_cooksD > (3 * mean(ACE_wls_cooksD, na.rm = TRUE)))]
ACE_wls_influential
ACE_wls_influential_names <- names(ACE_wls_influential)
ACE_wls_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_wls_influential_names,] # outliers only using of index values
ACE_wls_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wls_lev_bar.pdf")
barplot(hatvalues(ACE_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wls_lev.pdf")
leveragePlots(ACE_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_wt_hats <- as.data.frame(hatvalues(ACE_wls_wt))
ACE_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_wt_cooksD <- cooks.distance(ACE_wls_wt)
ACE_wls_wt_influential <- ACE_wls_wt_cooksD[(ACE_wls_wt_cooksD > (3 * mean(ACE_wls_wt_cooksD, na.rm = TRUE)))]
ACE_wls_wt_influential
ACE_wls_wt_influential_names <- names(ACE_wls_wt_influential)
ACE_wls_wt_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_wls_wt_influential_names,] # outliers only using of index values
ACE_wls_wt_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wls_wt_lev.pdf")
leveragePlots(ACE_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/LOI550_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "ACE"
element_title <- "inc_coh"
theme_set(theme_classic(10))
ACE <- ggplot(ACE_Subsample_xrf_matched_no_PB1, aes(x = inc_coh, y = LOI550)) +
  geom_errorbar(aes(ymin=LOI550-LOI550_err, ymax=LOI550+LOI550_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=inc_coh-inc_coh_sd, xmax=inc_coh+inc_coh_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_point(fill = "Site", color = "Site", size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/LOI550_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = ACE_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = ACE_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, "[XRF-CS]") , y = paste0("LOI 550C (%)")) +
  scale_y_continuous(labels = scales::comma) +
  xlim(3, 8) +
  ylim(0, 100) +
  ggtitle(paste("Site: ", site_title, ": LOI550 vs ", element_title))
ACE
# Define p value, OLS equation & R2 as a string to add to plot

ACE_LOI_lm_p <- function(ACE_LOI_lm) {
  f <- summary(ACE_LOI_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_LOI_lm_p(ACE_LOI_lm)

ACE_LOI_lm_eqn <- function(df){
  m <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_LOI_lm_p(ACE_LOI_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_wlm_p <- function(ACE_wlm) {
  f <- summary(ACE_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wlm_p(ACE_wlm)

ACE_wlm_eqn <- function(df){
  m <- lm(LOI550 ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weight = 1/(LOI550_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wlm_p(ACE_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_wls_p <- function(ACE_wls) {
  f <- summary(ACE_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_p(ACE_wls)

ACE_wls_eqn <- function(df){
  m <- ACE_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_p(ACE_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_wls_wt_p <- function(ACE_wls_wt) {
  f <- summary(ACE_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_wt_p(ACE_wls_wt)

ACE_wls_wt_eqn <- function(df){
  m <- ACE_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_wt_p(ACE_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_final <- ACE + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_LOI_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/FigS12_ACE_LOI_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ACE %Carbon OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
ACE_C_lm <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1)
summary(ACE_C_lm)
glance(ACE_C_lm)
model_performance(ACE_C_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_C_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_OLS_summary.txt")
summary(ACE_C_lm)
glance(ACE_C_lm)
model_performance(ACE_C_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_C_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_C_lm) # Performance package summary check for heteroscedasticity
icc(ACE_C_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_wlm <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weight = 1/(C_org_pc_err)^2)
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_OLS_wt_summary.txt")
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_C_lm) # Performance package summary check for heteroscedasticity
icc(ACE_C_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
ACE_model <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1) # define model
ACE_wt <- 1 / lm(abs(ACE_model$residuals) ~ ACE_model$fitted.values)$fitted.values^2 #define weights to use
ACE_wls <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weights=ACE_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls) # summary stats
glance(ACE_wls) # summary stats including AIC
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_WLS_summary.txt")
summary(ACE_wls)
glance(ACE_wls)
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls) # Performance package summary check for heteroscedasticity
icc(ACE_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
ACE_model_wt <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weight = 1/C_org_pc_err^2) # define model
ACE_wt_wt <- 1 / lm(abs(ACE_model_wt$residuals) ~ ACE_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_wls_wt <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weights=ACE_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls_wt) # summary stats
glance(ACE_wls_wt) # summary stats including AIC
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_WLS_wt_summary.txt")
summary(ACE_wls_wt)
glance(ACE_wls_wt)
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_C_lm_hats <- as.data.frame(hatvalues(ACE_C_lm))
ACE_C_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_C_lm_cooksD <- cooks.distance(ACE_C_lm)
ACE_C_lm_influential <- ACE_C_lm_cooksD[(ACE_C_lm_cooksD > (3 * mean(ACE_C_lm_cooksD, na.rm = TRUE)))]
ACE_C_lm_influential
ACE_C_lm_influential_names <- names(ACE_C_lm_influential)
ACE_C_lm_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_C_lm_influential_names,] # outliers only using of index values
ACE_C_lm_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_C_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_C_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_C_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_lm_lev.pdf")
leveragePlots(ACE_C_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_C_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wlm_hats <- as.data.frame(hatvalues(ACE_wlm))
ACE_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wlm_cooksD <- cooks.distance(ACE_wlm)
ACE_wlm_influential <- ACE_wlm_cooksD[(ACE_wlm_cooksD > (3 * mean(ACE_wlm_cooksD, na.rm = TRUE)))]
ACE_wlm_influential
ACE_wlm_influential_names <- names(ACE_wlm_influential)
ACE_wlm_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_wlm_influential_names,] # outliers only using of index values
ACE_wlm_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wlm_lev.pdf")
leveragePlots(ACE_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_hats <- as.data.frame(hatvalues(ACE_wls))
ACE_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_cooksD <- cooks.distance(ACE_wls)
ACE_wls_influential <- ACE_wls_cooksD[(ACE_wls_cooksD > (3 * mean(ACE_wls_cooksD, na.rm = TRUE)))]
ACE_wls_influential
ACE_wls_influential_names <- names(ACE_wls_influential)
ACE_wls_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_wls_influential_names,] # outliers only using of index values
ACE_wls_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wls_lev_bar.pdf")
barplot(hatvalues(ACE_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wls_lev.pdf")
leveragePlots(ACE_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_wt_hats <- as.data.frame(hatvalues(ACE_wls_wt))
ACE_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_wt_cooksD <- cooks.distance(ACE_wls_wt)
ACE_wls_wt_influential <- ACE_wls_wt_cooksD[(ACE_wls_wt_cooksD > (3 * mean(ACE_wls_wt_cooksD, na.rm = TRUE)))]
ACE_wls_wt_influential
ACE_wls_wt_influential_names <- names(ACE_wls_wt_influential)
ACE_wls_wt_outliers <- ACE_Subsample_xrf_matched_no_PB1[ACE_wls_wt_influential_names,] # outliers only using of index values
ACE_wls_wt_no_outliers <- ACE_Subsample_xrf_matched_no_PB1 %>% anti_join(ACE_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wls_wt_lev.pdf")
leveragePlots(ACE_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "ACE"
element_title <- "inc_coh"
theme_set(theme_classic(10))
ACE <- ggplot(ACE_Subsample_xrf_matched_no_PB1, aes(x = inc_coh, y = C_org_pc)) +
  geom_errorbar(aes(ymin=C_org_pc-C_org_pc_err, ymax=C_org_pc+C_org_pc_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=inc_coh-inc_coh_sd, xmax=inc_coh+inc_coh_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_point(fill = "Site", color = "Site", size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/C_org_pc_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = ACE_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = ACE_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, "[XRF-CS]") , y = paste0("Corg (%)")) +
  scale_y_continuous(labels = scales::comma) +
  xlim(3, 8) +
  #ylim(0, 100) +
  ggtitle(paste("Site: ", site_title, ": Organic Carbon vs ", element_title))
ACE
# Define p value, OLS equation & R2 as a string to add to plot

ACE_C_lm_p <- function(ACE_C_lm) {
  f <- summary(ACE_C_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_C_lm_p(ACE_C_lm)

ACE_C_lm_eqn <- function(df){
  m <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_C_lm_p(ACE_C_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_wlm_p <- function(ACE_wlm) {
  f <- summary(ACE_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wlm_p(ACE_wlm)

ACE_wlm_eqn <- function(df){
  m <- lm(C_org_pc ~ inc_coh, data = ACE_Subsample_xrf_matched_no_PB1, weight = 1/(C_org_pc_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wlm_p(ACE_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_wls_p <- function(ACE_wls) {
  f <- summary(ACE_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_p(ACE_wls)

ACE_wls_eqn <- function(df){
  m <- ACE_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_p(ACE_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_wls_wt_p <- function(ACE_wls_wt) {
  f <- summary(ACE_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_wt_p(ACE_wls_wt)

ACE_wls_wt_eqn <- function(df){
  m <- ACE_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_wt_p(ACE_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_final <- ACE + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_C_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/FigS12_ACE_Corg_pc_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# 6) Dry mass & Corg regression ------------------------------------------
# Remove PB1 from ACE matched dataset - no PB1 LOI550 data available -----------
ACE_Subsample_no_PB1 <- ACE_Subsample %>% 
  filter(!Site == "PB1") %>% 
  print()
write.csv(ACE_Subsample_no_PB1,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/ACE/ACE_Subsample_no_PB1.csv", 
          row.names = FALSE)

# ACE no PB1 - %Carbon vs DM OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
ACE_C_lm <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1)
summary(ACE_C_lm)
glance(ACE_C_lm)
model_performance(ACE_C_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_C_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_OLS_summary.txt")
summary(ACE_C_lm)
glance(ACE_C_lm)
model_performance(ACE_C_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_C_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_C_lm) # Performance package summary check for heteroscedasticity
icc(ACE_C_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_wlm <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1, weight = 1/(C_org_pc_err)^2)
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_OLS_wt_summary.txt")
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_C_lm) # Performance package summary check for heteroscedasticity
icc(ACE_C_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
ACE_model <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1) # define model
ACE_wt <- 1 / lm(abs(ACE_model$residuals) ~ ACE_model$fitted.values)$fitted.values^2 #define weights to use
ACE_wls <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1, weights=ACE_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls) # summary stats
glance(ACE_wls) # summary stats including AIC
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_WLS_summary.txt")
summary(ACE_wls)
glance(ACE_wls)
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls) # Performance package summary check for heteroscedasticity
icc(ACE_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
ACE_model_wt <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1, weight = 1/C_org_pc_err^2) # define model
ACE_wt_wt <- 1 / lm(abs(ACE_model_wt$residuals) ~ ACE_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_wls_wt <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1, weights=ACE_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls_wt) # summary stats
glance(ACE_wls_wt) # summary stats including AIC
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_WLS_wt_summary.txt")
summary(ACE_wls_wt)
glance(ACE_wls_wt)
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Leverage & Cooks distance
# 1) Unweighted OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_C_lm_hats <- as.data.frame(hatvalues(ACE_C_lm))
ACE_C_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_C_lm_cooksD <- cooks.distance(ACE_C_lm)
ACE_C_lm_influential <- ACE_C_lm_cooksD[(ACE_C_lm_cooksD > (3 * mean(ACE_C_lm_cooksD, na.rm = TRUE)))]
ACE_C_lm_influential
ACE_C_lm_influential_names <- names(ACE_C_lm_influential)
ACE_C_lm_outliers <- ACE_Subsample_no_PB1[ACE_C_lm_influential_names,] # outliers only using of index values
ACE_C_lm_no_outliers <- ACE_Subsample_no_PB1 %>% anti_join(ACE_C_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_C_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_C_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_lm_lev.pdf")
leveragePlots(ACE_C_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_C_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wlm_hats <- as.data.frame(hatvalues(ACE_wlm))
ACE_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wlm_cooksD <- cooks.distance(ACE_wlm)
ACE_wlm_influential <- ACE_wlm_cooksD[(ACE_wlm_cooksD > (3 * mean(ACE_wlm_cooksD, na.rm = TRUE)))]
ACE_wlm_influential
ACE_wlm_influential_names <- names(ACE_wlm_influential)
ACE_wlm_outliers <- ACE_Subsample_no_PB1[ACE_wlm_influential_names,] # outliers only using of index values
ACE_wlm_no_outliers <- ACE_Subsample_no_PB1 %>% anti_join(ACE_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wlm_lev.pdf")
leveragePlots(ACE_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_hats <- as.data.frame(hatvalues(ACE_wls))
ACE_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_cooksD <- cooks.distance(ACE_wls)
ACE_wls_influential <- ACE_wls_cooksD[(ACE_wls_cooksD > (3 * mean(ACE_wls_cooksD, na.rm = TRUE)))]
ACE_wls_influential
ACE_wls_influential_names <- names(ACE_wls_influential)
ACE_wls_outliers <- ACE_Subsample_no_PB1[ACE_wls_influential_names,] # outliers only using of index values
ACE_wls_no_outliers <- ACE_Subsample_no_PB1 %>% anti_join(ACE_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wls_lev_bar.pdf")
barplot(hatvalues(ACE_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wls_lev.pdf")
leveragePlots(ACE_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_wt_hats <- as.data.frame(hatvalues(ACE_wls_wt))
ACE_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_wt_cooksD <- cooks.distance(ACE_wls_wt)
ACE_wls_wt_influential <- ACE_wls_wt_cooksD[(ACE_wls_wt_cooksD > (3 * mean(ACE_wls_wt_cooksD, na.rm = TRUE)))]
ACE_wls_wt_influential
ACE_wls_wt_influential_names <- names(ACE_wls_wt_influential)
ACE_wls_wt_outliers <- ACE_Subsample_no_PB1[ACE_wls_wt_influential_names,] # outliers only using of index values
ACE_wls_wt_no_outliers <- ACE_Subsample_no_PB1 %>% anti_join(ACE_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wls_wt_lev.pdf")
leveragePlots(ACE_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/Carbon_DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "ACE"
element_title <- "Dry_mass"
theme_set(theme_classic(10))
ACE <- ggplot(ACE_Subsample_no_PB1, aes(x = Dry_mass, y = C_org_pc)) +
  geom_errorbar(aes(ymin=C_org_pc-C_org_pc_err, ymax=C_org_pc+C_org_pc_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Dry_mass-Dry_mass_err, xmax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_point(fill = "Site", color = "Site", size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/C_org_pc_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = ACE_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = ACE_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, "[XRF-CS]") , y = paste0("Corg (%)")) +
  scale_y_continuous(labels = scales::comma) +
  xlim(3, 8) +
  #ylim(0, 100) +
  ggtitle(paste("Site: ", site_title, ": Organic Carbon vs ", element_title))
ACE
# Define p value, OLS equation & R2 as a string to add to plot

ACE_C_lm_p <- function(ACE_C_lm) {
  f <- summary(ACE_C_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_C_lm_p(ACE_C_lm)

ACE_C_lm_eqn <- function(df){
  m <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_C_lm_p(ACE_C_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_wlm_p <- function(ACE_wlm) {
  f <- summary(ACE_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wlm_p(ACE_wlm)

ACE_wlm_eqn <- function(df){
  m <- lm(C_org_pc ~ Dry_mass, data = ACE_Subsample_no_PB1, weight = 1/(C_org_pc_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wlm_p(ACE_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_wls_p <- function(ACE_wls) {
  f <- summary(ACE_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_p(ACE_wls)

ACE_wls_eqn <- function(df){
  m <- ACE_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_p(ACE_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_wls_wt_p <- function(ACE_wls_wt) {
  f <- summary(ACE_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_wt_p(ACE_wls_wt)

ACE_wls_wt_eqn <- function(df){
  m <- ACE_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_wt_p(ACE_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_final <- ACE + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_C_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/FigS12_ACE_Carbon_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")



# END

# END --------------------------------------------------------------------------

