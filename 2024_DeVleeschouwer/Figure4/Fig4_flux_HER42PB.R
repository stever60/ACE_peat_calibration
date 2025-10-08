# Figure 4 - Ti flux calculations 

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

# Load libraries ----------------------------------------------------------
packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'purrr', 'readr', 'ggpubr', 'patchwork',
              'gridExtra', 'cowplot', 'vegan', 'rioja', 'ellipse', 'factoextra',
              'reshape2', 'GGally', 'ggsci', 'ggdendro', 'dendextend', 'dynamicTreeCut',
              'colorspace', 'cluster', 'magrittr', 'mgcv', 'gtable', 'repr',
              'bestNormalize','sjmisc', 'chemometrics', 'compositions', 
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis' )
lapply(packages, library, character.only=TRUE)



# Define parameters for plotting--------------------------------------------------

# XRF-CS
# elements defined by ITRAX acf and matched to Francois ICPMS element list
# XRF-CS acf elements matched to Francois ICPMS elements
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

# ICPMS
# elements defined by Francois & by ITRAX acf
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", 
                      "Co_ICP", "Ni_ICP", "Cu_ICP", "Zn_ICP", "As_ICP", "Rb_ICP", 
                      "Sr_ICP", "Zr_ICP", "Pb_ICP", "Dry_mass")
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", 
                      "Ni_ICP", "Cu_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", 
                         "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# Matched XRF and ICPMS elements
xrf_icp_elements <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                       "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP", "Mo_coh")
xrf_icp_elements1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                      "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                      "Zr", "Zr_ICP", "Mo_inc",  "Mo_coh") # Mo_inc included

# MSCL
mscl_param  <- c("SH20_mean_age", "SH20_mean_95CI", "accum_rate", "accum_rate_err",
                 "Den1_SAT", "Den1_SAT_err", "MS1_SAT", "DCMS1_SAT", "Impedance_SAT", 
                 "Fract_Porosity_SAT", "Resistivity_SAT", "WMAR_SAT", "WMAR_SAT_err")

# Subsample parameters
subsample_param <- c("Water_Content", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density","Dry_density",	"Dry_density_err", 
                     "DMAR", "DMAR_err", "DMAR_err_pc")

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density")

# Ti ppm & flux conversions ----------------------------------------------------

# Import ACE Ti XRF-CS log_inc data a
ACE_xrf_log_inc <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_log_inc.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>% 
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc) %>%
  filter(qc == "TRUE") %>% 
  filter(validity == "1") %>% 
  mutate(across(acf_icp_Elements_key, ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_log_inc

# Convert to ppm with +/- RMSE error from PLS ----------------------------------




# Flux calculations 
library(dplyr)
library(purrr)

# Import ACE-SHW database and create edited ACE database 
# Import ACE-SHW database and create edited ACE database 
ACE_MSCL <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_SHW_MSCL_Composite.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>%
  select(c(Location:SH20_median_age, mscl_param))
write.csv(ACE_MSCL,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_MSCL.csv", row.names = FALSE)

HER42PB_MSCL <- ACE_MSCL %>%
  filter(Site=="HER42PB")

# Interpolate 0.5 cm intervals to 0.1 cm intervals for matching & calibrating wet GRD vs MSCL dry density 


# Example dataframe with multiple parameters at 0.5 cm intervals

set.seed(123)

depth <- seq(0, 5, by = 0.5)

param1 <- sin(depth)

param2 <- cos(depth)

param3 <- depth^2

df <- data.frame(depth, param1, param2, param3)


# Define new depth sequence at 0.1 cm intervals

new_depth <- seq(min(df$depth), max(df$depth), by = 0.1)


# Interpolate each parameter column

# We'll skip the 'depth' column in interpolation

interpolated_df <- data.frame(depth = new_depth)


# Apply linear interpolation (approx) to each parameter

for (colname in names(df)[-1]) {
  
  interpolated_df[[colname]] <- approx(x = df$depth, y = df[[colname]], xout = new_depth)$y
  
}


# Result: 'interpolated_df' contains all parameters at 0.1 cm intervals

print(head(interpolated_df))




# Approach (Tidyverse style)



library(dplyr)

library(purrr)



# Interpolation function for one column

interp_fun <- function(y, depth, new_depth) {
  
  approx(x = depth, y = y, xout = new_depth)$y
  
}



# Apply to all columns except 'depth'

interpolated_df <- tibble(depth = new_depth)

params <- names(df)[-1]



interpolated_df <- interpolated_df %>%
  
  bind_cols(map_dfc(params, ~ interp_fun(df[[.x]], df$depth, new_depth) %>%
                      
                      setNames(.x)))




print(head(interpolated_df))






# BI10 ICPMS & Subsample matching ---------------------------------------

# Select BI10 site data for both datasets
ACE_ICPMS_BI10 <- ACE_ICPMS %>% 
  filter(Site == "BI10") %>%
  print()

BI10_Subsample <- ACE_Subsample %>% 
  filter(Site == "BI10") %>% 
  print()

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
BI10_Subsample_ICPMS_matched <- itrax_reduce(ACE_ICPMS_BI10, names = BI10_Subsample$sample,
                                             breaks_lower = BI10_Subsample$top,
                                             breaks_upper = BI10_Subsample$bottom,
                                             fun = mean,
                                             edges = c(">=", "<=")) %>%
  select(resample_names, P_ICP:Th_ICP_sd) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., BI10_Subsample, by = "sample") %>%
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section) %>% 
  na.omit()
BI10_Subsample_ICPMS_matched


write.csv(BI10_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/BI10/BI10_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# HER42PB ICPMS & Subsample matching ---------------------------------------

# Select HER42PB site data for both datasets
ACE_ICPMS_HER42PB <- ACE_ICPMS %>% 
  filter(Site == "HER42PB") %>%
  print()

HER42PB_Subsample <- ACE_Subsample %>% 
  filter(Site == "HER42PB") %>% 
  print()

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
HER42PB_Subsample_ICPMS_matched <- itrax_reduce(ACE_ICPMS_HER42PB, names = HER42PB_Subsample$sample,
                                                breaks_lower = HER42PB_Subsample$top,
                                                breaks_upper = HER42PB_Subsample$bottom,
                                                fun = mean,
                                                edges = c(">=", "<=")) %>%
  select(resample_names, P_ICP:Th_ICP_sd) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., HER42PB_Subsample, by = "sample") %>%
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section) %>% 
  na.omit()
HER42PB_Subsample_ICPMS_matched
write.csv(HER42PB_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# KER1 ICPMS & Subsample matching ---------------------------------------

# Select KER1 site data for both datasets
ACE_ICPMS_KER1 <- ACE_ICPMS %>% 
  filter(Site == "KER1") %>%
  print()

KER1_Subsample <- ACE_Subsample %>% 
  filter(Site == "KER1") %>% 
  print()

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
KER1_Subsample_ICPMS_matched <- itrax_reduce(ACE_ICPMS_KER1, names = KER1_Subsample$sample,
                                             breaks_lower = KER1_Subsample$top,
                                             breaks_upper = KER1_Subsample$bottom,
                                             fun = mean,
                                             edges = c(">=", "<=")) %>%
  select(resample_names, P_ICP:Th_ICP_sd) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., KER1_Subsample, by = "sample") %>%
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section) %>% 
  na.omit()
KER1_Subsample_ICPMS_matched
write.csv(KER1_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/KER1/KER1_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# KER3 ICPMS & Subsample matching ---------------------------------------

# Select KER3 site data for both datasets
ACE_ICPMS_KER3 <- ACE_ICPMS %>% 
  filter(Site == "KER3") %>%
  print()

KER3_Subsample <- ACE_Subsample %>% 
  filter(Site == "KER3") %>% 
  print()

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
KER3_Subsample_ICPMS_matched <- itrax_reduce(ACE_ICPMS_KER3, names = KER3_Subsample$sample,
                                             breaks_lower = KER3_Subsample$top,
                                             breaks_upper = KER3_Subsample$bottom,
                                             fun = mean,
                                             edges = c(">=", "<=")) %>%
  select(resample_names, P_ICP:Th_ICP_sd) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., KER3_Subsample, by = "sample") %>%
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section) %>% 
  na.omit()
KER3_Subsample_ICPMS_matched
write.csv(KER3_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/KER3/KER3_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)


# PB1 ICPMS & Subsample matching ---------------------------------------

# Select PB1 site data for both datasets
ACE_ICPMS_PB1 <- ACE_ICPMS %>% 
  filter(Site == "PB1") %>%
  print()

PB1_Subsample <- ACE_Subsample %>% 
  filter(Site == "PB1") %>% 
  print()

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
PB1_Subsample_ICPMS_matched <- itrax_reduce(ACE_ICPMS_PB1, names = PB1_Subsample$sample,
                                            breaks_lower = PB1_Subsample$top,
                                            breaks_upper = PB1_Subsample$bottom,
                                            fun = mean,
                                            edges = c(">=", "<=")) %>%
  select(resample_names, P_ICP:Th_ICP_sd) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., PB1_Subsample, by = "sample") %>%
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section) %>% 
  na.omit()
PB1_Subsample_ICPMS_matched
write.csv(PB1_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/PB1/PB1_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# Create ACE Subsample-ICPMS matched dataset ----------------------------

ACE_Subsample_ICPMS_matched<- bind_rows(BI10_Subsample_ICPMS_matched,
                                        HER42PB_Subsample_ICPMS_matched,
                                        KER1_Subsample_ICPMS_matched,
                                        KER3_Subsample_ICPMS_matched,
                                        PB1_Subsample_ICPMS_matched)
ACE_Subsample_ICPMS_matched
tail(ACE_Subsample_ICPMS_matched)
write.csv(ACE_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)




