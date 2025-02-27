# ITRAX-ICPMS Error analysis - cps

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

# key elements to simplify plotting 
xrf_icp_Elements_key <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP",
                          "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "dry_mass_pc")

# key elements_reduced to simplify plotting further  
xrf_icp_Elements_key_reduced <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")



# Import matched itrax-ICPMS datafiles from each site -------------------

ACE_itrax_BI10 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/BI10/cps/BI10_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_HER42PB <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/HER42PB/cps/HER42PB_xrf_icp_matched_cps.csv", 
                              col_names = TRUE, skip = 0)
ACE_itrax_KER1 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER1/cps/KER1_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_KER3 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/cps/KER3_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_PB1 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/PB1/cps/PB1_xrf_icp_matched_cps.csv", 
                          col_names = TRUE, skip = 0)
ACE_itrax_POB4 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/POB4/cps/POB4_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)

# Combine matched output from each site
ACE_xrf_icp_matched_input <- bind_rows(ACE_itrax_BI10, 
                                 ACE_itrax_HER42PB, 
                                 ACE_itrax_KER1, 
                                 ACE_itrax_KER3, 
                                 ACE_itrax_PB1, 
                                 ACE_itrax_POB4) %>% 
  print()

write.csv(ACE_xrf_icp_matched_input,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/cps/ACE_xrf_icp_matched_cps.csv", row.names = FALSE)

# Or import pre-made ACE cps matched file ----------------------------------------------

ACE_xrf_icp_matched_input <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/cps/ACE_xrf_icp_matched_cps.csv")
is.na(ACE_xrf_icp_matched_input)<-sapply(ACE_xrf_icp_matched_input, is.infinite) # replace any infinite values with NA

# Replace zeros with half minimum value for each element 
ACE_xrf_icp_matched_0replace <- ACE_xrf_icp_matched_input %>% 
  filter(!Site =="POB4") %>% #remove POB4 data
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
write.csv(ACE_xrf_icp_matched_0replace,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Matching_err/cps/ACE_xrf_icp_matched_cps_zero_replace.csv", row.names = FALSE)

# Calculate standard error - to test 0.05% signifcance null hypothesis
n_SE = 274^0.5
ACE_xrf_matched_errors <- ACE_xrf_icp_matched_0replace %>% 
  mutate(K_sd_pc = ((K_sd/K)*100)) %>%
  mutate(K_se_pc = ((K_sd/K) / n_SE)*100) %>%
  relocate(c(K_sd,K_sd_pc, K_se_pc), .after = K) %>%
  mutate(Ca_sd_pc = ((Ca_sd/Ca)*100)) %>%
  mutate(Ca_se_pc = ((Ca_sd/Ca) / n_SE)*100) %>%
  relocate(c(Ca_sd,Ca_sd_pc,Ca_se_pc), .after = Ca) %>%
  mutate(Ti_sd_pc = ((Ti_sd/Ti)*100)) %>%
  mutate(Ti_se_pc = ((Ti_sd/Ti) / n_SE)*100) %>%
  relocate(c(Ti_sd,Ti_sd_pc,Ti_se_pc), .after = Ti) %>%
  mutate(Mn_sd_pc = ((Mn_sd/Mn)*100)) %>%
  mutate(Mn_SE = ((Mn_sd/Mn) / n_SE)*100) %>%
  relocate(c(Mn_sd,Mn_sd_pc,Mn_SE), .after = Mn) %>%
  mutate(Fe_sd_pc = ((Fe_sd/Fe)*100)) %>%
  mutate(Fe_se_pc = ((Fe_sd/Fe) / n_SE)*100) %>%
  relocate(c(Fe_sd,Fe_sd_pc,Fe_se_pc), .after = Fe) %>%
  mutate(Co_sd_pc = ((Co_sd/Co)*100)) %>%
  mutate(Co_se_pc = ((Co_sd/Co) / n_SE)*100) %>%
  relocate(c(Co_sd,Co_sd_pc,Co_se_pc), .after = Co) %>%
  mutate(Ni_sd_pc = ((Ni_sd/Ni)*100)) %>%
  mutate(Ni_se_pc = ((Ni_sd/Ni) / n_SE)*100) %>%
  relocate(c(Ni_sd,Ni_sd_pc,Ni_se_pc), .after = Ni) %>%
  mutate(Cu_sd_pc = ((Cu_sd/Cu)*100)) %>%
  mutate(Cu_se_pc = ((Cu_sd/Cu) / n_SE)*100) %>%
  relocate(c(Cu_sd,Cu_sd_pc,Cu_se_pc), .after = Cu) %>%
  mutate(Zn_sd_pc = ((Zn_sd/Zn)*100)) %>%
  mutate(Zn_SE = ((Zn_sd/Zn) / n_SE)*100) %>%
  relocate(c(Zn_sd,Zn_sd_pc,Zn_SE), .after = Zn) %>%
  mutate(Rb_sd_pc = ((Rb_sd/Rb)*100)) %>%
  mutate(Rb_se_pc = ((Rb_sd/Rb) / n_SE)*100) %>%
  relocate(c(Rb_sd,Rb_sd_pc,Rb_se_pc), .after = Rb) %>%
  mutate(Sr_sd_pc = ((Sr_sd/Sr)*100)) %>%
  mutate(Sr_se_pc = ((Sr_sd/Sr) / n_SE)*100) %>%
  relocate(c(Sr_sd,Sr_sd_pc,Sr_se_pc), .after = Sr) %>%
  mutate(Zr_sd_pc = ((Zr_sd/Zr)*100)) %>%
  mutate(Zr_se_pc = ((Zr_sd/Zr) / n_SE)*100) %>%
  relocate(c(Zr_sd,Zr_sd_pc,Zr_se_pc), .after = Zr) %>%
  mutate(Mo_inc_sd_pc = ((Mo_inc_sd/Mo_inc)*100)) %>%
  mutate(Mo_inc_se_pc = ((Mo_inc_sd/Mo_inc) / n_SE)*100) %>%
  relocate(c(Mo_inc_sd,Mo_inc_sd_pc,Mo_inc_se_pc), .after = Mo_inc) %>%
  mutate(Mo_coh_sd_pc = ((Mo_inc_sd/Mo_inc)*100)) %>%
  mutate(Mo_coh_se_pc = ((Mo_coh_sd/Mo_coh) / n_SE)*100) %>%
  relocate(c(Mo_coh_sd,Mo_coh_sd_pc,Mo_coh_se_pc), .after = Mo_coh) %>%
# ICPMS data - change errors supplied by FDV originally to CRM % values in Table 2 
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
  relocate(dry_mass_err, .after = dry_mass_pc) %>% 
  select(!(K_ICP_sd_FDV2023:Zr_ICP_sd_FDV2023)) #remove old errors 
write.csv(ACE_xrf_matched_errors,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Matching_err/cps/ACE_matched_xrf_icp_cps_errors.csv", row.names = FALSE)

# Summary stats for cps and errors - all elements and scatter parameters
ACE_matched_stats <- ACE_xrf_matched_errors %>%
  select(kcps:dry_mass_err) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(ACE_matched_stats)
write.csv(ACE_matched_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Matching_err/cps/ACE_matched_xrf_icp_cps_errors_stats.csv", row.names = FALSE)


