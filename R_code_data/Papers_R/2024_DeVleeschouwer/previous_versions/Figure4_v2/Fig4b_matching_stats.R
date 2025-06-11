# Figure 4 - inc/coh vs %C matching and stats

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
packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'readr', 'ggpubr', 'patchwork',
              'gridExtra', 'cowplot', 'vegan', 'rioja', 'ellipse', 'factoextra',
              'reshape2', 'GGally', 'ggsci', 'ggdendro', 'dendextend', 'dynamicTreeCut',
              'colorspace', 'cluster', 'magrittr', 'mgcv', 'gtable', 'repr',
              'bestNormalize','sjmisc', 'chemometrics', 'compositions', 
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis', 'itraxR' )
lapply(packages, library, character.only=TRUE)

# Import subsample dataset ----------------------------------------------------------
ACE_Subsample <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_Subsample.csv") %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all inf to NA
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  select(c(Location:Section, Sample, strat_depth_top: accrate_err, Water_Content_pc: DMAR_err_pc)) %>% 
  select(-Ash_Mass) %>% 
  rename(sample = Sample, top = strat_depth_top, bottom = strat_depth_bottom, 
         mp = strat_depth_mid) #rename columns for matching 
ACE_Subsample
write.csv(ACE_Subsample,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_Subsample.csv", 
          row.names = FALSE)

HER42PB_Subsample <- ACE_Subsample %>% 
  filter(Site == "HER42PB") %>% 
  print()
write.csv(HER42PB_Subsample,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB_Subsample.csv", 
          row.names = FALSE)

BI10_Subsample <- ACE_Subsample %>% 
  filter(Site == "BI10") %>% 
  print()
write.csv(BI10_Subsample,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/BI10_Subsample.csv", 
          row.names = FALSE)

# Import qc and element filtered ITRAX cps dataset -------------------------------------------------------------
ACE_itrax <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_cps.csv")
ACE_itrax

# Select HER42PB site data - use find/replace for other sites
ACE_itrax_HER42PB_cps <- ACE_itrax %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "HER42PB") %>% 
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:Total_scatter) %>% 
  print()
write.csv(ACE_itrax_HER42PB_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_itrax_HER42PB_cps.csv", 
          row.names = FALSE)

# Select BI10 site data - use find/replace for other sites
ACE_itrax_BI10_cps <- ACE_itrax %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "BI10") %>% 
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:Total_scatter) %>% 
  print()
write.csv(ACE_itrax_BI10_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_itrax_BI10_cps.csv", 
          row.names = FALSE)


# Matching Subsample - XRF ----------------------------------------------------------------


# HER42PB ----------------------------------------------------------------

# Choose and rename datasets to use for matching
ACE_Subsample_HER42PB <- HER42PB_Subsample
ACE_Subsample_HER42PB

ACE_itrax_HER42PB <- ACE_itrax_HER42PB_cps %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% 
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc)
ACE_itrax_HER42PB

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
  relocate(sample, .after = Section) %>% 
  mutate(Ln_C_org = log(C_org)) %>%
  mutate(Ln_C_org_err = (C_org_err/C_org)*Ln_C_org) %>% 
  relocate(Ln_C_org, .after = C_org_err) %>% 
  relocate(Ln_C_org_err, .after = Ln_C_org)
HER42PB_Subsample_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure

# XRF elements defined by Francois & ITRAX acf
Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe",
                  "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd",
                     "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

HER42PB_Subsample_xrf_matched <- HER42PB_Subsample_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(Elements_min, Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
HER42PB_Subsample_xrf_matched
HER42PB_Subsample_xrf_matched$Zn # check no 0 values present

write.csv(HER42PB_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB_Subsample_xrf_matched.csv", 
          row.names = FALSE)

# BI10 ----------------------------------------------------------------

# Choose and rename datasets to use for matching
ACE_Subsample_BI10 <- BI10_Subsample
ACE_Subsample_BI10

ACE_itrax_BI10 <- ACE_itrax_BI10_cps %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% 
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc)
ACE_itrax_BI10

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
  mutate(Ln_C_org = log(C_org)) %>%
  mutate(Ln_C_org_err = (C_org_err/C_org)*Ln_C_org) %>% 
  relocate(Ln_C_org, .after = C_org_err) %>% 
  relocate(Ln_C_org_err, .after = Ln_C_org) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
BI10_Subsample_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure

# XRF elements defined by Francois & ITRAX acf
Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe",
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd",
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

BI10_Subsample_xrf_matched <- BI10_Subsample_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(Elements_min, Elements_min_sd))), 
          ~ (. == 0) * min(.[. != 0])/2 + .)
BI10_Subsample_xrf_matched
BI10_Subsample_xrf_matched$Zn # check no 0 values present

write.csv(BI10_Subsample_xrf_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/BI10_Subsample_xrf_matched.csv", 
          row.names = FALSE)





