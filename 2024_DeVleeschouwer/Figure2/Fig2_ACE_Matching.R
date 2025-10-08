# Fig2_Subsample_ICPMS_XRF_matching

# 10/6/25: depth matching now defined by Subsample-ICPMS matched composite dataset
# n = 265 reduced from n = 272 compared to previous
# Code adapted from "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/ACE/ACE_matching_cps.R" 

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

# 1) Subsample & ICPMS dataset matching ----------------------------------------
# Import Subsample data  ---------------------------
subsample_param <- c("Water_Content", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density","Dry_density",	"Dry_density_err", 
                     "DMAR", "DMAR_err", "DMAR_err_pc")

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density")

ACE_Subsample <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/ACE_Subsample.csv") %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all inf to NA
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  select(c(Location:Section, Sample, strat_depth_top: accrate_err, Water_Content: DMAR_err_pc)) %>% 
  select(-Ash_Mass) %>%   
  rename(sample = Sample, top = strat_depth_top, bottom = strat_depth_bottom, 
         mp = strat_depth_mid) #rename columns for matching 
ACE_Subsample
head(ACE_Subsample)
tail(ACE_Subsample)

# Import ICPMS composite datafile --------------------------------------------------

# Define ICPMS elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

ACE_ICPMS <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/ACE_SHW_ICPMS_Composite.csv") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  arrange(Sample_ID) %>% 
  rename(sample = Sample_ID) %>% #rename Sample_ID column to match sample in ACE_Subsample
  print()
ACE_ICPMS
head(ACE_ICPMS)
tail(ACE_ICPMS)

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


# 2) XRF dataset & ICPMS+Subsample dataset matching ----------------------------
# Define lists of elements to use ----------------------------------------------

# XRF-CS acf elements matched to Francois ICPMS elements
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")
# Scatter parameters
scatter_param <- c("Mo_inc", "Mo_coh", "inc_coh", "Ln_inc_coh", "coh_inc", "Ln_coh_inc", "Total_scatter",
                   "Mo_inc_sd", "Mo_coh_sd", "inc_coh_sd", "Ln_inc_coh_sd", "coh_inc_sd", "Ln_coh_inc_sd", "Total_scatter_sd")

# Import ITRAX XRF-CS cps qc dataset -------------------------------------------
ACE_itrax <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/ACE_xrf_qc.csv") %>% 
# originally from "Papers_R/2024_DeVleeschouwer/ACE_matching/Output/itrax_Composite/Matching_mean/ACE_xrf_qc.csv"
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

# BI10 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_BI10 <- BI10_Subsample_ICPMS_matched

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
BI10_ss_icpms_xrf_matched <- itrax_reduce(ACE_itrax_BI10, names = ACE_Subsample_ICPMS_BI10$sample,
                                    breaks_lower = ACE_Subsample_ICPMS_BI10$top,
                                    breaks_upper = ACE_Subsample_ICPMS_BI10$bottom,
                                    fun = mean,
                                    edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_BI10, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
BI10_ss_icpms_xrf_matched

# Calculate stdev for ITRAX matching data
BI10_ss_icpms_xrf_matched_sd <- itrax_reduce(ACE_itrax_BI10, names = ACE_Subsample_ICPMS_BI10$sample,
                                       breaks_lower = ACE_Subsample_ICPMS_BI10$top,
                                       breaks_upper = ACE_Subsample_ICPMS_BI10$bottom,
                                       fun = sd,
                                       edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_BI10, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
BI10_ss_icpms_xrf_matched_sd

# Merge mean & sd data into the same dataframe
BI10_Subsample_ICPMS_xrf_matched0 <-  BI10_ss_icpms_xrf_matched_sd %>% 
  inner_join(., BI10_ss_icpms_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
BI10_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
BI10_Subsample_ICPMS_XRF_matched <- BI10_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
BI10_Subsample_ICPMS_XRF_matched
BI10_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
BI10_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present

write.csv(BI10_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/BI10/BI10_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# HER42PB Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_HER42PB <- HER42PB_Subsample_ICPMS_matched

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
HER42PB_ss_icpms_xrf_matched <- itrax_reduce(ACE_itrax_HER42PB, names = ACE_Subsample_ICPMS_HER42PB$sample,
                                          breaks_lower = ACE_Subsample_ICPMS_HER42PB$top,
                                          breaks_upper = ACE_Subsample_ICPMS_HER42PB$bottom,
                                          fun = mean,
                                          edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_HER42PB, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
HER42PB_ss_icpms_xrf_matched

# Calculate stdev for ITRAX matching data
HER42PB_ss_icpms_xrf_matched_sd <- itrax_reduce(ACE_itrax_HER42PB, names = ACE_Subsample_ICPMS_HER42PB$sample,
                                             breaks_lower = ACE_Subsample_ICPMS_HER42PB$top,
                                             breaks_upper = ACE_Subsample_ICPMS_HER42PB$bottom,
                                             fun = sd,
                                             edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_HER42PB, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
HER42PB_ss_icpms_xrf_matched_sd

# Merge mean & sd data into the same dataframe
HER42PB_Subsample_ICPMS_xrf_matched0 <-  HER42PB_ss_icpms_xrf_matched_sd %>% 
  inner_join(., HER42PB_ss_icpms_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
HER42PB_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
HER42PB_Subsample_ICPMS_XRF_matched <- HER42PB_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
HER42PB_Subsample_ICPMS_XRF_matched
HER42PB_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
HER42PB_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(HER42PB_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# KER1 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_KER1 <- KER1_Subsample_ICPMS_matched

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
KER1_ss_icpms_xrf_matched <- itrax_reduce(ACE_itrax_KER1, names = ACE_Subsample_ICPMS_KER1$sample,
                                             breaks_lower = ACE_Subsample_ICPMS_KER1$top,
                                             breaks_upper = ACE_Subsample_ICPMS_KER1$bottom,
                                             fun = mean,
                                             edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_KER1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
KER1_ss_icpms_xrf_matched

# Calculate stdev for ITRAX matching data
KER1_ss_icpms_xrf_matched_sd <- itrax_reduce(ACE_itrax_KER1, names = ACE_Subsample_ICPMS_KER1$sample,
                                                breaks_lower = ACE_Subsample_ICPMS_KER1$top,
                                                breaks_upper = ACE_Subsample_ICPMS_KER1$bottom,
                                                fun = sd,
                                                edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_KER1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
KER1_ss_icpms_xrf_matched_sd

# Merge mean & sd data into the same dataframe
KER1_Subsample_ICPMS_xrf_matched0 <-  KER1_ss_icpms_xrf_matched_sd %>% 
  inner_join(., KER1_ss_icpms_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
KER1_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
KER1_Subsample_ICPMS_XRF_matched <- KER1_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
KER1_Subsample_ICPMS_XRF_matched
KER1_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
KER1_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(KER1_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/KER1/KER1_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# KER3 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_KER3 <- KER3_Subsample_ICPMS_matched

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
KER3_ss_icpms_xrf_matched <- itrax_reduce(ACE_itrax_KER3, names = ACE_Subsample_ICPMS_KER3$sample,
                                          breaks_lower = ACE_Subsample_ICPMS_KER3$top,
                                          breaks_upper = ACE_Subsample_ICPMS_KER3$bottom,
                                          fun = mean,
                                          edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_KER3, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
KER3_ss_icpms_xrf_matched

# Calculate stdev for ITRAX matching data
KER3_ss_icpms_xrf_matched_sd <- itrax_reduce(ACE_itrax_KER3, names = ACE_Subsample_ICPMS_KER3$sample,
                                             breaks_lower = ACE_Subsample_ICPMS_KER3$top,
                                             breaks_upper = ACE_Subsample_ICPMS_KER3$bottom,
                                             fun = sd,
                                             edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_KER3, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
KER3_ss_icpms_xrf_matched_sd

# Merge mean & sd data into the same dataframe
KER3_Subsample_ICPMS_xrf_matched0 <-  KER3_ss_icpms_xrf_matched_sd %>% 
  inner_join(., KER3_ss_icpms_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
KER3_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf

KER3_Subsample_ICPMS_XRF_matched <- KER3_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
KER3_Subsample_ICPMS_XRF_matched
KER3_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
KER3_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(KER3_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/KER3/KER3_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# PB1 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_PB1 <- PB1_Subsample_ICPMS_matched

# Match ITRAX cps to ICP dataframe using itrax_reduce in itrax.R & mean function
PB1_ss_icpms_xrf_matched <- itrax_reduce(ACE_itrax_PB1, names = ACE_Subsample_ICPMS_PB1$sample,
                                          breaks_lower = ACE_Subsample_ICPMS_PB1$top,
                                          breaks_upper = ACE_Subsample_ICPMS_PB1$bottom,
                                          fun = mean,
                                          edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_PB1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>% 
  relocate(Location:accrate_err, .before = sample) %>% 
  relocate(sample, .after = Section)
PB1_ss_icpms_xrf_matched

# Calculate stdev for ITRAX matching data
PB1_ss_icpms_xrf_matched_sd <- itrax_reduce(ACE_itrax_PB1, names = ACE_Subsample_ICPMS_PB1$sample,
                                             breaks_lower = ACE_Subsample_ICPMS_PB1$top,
                                             breaks_upper = ACE_Subsample_ICPMS_PB1$bottom,
                                             fun = sd,
                                             edges = c(">=", "<=")) %>%
  select(resample_names, kcps:Total_scatter) %>%
  rename(sample = resample_names)  %>% 
  inner_join(., ACE_Subsample_ICPMS_PB1, by = "sample") %>%
  rename(min_depth = top, max_depth = bottom, midpoint = mp) %>%
  rename_with(.fn = function(.x){paste0(.x,"_sd")}) %>% 
  rename(sample = sample_sd) %>% 
  select(sample, kcps_sd:Total_scatter_sd)
PB1_ss_icpms_xrf_matched_sd

# Merge mean & sd data into the same dataframe
PB1_Subsample_ICPMS_xrf_matched0 <-  PB1_ss_icpms_xrf_matched_sd %>% 
  inner_join(., PB1_ss_icpms_xrf_matched, by = "sample") %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows with NAs
  relocate(kcps_sd:Total_scatter_sd, .after = Total_scatter) %>% 
  relocate(sample, .after = Section)
PB1_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
PB1_Subsample_ICPMS_XRF_matched <- PB1_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
PB1_Subsample_ICPMS_XRF_matched
PB1_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
PB1_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(PB1_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/PB1/PB1_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# Create ACE_Subsample_ICPMS_XRF_matched dataset n = 265 -----------------------

ACE_Subsample_ICPMS_XRF_matched<- bind_rows(BI10_Subsample_ICPMS_XRF_matched,
                                        HER42PB_Subsample_ICPMS_XRF_matched,
                                        KER1_Subsample_ICPMS_XRF_matched,
                                        KER3_Subsample_ICPMS_XRF_matched,
                                        PB1_Subsample_ICPMS_XRF_matched)
ACE_Subsample_ICPMS_XRF_matched

write_csv(ACE_Subsample_ICPMS_XRF_matched, "Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_Subsample_ICPMS_XRF_matched.csv")
write_csv(ACE_Subsample_ICPMS_XRF_matched, "Papers_R/2024_DeVleeschouwer/Figure3/Data/Input/ACE_Subsample_ICPMS_XRF_matched.csv")
write_csv(ACE_Subsample_ICPMS_XRF_matched, "Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/ACE_Subsample_ICPMS_XRF_matched.csv")
write.csv(ACE_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)
write.csv(ACE_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/ACE/ACE_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)





