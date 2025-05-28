# Figure S12: ICPMS - Subsample - XRF dataset matching 

# n = 266 defined by size of ICPMS dataset
# reduced from n = 272 due to new subsample composite dataset

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
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis', 
              'ggrepel', 'itraxR', 'PeriodicTable', 'errors', 'forecast', 'broom',
              'directlabels', 'performance', 'lmtest', 'ggpmisc', 'cowplot', 'Hmisc','car')
lapply(packages, library, character.only=TRUE)


# 1) ICPMS dataset & Subsample dataset matching  -----------------------------
# Import Subsample data  ---------------------------
subsample_param <- c("Water_Content_pc", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density_g_cm3","Dry_density_g_cm3",	"Dry_density_err_g_cm3", 
                     "DMAR_g_cm2_yr", "DMAR_err_g_cm2_yr"	)

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density_g_cm3")

ACE_Subsample <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_Subsample.csv") %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all inf to NA
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  select(c(Location:Section, Sample, strat_depth_top: accrate_err, Water_Content_pc: DMAR_err_pc)) %>% 
  select(-Ash_Mass) %>% 
  rename(sample = Sample, top = strat_depth_top, bottom = strat_depth_bottom, 
         mp = strat_depth_mid) #rename columns for matching 
ACE_Subsample
head(ACE_Subsample)
tail(ACE_Subsample)

# Import ICPMS composite data --------------------------------------------------
ACE_ICPMS <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_SHW_ICPMS_Composite.csv")
ACE_ICPMS
head(ACE_ICPMS)
tail(ACE_ICPMS)

# BI10 ICPMS & Subsample matching ---------------------------------------

# Select BI10 site data for both datasets
ACE_ICPMS_BI10 <- ACE_ICPMS %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "BI10") %>%
  arrange(Sample_ID) %>% 
  rename(sample = Sample_ID) %>% 
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
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section)
BI10_Subsample_ICPMS_matched
write.csv(BI10_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/BI10/BI10_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# HER42PB ICPMS & Subsample matching ---------------------------------------

# Select HER42PB site data for both datasets
ACE_ICPMS_HER42PB <- ACE_ICPMS %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "HER42PB") %>%
  arrange(Sample_ID) %>% 
  rename(sample = Sample_ID) %>% 
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
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section)
HER42PB_Subsample_ICPMS_matched
write.csv(HER42PB_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/HER42PB/HER42PB_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# KER1 ICPMS & Subsample matching ---------------------------------------

# Select KER1 site data for both datasets
ACE_ICPMS_KER1 <- ACE_ICPMS %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "KER1") %>%
  arrange(Sample_ID) %>% 
  rename(sample = Sample_ID) %>% 
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
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section)
KER1_Subsample_ICPMS_matched
write.csv(KER1_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER1/KER1_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# KER3 ICPMS & Subsample matching ---------------------------------------

# Select KER3 site data for both datasets
ACE_ICPMS_KER3 <- ACE_ICPMS %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "KER3") %>%
  arrange(Sample_ID) %>% 
  rename(sample = Sample_ID) %>% 
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
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section)
KER3_Subsample_ICPMS_matched
write.csv(KER3_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER3/KER3_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)


# PB1 ICPMS & Subsample matching ---------------------------------------

# Select PB1 site data for both datasets
ACE_ICPMS_PB1 <- ACE_ICPMS %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "PB1") %>%
  arrange(Sample_ID) %>% 
  rename(sample = Sample_ID) %>% 
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
  filter(!if_any(everything(), is.na)) %>% #remove rows will all NAs
  relocate(Location:DMAR_err_pc, .after = sample) %>% 
  relocate(sample, .after = Section)
PB1_Subsample_ICPMS_matched
write.csv(PB1_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/PB1/PB1_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)

# ACE ----------------------------

ACE_Subsample_ICPMS_matched<- bind_rows(BI10_Subsample_ICPMS_matched,
                                      HER42PB_Subsample_ICPMS_matched,
                                      KER1_Subsample_ICPMS_matched,
                                      KER3_Subsample_ICPMS_matched,
                                      PB1_Subsample_ICPMS_matched)
ACE_Subsample_ICPMS_matched

write.csv(ACE_Subsample_ICPMS_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/ACE/ACE_Subsample_ICPMS_matched.csv", 
          row.names = FALSE)






# 2) XRF dataset & ICPMS+Subsample dataset matching ----------------------------
# Import XRF cps dataset (qc and element filtered) -----------------------------
ACE_itrax <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_cps.csv")
ACE_itrax

# Select BI10 site data - use find/replace for other sites
ACE_itrax_BI10_cps <- ACE_itrax %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "BI10") %>% 
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:Total_scatter) %>% 
  print()
write.csv(ACE_itrax_BI10_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/BI10/ACE_itrax_BI10_cps.csv", 
          row.names = FALSE)

# Select HER42PB site data - use find/replace for other sites
ACE_itrax_HER42PB_cps <- ACE_itrax %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "HER42PB") %>% 
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:Total_scatter) %>% 
  print()
write.csv(ACE_itrax_HER42PB_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/HER42PB/ACE_itrax_HER42PB_cps.csv", 
          row.names = FALSE)

# Select KER1 site data - use find/replace for other sites
ACE_itrax_KER1_cps <- ACE_itrax %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "KER1") %>% 
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:Total_scatter) %>% 
  print()
write.csv(ACE_itrax_KER1_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER1/ACE_itrax_KER1_cps.csv", 
          row.names = FALSE)

# Select KER3 site data - use find/replace for other sites
ACE_itrax_KER3_cps <- ACE_itrax %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "KER3") %>% 
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:Total_scatter) %>% 
  print()
write.csv(ACE_itrax_KER3_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER3/ACE_itrax_KER3_cps.csv", 
          row.names = FALSE)

# Select PB1 site data - use find/replace for other sites
ACE_itrax_PB1_cps <- ACE_itrax %>% 
  filter(!if_any(everything(), is.na)) %>% # remove rows will all NAs
  filter(Site == "PB1") %>% 
  filter(qc == TRUE) %>%
  select(depth:surface, kcps:Total_scatter) %>% 
  print()
write.csv(ACE_itrax_PB1_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/PB1/ACE_itrax_PB1_cps.csv", 
          row.names = FALSE)

# BI10 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_BI10 <- BI10_Subsample_ICPMS_matched

# or Import from existing file
#ACE_Subsample_ICPMS_BI10 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/BI10/BI10_Subsample_ICPMS_matched.csv", row.names = FALSE)

ACE_itrax_BI10 <- ACE_itrax_BI10_cps %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% 
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc)
ACE_itrax_BI10

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
  mutate(Ln_C_org = log(C_org_pc)) %>%
  mutate(Ln_C_org_err = (C_org_pc_err/C_org_pc)*Ln_C_org) %>% 
  relocate(Ln_C_org, .after = C_org_pc_err) %>% 
  relocate(Ln_C_org_err, .after = Ln_C_org) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
BI10_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe",
                  "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd",
                     "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

BI10_Subsample_ICPMS_XRF_matched <- BI10_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(Elements_min, Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
BI10_Subsample_ICPMS_XRF_matched
BI10_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
BI10_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present

write.csv(BI10_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/BI10/BI10_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# HER42PB Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_HER42PB <- HER42PB_Subsample_ICPMS_matched

# or Import from existing file
#ACE_Subsample_ICPMS_HER42PB <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/HER42PB/HER42PB_Subsample_ICPMS_matched.csv", row.names = FALSE)

ACE_itrax_HER42PB <- ACE_itrax_HER42PB_cps %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% 
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc)
ACE_itrax_HER42PB

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
  relocate(sample, .after = Section) %>% 
  mutate(Ln_C_org = log(C_org_pc)) %>%
  mutate(Ln_C_org_err = (C_org_pc_err/C_org_pc)*Ln_C_org) %>% 
  relocate(Ln_C_org, .after = C_org_pc_err) %>% 
  relocate(Ln_C_org_err, .after = Ln_C_org) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
HER42PB_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe",
                  "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd",
                     "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

HER42PB_Subsample_ICPMS_XRF_matched <- HER42PB_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(Elements_min, Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
HER42PB_Subsample_ICPMS_XRF_matched
HER42PB_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
HER42PB_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(HER42PB_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/HER42PB/HER42PB_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# KER1 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_KER1 <- KER1_Subsample_ICPMS_matched

# or Import from existing file
#ACE_Subsample_ICPMS_KER1 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER1/KER1_Subsample_ICPMS_matched.csv", row.names = FALSE)

ACE_itrax_KER1 <- ACE_itrax_KER1_cps %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% 
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc)
ACE_itrax_KER1

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
  relocate(sample, .after = Section) %>% 
  mutate(Ln_C_org = log(C_org_pc)) %>%
  mutate(Ln_C_org_err = (C_org_pc_err/C_org_pc)*Ln_C_org) %>% 
  relocate(Ln_C_org, .after = C_org_pc_err) %>% 
  relocate(Ln_C_org_err, .after = Ln_C_org) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
KER1_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe",
                  "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd",
                     "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

KER1_Subsample_ICPMS_XRF_matched <- KER1_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(Elements_min, Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
KER1_Subsample_ICPMS_XRF_matched
KER1_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
KER1_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(KER1_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER1/KER1_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# KER3 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_KER3 <- KER3_Subsample_ICPMS_matched

# or Import from existing file
#ACE_Subsample_ICPMS_KER3 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER3/KER3_Subsample_ICPMS_matched.csv", row.names = FALSE)

ACE_itrax_KER3 <- ACE_itrax_KER3_cps %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% 
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc)
ACE_itrax_KER3

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
  relocate(sample, .after = Section) %>% 
  mutate(Ln_C_org = log(C_org_pc)) %>%
  mutate(Ln_C_org_err = (C_org_pc_err/C_org_pc)*Ln_C_org) %>% 
  relocate(Ln_C_org, .after = C_org_pc_err) %>% 
  relocate(Ln_C_org_err, .after = Ln_C_org) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
KER3_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe",
                  "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd",
                     "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

KER3_Subsample_ICPMS_XRF_matched <- KER3_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(Elements_min, Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
KER3_Subsample_ICPMS_XRF_matched
KER3_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
KER3_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(KER3_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/KER3/KER3_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# PB1 Subsample-ICPMS & XRF matching ---------------------------------------

# Rename Subsample-ICPMS combined dataset 
ACE_Subsample_ICPMS_PB1 <- PB1_Subsample_ICPMS_matched

# or Import from existing file
#ACE_Subsample_ICPMS_PB1 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/PB1/PB1_Subsample_ICPMS_matched.csv", row.names = FALSE)

ACE_itrax_PB1 <- ACE_itrax_PB1_cps %>% 
  mutate(Ln_inc_coh = log(inc_coh)) %>% 
  relocate(Ln_inc_coh, .after = inc_coh) %>%
  mutate(Ln_coh_inc = log(coh_inc)) %>% 
  relocate(Ln_coh_inc, .after = coh_inc)
ACE_itrax_PB1

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
  relocate(sample, .after = Section) %>% 
  mutate(Ln_C_org = log(C_org_pc)) %>%
  mutate(Ln_C_org_err = (C_org_pc_err/C_org_pc)*Ln_C_org) %>% 
  relocate(Ln_C_org, .after = C_org_pc_err) %>% 
  relocate(Ln_C_org_err, .after = Ln_C_org) %>% 
  filter(!(sample == '475')) # remove LOI data marked as 'stone'
PB1_Subsample_ICPMS_xrf_matched0

# Replace 0 values in each column with half column min to allow lm and log to function
# Recommended procedure in Bertrand et al. (submitted) to retain dataframe structure
# XRF elements defined by Francois & ITRAX acf
Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe",
                  "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd",
                     "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

PB1_Subsample_ICPMS_XRF_matched <- PB1_Subsample_ICPMS_xrf_matched0 %>% 
  mutate_at(vars(all_of(c(Elements_min, Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .)
PB1_Subsample_ICPMS_XRF_matched
PB1_Subsample_ICPMS_XRF_matched$Zn # check no 0 values present
PB1_Subsample_ICPMS_XRF_matched$coh_inc # check coh/inc values are present
write.csv(PB1_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/PB1/PB1_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)

# ACE ----------------------------

ACE_Subsample_ICPMS_XRF_matched<- bind_rows(BI10_Subsample_ICPMS_XRF_matched,
                                        HER42PB_Subsample_ICPMS_XRF_matched,
                                        KER1_Subsample_ICPMS_XRF_matched,
                                        KER3_Subsample_ICPMS_XRF_matched,
                                        PB1_Subsample_ICPMS_XRF_matched)
ACE_Subsample_ICPMS_XRF_matched

write.csv(ACE_Subsample_ICPMS_XRF_matched,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/ACE/ACE_Subsample_ICPMS_XRF_matched.csv", 
          row.names = FALSE)




# 3) Corrrelation - Dry Mass vs coh/inc-------------------------------------------------------------
# ACE OLS & WLS ------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
ACE_DM_lm <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched)
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_OLS_summary.txt")
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_wlm <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weight = 1/(Dry_mass_err)^2)
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_OLS_wt_summary.txt")
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
ACE_model <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched) # define model
ACE_wt <- 1 / lm(abs(ACE_model$residuals) ~ ACE_model$fitted.values)$fitted.values^2 #define weights to use
ACE_wls <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weights=ACE_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls) # summary stats
glance(ACE_wls) # summary stats including AIC
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_WLS_summary.txt")
summary(ACE_wls)
glance(ACE_wls)
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls) # Performance package summary check for heteroscedasticity
icc(ACE_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
ACE_model_wt <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weight = 1/Dry_mass_err^2) # define model
ACE_wt_wt <- 1 / lm(abs(ACE_model_wt$residuals) ~ ACE_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_wls_wt <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weights=ACE_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls_wt) # summary stats
glance(ACE_wls_wt) # summary stats including AIC
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_WLS_wt_summary.txt")
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
ACE_DM_lm_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_DM_lm_influential_names,] # outliers only using of index values
ACE_DM_lm_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_DM_lm_lev.pdf")
leveragePlots(ACE_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_DM_lm_lev_cooks.pdf")
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
ACE_wlm_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_wlm_influential_names,] # outliers only using of index values
ACE_wlm_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wlm_lev.pdf")
leveragePlots(ACE_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wlm_lev_cooks.pdf")
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
ACE_wls_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_wls_influential_names,] # outliers only using of index values
ACE_wls_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wls_lev_bar.pdf")
barplot(hatvalues(ACE_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wls_lev.pdf")
leveragePlots(ACE_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wls_lev_cooks.pdf")
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
ACE_wls_wt_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_wls_wt_influential_names,] # outliers only using of index values
ACE_wls_wt_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wls_wt_lev.pdf")
leveragePlots(ACE_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "ACE"
element_title <- "coh_inc"
theme_set(theme_classic(10))
ACE <- ggplot(ACE_Subsample_ICPMS_XRF_matched, aes(x = coh_inc, y = Dry_mass)) +
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
  #xlim(0.15, 0.4) +
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
  m <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched);
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
  m <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weight = 1/(Dry_mass_err)^2);
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/FigS12/OLS_WLS/ACE_267/FigS12_ACE_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")




