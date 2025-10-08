# ICPMS and XRF error corrections for as ppm plots

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



# ICPMS and XRF error corrections for as ppm plots ----------------------------------------------------------

# Import ICPMS data - change errors to Table 2 (previously +/-5%)
ACE_ICP_ppm <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_matched_xrf_icp_cps.csv") %>%
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
  select(Site:Zr_ICP_sd) %>% 
  
  # Add 10% max error to XRF-CS cps data
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
write.csv(ACE_ICP_ppm,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/Fig4/ACE/ACE_matched_xrf_icp_cps.csv", row.names = FALSE)
