# Figure 4

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
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis' )
lapply(packages, library, character.only=TRUE)

# Set plotting parameters & universal plot size for base R  --------------------
 
# cm graph plot size *from* cm *to* inches:
plotinch_x <- 2 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch_y <- 7 / cm(1) # -> 8 cm  is  3.149606 inches
# or use aspect ratio
aspect_ratio <- 2.5
options(repr.plot.width=plotinch_x, repr.plot.height=plotinch_y)

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
mscl_param  <- c("Den1_SAT", "MS1_SAT", "DCMS1_SAT", "Impedance_SAT", 
                 "Fract_Porosity_SAT", "Resistivity_SAT")

# Subsample parameters
subsample_param <- c("Water_Content", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density","Dry_density",	"Dry_density_err", 
                     "DMAR", "DMAR_err", "DMAR_err_pc")

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density")


# Import datasets, standardise & centre (Z-scores) -----------------------------
# XRF-CS -----------------------------------------------------------------------
# cps --------------------------------------------------------------------------

ACE_xrf_cps <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_cps.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1")
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc)
ACE_xrf_cps

# cps Z-scores
ACE_xrf_cps <- ACE_xrf_cps 
ACE_xrf_cps[, acf_icp_Elements_key] <- scale(ACE_xrf_cps[, acf_icp_Elements_key], center = TRUE, scale = TRUE)
ACE_xrf_cps
write.csv(ACE_xrf_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_cps_Z.csv", row.names = FALSE)

# cps - convert to long format for facet plotting
acf_icp_Elements_key3 <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh") #no Zn, Mn for plotting
ACE_xrf_cps_long <- ACE_xrf_cps %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_cps_long

ACE_xrf_cps_long <- ACE_xrf_cps %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_cps_long

# log_inc -----------------------------------------------------------------

ACE_xrf_log_inc <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_log_inc.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1")
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc) %>%
  mutate(across(acf_icp_Elements_key, ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_log_inc

HER42PB_xrf_log_inc <- ACE_xrf_log_inc %>% 
  filter(Site=="HER42PB")
write.csv(HER42PB_xrf_log_inc,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_xrf_log_inc.csv", row.names = FALSE)

# log_inc Z-scores
ACE_xrf_log_inc_Z <- ACE_xrf_log_inc 
ACE_xrf_log_inc_Z[, acf_icp_Elements_key] <- scale(ACE_xrf_log_inc[, acf_icp_Elements_key], center = TRUE, scale = TRUE)
ACE_xrf_log_inc_Z
write.csv(ACE_xrf_log_inc_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_log_inc_Z.csv", row.names = FALSE)

HER42PB_xrf_log_inc_Z <- ACE_xrf_log_inc_Z %>% 
  filter(Site=="HER42PB")
write.csv(HER42PB_xrf_log_inc_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_xrf_log_inc_Z.csv", row.names = FALSE)

# log_inc - convert to long format for facet plotting
acf_icp_Elements_key3 <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh") #no Zn, Mn for plotting
ACE_xrf_log_inc_long <- ACE_xrf_log_inc %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_log_inc_long

ACE_xrf_log_inc_Z_long <- ACE_xrf_log_inc_Z %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_log_inc_Z_long

HER42PB_xrf_log_inc_long <- HER42PB_xrf_log_inc %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_xrf_log_inc_long

HER42PB_xrf_log_inc_Z_long <- HER42PB_xrf_log_inc_Z %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_xrf_log_inc_Z_long

# Plot vs Depth - key XRF_CS log_inc data
HER42PB_XRFCS_log_inc_depth <- ggplot(HER42PB_xrf_log_inc_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 0.2) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "HER42PB XRF-CS log_inc data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_XRFCS_log_inc_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/log_inc/Fig4.1_XRFCS_log_inc_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# clr ---------------------------------------------------------------------

ACE_xrf_clr <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_clr.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1")
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc)
ACE_xrf_clr

HER42PB_xrf_clr <- ACE_xrf_clr %>% 
  filter(Site=="HER42PB")
HER42PB_xrf_clr
write.csv(HER42PB_xrf_clr,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_xrf_clr.csv", row.names = FALSE)

# clr Z-scores 
ACE_xrf_clr_Z <- ACE_xrf_clr 
ACE_xrf_clr_Z[, acf_icp_Elements_key] <- scale(ACE_xrf_clr[, acf_icp_Elements_key], center = TRUE, scale = TRUE)
ACE_xrf_clr_Z
write.csv(ACE_xrf_clr_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_xrf_clr_Z.csv", row.names = FALSE)

HER42PB_xrf_clr_Z <- ACE_xrf_clr_Z %>% 
  filter(Site=="HER42PB")
write.csv(HER42PB_xrf_clr_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_xrf_clr_Z.csv", row.names = FALSE)

# clr - convert to long format for facet plotting
acf_icp_Elements_key3 <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr") #no "Mo_coh" for plotting - not included in clr elements
ACE_xrf_clr_long <- ACE_xrf_clr %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_clr_long

ACE_xrf_clr_Z_long <- ACE_xrf_clr_Z %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_clr_Z_long

HER42PB_xrf_clr_long <- HER42PB_xrf_clr %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_xrf_clr_long

HER42PB_xrf_clr_Z_long <- HER42PB_xrf_clr_Z %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_xrf_clr_Z_long

# Plot vs Depth - key XRF_CS clr data
HER42PB_XRFCS_clr_depth <- ggplot(HER42PB_xrf_clr_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 0.2) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "HER42PB XRF-CS clr data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_XRFCS_clr_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/clr/Fig4.1_XRFCS_clr_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")


# Matched XRF & ICPMS log data ----------------------

ACE_matched_xrf_icp_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_xrf_icp_matched_log_inc.csv") 
is.na(ACE_matched_xrf_icp_log_inc)<-sapply(ACE_matched_xrf_icp_log_inc, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_log_inc

# Standardise and centre dataframe - Z-scores
ACE_matched_xrf_icp_log_inc_Z <- ACE_matched_xrf_icp_log_inc 
ACE_matched_xrf_icp_log_inc_Z[, xrf_icp_elements] <- scale(ACE_matched_xrf_icp_log_inc[, xrf_icp_elements], center = TRUE, scale = TRUE)
ACE_matched_xrf_icp_log_inc_Z
write.csv(ACE_matched_xrf_icp_log_inc_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_subsample_xrf_icp_matched_log_inc_Z.csv", row.names = FALSE)

# matched - convert to long format for facet plotting
xrf_icp_elements3 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP",
                       "Fe", "Fe_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP", "Mo_coh")
xrf_icp_elements3_xrf <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")
xrf_icp_elements3_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP","Fe_ICP", "Sr_ICP", "Zr_ICP", "Mo_coh")

ACE_matched_xrf_icp_log_inc_long <- ACE_matched_xrf_icp_log_inc %>% 
  select(c(all_of(xrf_icp_elements3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`xrf_icp_elements3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_long

ACE_matched_xrf_icp_log_inc_Z_long <- ACE_matched_xrf_icp_log_inc_Z %>% 
  select(c(all_of(xrf_icp_elements3), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`xrf_icp_elements3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_Z_long

ACE_matched_xrf_icp_log_inc_long_xrf <- ACE_matched_xrf_icp_log_inc %>% 
  select(c(all_of(xrf_icp_elements3_xrf), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_xrf`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_long_xrf

ACE_matched_xrf_icp_log_inc_Z_long_xrf <- ACE_matched_xrf_icp_log_inc_Z %>% 
  select(c(all_of(xrf_icp_elements3_xrf), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_xrf`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_Z_long_xrf

ACE_matched_xrf_icp_log_inc_long_icp <- ACE_matched_xrf_icp_log_inc %>% 
  select(c(all_of(xrf_icp_elements3_icp), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_icp`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_long_icp

ACE_matched_xrf_icp_log_inc_Z_long_icp <- ACE_matched_xrf_icp_log_inc_Z %>% 
  select(c(all_of(xrf_icp_elements3_icp), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_icp`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_Z_long_icp

# Summary depth plots 
library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo
# Plot vs Depth - key XRF_CS log_inc data
HER42PB_XRFCS_log_inc_depth <- ggplot(ACE_matched_xrf_icp_log_inc_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 0.2) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "HER42PB XRF-CS log_inc data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_XRFCS_log_inc_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/clr/Fig4.1_XRFCS_log_inc_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")




# MSCL data ---------------------------------------------------------------

ACE_MSCL <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_SHW_MSCL_Composite.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1")
ACE_MSCL
write.csv(ACE_MSCL,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_MSCL.csv", row.names = FALSE)

HER42PB_MSCL <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_SHW_MSCL_Composite.csv") %>% 
  filter(Site=="HER42PB")
HER42PB_MSCL
write.csv(HER42PB_MSCL,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_MSCL.csv", row.names = FALSE)

# Standardise and centre dataframe - Z-scores
ACE_MSCL_Z <- ACE_MSCL
ACE_MSCL_Z[, mscl_param] <- scale(ACE_MSCL[, mscl_param], center = TRUE, scale = TRUE)
ACE_MSCL_Z
write.csv(ACE_MSCL_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_MSCL_Z.csv", row.names = FALSE)

HER42PB_MSCL_Z <- HER42PB_MSCL
HER42PB_MSCL_Z[, mscl_param] <- scale(HER42PB_MSCL[, mscl_param], center = TRUE, scale = TRUE)
HER42PB_MSCL_Z
write.csv(HER42PB_MSCL_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_MSCL_Z.csv", row.names = FALSE)

# Convert to long format
ACE_MSCL_long <- ACE_MSCL %>% 
  select(all_of(mscl_param), Site, depth, SH20_mean_age) %>%
  pivot_longer(all_of(`mscl_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_MSCL_long

ACE_MSCL_long_Z <- ACE_MSCL %>% 
  select(all_of(mscl_param), Site, depth, SH20_mean_age) %>%
  pivot_longer(all_of(`mscl_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_MSCL_long

HER42PB_MSCL_long <- HER42PB_MSCL %>% 
  select(all_of(mscl_param), Site, depth, SH20_mean_age) %>%
  pivot_longer(c(`mscl_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_MSCL_long

HER42PB_MSCL_long_Z <- HER42PB_MSCL_Z %>% 
  select(all_of(mscl_param), Site, depth, SH20_mean_age) %>%
  pivot_longer(c(`mscl_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_MSCL_long_Z

library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Plot vs Depth - key MSCL data
HER42PB_MSCL_depth <- ggplot(HER42PB_MSCL_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "white", color = "black", size = 1) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "HER42PB MSCL data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_MSCL_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/MSCL/Fig4.1_MSCL_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")


# Subsample data ----------------------------------------------------------

ACE_Subsample <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_Subsample.csv")
ACE_Subsample
 
HER42PB_Subsample <- ACE_Subsample%>% 
  filter(Site == "HER42PB")
HER42PB_Subsample

# Standardise and centre dataframe - Z-scores
HER42PB_Subsample_Z <- HER42PB_Subsample
HER42PB_Subsample_Z[, subsample_param1] <- scale(HER42PB_Subsample[, subsample_param1], center = TRUE, scale = TRUE)
HER42PB_Subsample_Z
write.csv(HER42PB_Subsample_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB/HER42PB_Subsample_comp_Z.csv", row.names = FALSE)

HER42PB_Subsample_long <- HER42PB_Subsample %>% 
  select(c(all_of(subsample_param1), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param1`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_Subsample_long

HER42PB_Subsample_long_Z <- HER42PB_Subsample_Z %>% 
  select(c(all_of(subsample_param1), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param1`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
HER42PB_Subsample_long_Z

library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo
# Plot vs Depth - all subsample data
p1_HER42PB_Subsample_depth <- ggplot(HER42PB_Subsample_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 1) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "HER42PB Subsample data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p1_HER42PB_Subsample_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Subsample/Fig4.1b_HER42PB_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")





