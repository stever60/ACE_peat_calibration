# Figure 4 - conversion to ppm and indiviual plots for comparison

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

# ------------------------------------------------------------------------------
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

# Key XRF and ICPMS elements
xrf_icp_elements <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                      "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                      "Zr", "Zr_ICP", "Mo_coh")

xrf_icp_elements1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                       "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP", "Mo_inc",  "Mo_coh", "inc_coh", "coh_inc", "Dry_mass") # Mo_inc included

xrf_icp_elements2 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                       "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP")

icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", 
                      "Ni_ICP", "Cu_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", 
                         "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# XRF-CS precitors & ICPMS response variables for PLS calibration model
main_elements_xrf <- c("K", "Ca", "Ti", "Mn","Fe", "Zn", "Rb", "Sr","Zr")
key_elements_xrf <- c("Ti", "Ca", "Mn", "Fe", "Sr", "Zr")
main_elements_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP","Fe_ICP", "Zn_ICP", "Rb_ICP", "Sr_ICP","Zr_ICP")
key_elements_icp <- c("Ti_ICP", "Ca_ICP", "Mn_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP")
final_elements_xrf <- c("Ti", "Ca", "Sr", "Zr")
final_elements_icp <- c("Ti_ICP", "Ca_ICP", "Sr_ICP", "Zr_ICP")

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


# ------------------------------------------------------------------------------
# Import datasets ----------------------------------------------------


ACE_xrf_log_inc <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_log_inc.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>% 
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc) %>%
  filter(qc == "TRUE") %>% 
  filter(validity == "1") %>% 
  mutate(across(acf_icp_Elements_key, ~ ifelse(. <=-10, NA, .))) # remove outliers > -10 log ICPMS value
ACE_xrf_log_inc

# WLS calibration: convert XRF-CS data to ppm log_inc +/- RMSE error -------------
ACE_WLS_Ti <- ACE_xrf_log_inc %>%
  select(Site, depth, SH20_age, Ti) %>%
  #filter (Site == "HER42PB") %>% 
  mutate(Ti_convert_WLS = 11+0.68*Ti) %>% # Ti Ln equation: y = 11+0.68 x  where y = Ln(Ti_ICP) and x = Ln(Ti_XRF)
  mutate(Ti_ppm_WLS = exp(Ti_convert_WLS)) %>% 
  mutate(Ti_upper_WLS = Ti+0.645) %>% # RMSE = 0.663
  mutate(Ti_convert_upper_WLS = 11+0.68*Ti_upper_WLS) %>%
  mutate(Ti_upper_RMSE_WLS = exp(Ti_convert_upper_WLS)) %>% 
  mutate(Ti_lower_WLS = Ti-0.645) %>% 
  mutate(Ti_convert_lower_WLS = 11+0.68*Ti_lower_WLS) %>%
  mutate(Ti_lower_RMSE_WLS = exp(Ti_convert_lower_WLS)) %>% 
  select(Site, depth, SH20_age, Ti, Ti_convert_WLS, Ti_ppm_WLS, Ti_lower_RMSE_WLS, Ti_upper_RMSE_WLS)
ACE_WLS_Ti
write.csv(ACE_WLS_Ti,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_WLS_Ti.csv", row.names = FALSE)

# PLS calibration: convert XRF-CS data to ppm log_inc +/- RMSE error -------------

ACE_PLS_ppm <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_PLS/Output/PLSR/log_inc/ACE_PLS_ppm_pred_log_inc.csv")

# Create WLS & PLS dataset for plotting & save to file
ACE_WLS_PLS <- ACE_WLS_Ti %>%
  select(-c(Site, depth, SH20_age, Ti)) %>% 
  bind_cols(ACE_PLS_ppm) %>% 
  relocate(c(Ti_convert_WLS:Ti_upper_RMSE_WLS), .after = Zr_up_RMSE)
write.csv(ACE_WLS_PLS,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_WLS_PLS.csv", row.names = FALSE)

# Import ICPMS Ti data for comparison / overlay
ACE_ICP <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_subsample_icp_xrf_matched_cps.csv")

# Import GEOTEK MSCL ACE dataset
ACE_MSCL <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_SHW_MSCL_Composite.csv") %>% 
  filter(Site == "BI10"| Site == "HER42PB" | Site == "KER1" | Site == "KER3" | Site == "PB1") %>%
  select(c(Location:SH20_median_age, mscl_param))
write.csv(ACE_MSCL,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE/ACE_MSCL.csv", row.names = FALSE)

# Import Subsample ACE dataset
ACE_Subsample <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_Subsample.csv")

# ------------------------------------------------------------------------------
# Plots ------------------------------------------------------------------------
# XRF-CS WLS - Ti plot -----------------------------------------------------------
HER42PB_WLS_PLS <- ACE_WLS_PLS %>% 
  filter(Site == "HER42PB")  
theme_set(theme_classic(12))
Site_title = "HER42PB"
Element_title = "Ti"

# Depth plot - XRF-CS
HER42PB_Ti <- 
  ggplot() +
  geom_lineh(data = HER42PB_WLS_PLS, aes(x=Ti_lower_RMSE_WLS, y=depth), color = "lightgrey") +
  geom_lineh(data = HER42PB_WLS_PLS, aes(x=Ti_upper_RMSE_WLS, y=depth), color = "lightgrey") +
  geom_lineh(data = HER42PB_WLS_PLS, aes(x=Ti_ppm_WLS, y=depth), color = "blue") +
  geom_lineh(data = ACE_ICP_ppm_Ti, aes(x=Ti_ICP, y=depth), color = "red") +
  #geom_errorbarh(data = ACE_ICP_ppm_Ti, aes(x=Ti_ICP, y=depth, xmin=Ti_ICP-Ti_ICP_sd, xmax=Ti_ICP+Ti_ICP_sd), 
  #               color = "red", height=0) +
  geom_point(data = ACE_ICP_ppm_Ti, aes(x=Ti_ICP, y=depth), fill = "white", color = "red", shape = 21, size = 2) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(0,7000) +
  labs(
    title = paste0(Site_title, ": WLS calibration for ", Element_title, " (log_inc model)"),
    x = bquote(.(Element_title) ~ "(mg kg"^{-1}*")"),
    y = "Depth (cm)"
  ) +
  #ggtitle(paste(Site_title, ": ", Element_title)) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
HER42PB_Ti
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.1_ICPMS_WLS_Ti.pdf",
        height = c(24), width = c(8), dpi = 600, units = "cm")



# XRF-CS PLS plot function with ICPMS overlay ----------------------------
plot_element_depth_log_inc_with_ICP <- function(site_name, element, pred_data, 
                                                icp_data, output_dir, age_model, 
                                                age_breaks = NULL) {
  
  # make sure tidypaleo is loaded for age secondary axis
  library(tidypaleo)
  
  # Define age breaks automatically if not provided
  if (is.null(age_breaks)) {
    age_range <- range(age_model$age, na.rm = TRUE)
    age_breaks <- pretty(age_range, n = 6)
  }
  
  # Dynamically construct column names
  element_pred <- sym(paste0(element, "_ICP_pred_ppm"))
  element_low  <- sym(paste0(element, "_low_RMSE"))
  element_up   <- sym(paste0(element, "_up_RMSE"))
  element_icp  <- sym(paste0(element, "_ICP"))
  
  # Build plot
  p <- ggplot(pred_data, aes(y = depth)) +
    geom_ribbon(aes(xmin = !!element_low, xmax = !!element_up), color = "grey80", fill = "grey80", na.rm = TRUE) + 
    #geom_point(aes(x = !!element_pred), color = "blue", size = 0.3) +
    geom_lineh(aes(x = !!element_pred), linewidth = 0.5, color = "blue", na.rm = TRUE) + 
    geom_lineh(data = icp_data, aes(x = !!element_icp, y = depth), color = "red", linewidth = 0.75, na.rm = TRUE) +
    geom_errorbarh(data = icp_data,aes_string(y = "depth", 
                                              xmin = paste0(element, "_ICP - ", element, "_ICP_sd"),
                                              xmax = paste0(element, "_ICP + ", element, "_ICP_sd")),
                   height = 0, color = "red", linewidth = 0.75) +
    geom_point(data = icp_data, aes(x = !!element_icp, y = depth), shape = 21, size = 2, color = "red", 
               fill = "white", stroke = 0.75, na.rm = TRUE) +  
    #scale_y_reverse() +
    scale_y_depth_age(
      age_model,
      name = "Depth (cm)",
      age_name = "Age (cal a BP)",
      age_breaks = age_breaks) +
    labs(
      title = paste0(site_name, ": PLSR ", element, " Predictions with Confidence Interval (log_inc model)"),
      x = bquote(.(element) ~ "(mg kg"^{-1}*")"),
      y = "Depth (cm)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(color = "black", size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.ticks.length = unit(0.3, "cm"),
      axis.line = element_line(size = 0.75),
      axis.ticks = element_line(size = 0.75)
    )
  
  # Save plot
  if (!is.null(output_dir)) {
    file_path <- file.path(output_dir, paste0(site_name, "_", element, "_log_inc_PLS_depth.pdf"))
    ggsave(file_path, plot = p, height = 24, width = 12, dpi = 600, units = "cm")
  }
  return(p)
}

# XRF-CS HER42PB PLS plot  ------------------------------------------------

# Inputs for depth function - save to file

HER42PB_ICP <- ACE_ICP %>% 
  filter(Site == "HER42PB")

# Define age-depth model in tidypaleo
HER42PB_adm <- age_depth_model(
  HER42PB_WLS_PLS, 
  depth = depth,
  age = SH20_age
)
# Plot to screen and save to names file directory
plot_element_depth_log_inc_with_ICP(
  site_name = "HER42PB",
  element = "Ti",
  pred_data = HER42PB_WLS_PLS,
  icp_data = HER42PB_ICP,
  output_dir = "Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/", 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
# Create site_depth for plot matrix
HER42PB_depth_Ti <- plot_element_depth_log_inc_with_ICP(
  site_name = "HER42PB",
  element = "Ti",
  pred_data = HER42PB_WLS_PLS,
  icp_data = HER42PB_ICP,
  output_dir = NULL, 
  age_model = HER42PB_adm,
  age_breaks = seq(0, 12000, by = 500)
)
HER42PB_depth_Ti



# MSCL wet density plot ---------------------------------------------------
HER42PB_Den_MSCL <- ACE_MSCL  %>% 
  filter(Site == "HER42PB") %>% 
  ggplot() +
  ggplot() +
  geom_lineh(aes(x=Den1_SAT, y=depth), color = "grey") +
  geom_point(aes(x=Den1_SAT, y=depth), fill = "grey", color = "grey", shape = 21, size = 1) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(1.2,1.5) +
  labs(x = paste("Wet GRD (g cm-3)") , y = paste0("Depth (cm)")) +
  ggtitle(paste(Site_title, ": MSCL")) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
HER42PB_Den_MSCL
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.2_Den_MSCL.pdf",
       height = c(24), width = c(8), dpi = 600, units = "cm")

# Subsample Dry Density plot --------------------------------------------
HER42PB_Den_Subsample <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_Subsample.csv") %>% 
  filter(Site == "HER42PB") %>% 
  ggplot() +
  geom_lineh(aes(x=Dry_density, y=depth), color = "black") +
  geom_point(aes(x=Dry_density, y=depth), fill = "white", color = "black", shape = 21, size = 2) +
  #geom_errorbarh(data = HER42PB_Subsample, aes(x=Dry_density, y=depth, xmin=Dry_density-Dry_density_err, xmax=Dry_density+Dry_density_err), color = "red", height=0) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(0,1) +
  labs(x = paste("Dry density (g cm-3)") , y = paste0("Depth (cm)")) +
  ggtitle(paste(Site_title, ": Subsample")) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
HER42PB_Den_Subsample
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.4_Den_subsample.pdf",
       height = c(24), width = c(8), dpi = 600, units = "cm")


# MSCL Density corrected MS (DCMS) plot --------------------------------------
HER42PB_DCMS1_SAT_MSCL <- ACE_MSCL %>% 
  filter(Site == "HER42PB")
  ggplot() +
  geom_lineh(aes(x=DCMS1_SAT, y=depth), color = "black") +
  geom_point(aes(x=DCMS1_SAT, y=depth), fill = "black", color = "black", shape = 21, size = 1) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(-2,10) +
  labs(x = paste("DC-MS (MSχ [κ/ρ]) x10-8 (m3 kg-1)") , y = paste0("Depth (cm)")) +
  ggtitle(paste(Site_title, ": MSCL")) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
HER42PB_DCMS1_SAT_MSCL
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.3_DCMS_MSCL.pdf",
       height = c(24), width = c(8), dpi = 600, units = "cm")


# Subsample Corg plot -------------------------------------------------------

HER42PB_C_Subsample <-
  ggplot() +
  geom_lineh(data = HER42PB_Subsample, aes(x=C_org, y=depth), color = "black") +
  geom_point(data = HER42PB_Subsample, aes(x=C_org, y=depth), fill = 'black', color = "black", shape = 21, size = 1) +
  #geom_errorbarh(data = HER42PB_Subsample, aes(x=C_org, y=depth, xmin=C_org-C_org_err, xmax=C_org+C_org_err), color = "red", height=0) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(0,50) +
  labs(x = paste("% Carbon") , y = paste0("Depth (cm)")) +
  ggtitle(paste(Site_title, ": Subsample")) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
HER42PB_C_Subsample
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.5_Corg.pdf",
       height = c(24), width = c(8), dpi = 600, units = "cm")

# XRF-CS inc/coh plot ----------------------------------------------------------
HER42PB_inc_coh <- ACE_xrf_log_inc %>% 
  filter(Site == "HER42PB") %>% 
  ggplot() +
  geom_lineh(data = HER42PB_xrf, aes(x=inc_coh, y=depth), color = "grey") +
  #geom_point(data = HER42PB_xrf, aes(x=inc_coh, y=depth), fill = "grey", color = "grey", shape = 21, size = 1) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(4,8) +
  labs(x = paste("XRF inc/coh") , y = paste0("Depth (cm)")) +
  ggtitle(paste(Site_title, ": inc/coh")) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
HER42PB_inc_coh
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.6_inc_coh_XRF.pdf",
       height = c(24), width = c(8), dpi = 600, units = "cm")


# Subsample dry mass plot -----------------------------------------------
HER42PB_DM_Subsample <-
  ggplot() +
  geom_lineh(data = HER42PB_Subsample, aes(x=Dry_mass, y=depth), color = "black") +
  geom_point(data = HER42PB_Subsample, aes(x=Dry_mass, y=depth), fill = 'black', color = "black", shape = 21, size = 1) +
  #geom_errorbarh(data = HER42PB_Subsample, aes(x=Dry_mass, y=depth, xmin=Dry_mass-Dry_mass_err, xmax=Dry_mass+Dry_mass_err), color = "red", height=0) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(0,50) +
  labs(x = paste("% Carbon") , y = paste0("Depth (cm)")) +
  ggtitle(paste(Site_title, ": Subsample")) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
HER42PB_C_Subsample
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.5_Corg.pdf",
       height = c(24), width = c(8), dpi = 600, units = "cm")


# XRF-CS coh/inc plot --------------------------------------------------
HER42PB_xrf <- ACE_xrf_log_inc %>% 
  filter(Site == "HER42PB")
p6_XRF_inc_coh <- 
  ggplot() +
  geom_lineh(data = HER42PB_xrf, aes(x=coh_inc, y=depth), color = "grey") +
  #geom_point(data = HER42PB_xrf, aes(x=inc_coh, y=depth), fill = "grey", color = "grey", shape = 21, size = 1) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(4,8) +
  labs(x = paste("XRF coh/inc") , y = paste0("Depth (cm)")) +
  ggtitle(paste(Site_title, ": inc/coh")) +
  theme(axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12, colour = "black"),
        title = element_text(size=12, colour = "black")) 
p6_XRF_inc_coh
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.6_inc_coh_XRF.pdf",
       height = c(24), width = c(8), dpi = 600, units = "cm")



# # COMBINE PLOTS ---------------------------------------------------------


ggarrange(HER42PB_C_Subsample, HER42PB_inc_coh,
          ncol = 4, nrow = 1, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4_carbon.pdf",
       height = c(30), width = c(40), dpi = 600, units = "cm")

ggarrange(HER42PB_Den_MSCL, HER42PB_Den_Subsample, HER42PB_Ti, HER42PB_DCMS1_SAT_MSCL,
          ncol = 4, nrow = 1, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Sites/HER42PB/Fig4.pdf",
       height = c(15), width = c(40), dpi = 600, units = "cm")

# ------------------------------------------------------------------------------
# ********** END ************
# ------------------------------------------------------------------------------
# Other data Plots to copy 
# BI10 subsample -------------------------------------------------------

BI10_Subsample <- ACE_Subsample%>% 
  filter(Site == "BI10")
BI10_Subsample

subsample_param <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density_g_cm3")

# Standardise and centre dataframe - Z-scores
BI10_Subsample_Z <- BI10_Subsample
BI10_Subsample_Z[, subsample_param] <- scale(BI10_Subsample[, subsample_param], center = TRUE, scale = TRUE)
BI10_Subsample_Z
write.csv(BI10_Subsample_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/BI10/BI10/BI10_Subsample_comp_Z.csv", row.names = FALSE)

BI10_Subsample_long <- BI10_Subsample %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
BI10_Subsample_long

BI10_Subsample_long_Z <- BI10_Subsample_Z %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
BI10_Subsample_long_Z

library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Plot vs Depth - all subsample data
p1_BI10_Subsample_depth <- ggplot(BI10_Subsample_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 1) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "BI10 Subsample data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p1_BI10_Subsample_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Subsample/Fig4.1b_BI10_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# KER1 subsample -------------------------------------------------------
  
  KER1_Subsample <- ACE_Subsample%>% 
  filter(Site == "KER1")
KER1_Subsample

subsample_param <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density_g_cm3")

# Standardise and centre dataframe - Z-scores
KER1_Subsample_Z <- KER1_Subsample
KER1_Subsample_Z[, subsample_param] <- scale(KER1_Subsample[, subsample_param], center = TRUE, scale = TRUE)
KER1_Subsample_Z
write.csv(KER1_Subsample_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/KER1/KER1_Subsample_comp_Z.csv", row.names = FALSE)

KER1_Subsample_long <- KER1_Subsample %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
KER1_Subsample_long

KER1_Subsample_long_Z <- KER1_Subsample_Z %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
KER1_Subsample_long_Z

library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Plot vs Depth - all subsample data
p1_KER1_Subsample_depth <- ggplot(KER1_Subsample_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 1) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "KER1 Subsample data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p1_KER1_Subsample_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Subsample/Fig4.1b_KER1_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# KER1 subsample -------------------------------------------------------

KER1_Subsample <- ACE_Subsample%>% 
  filter(Site == "KER1")
KER1_Subsample

subsample_param <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density_g_cm3")

# Standardise and centre dataframe - Z-scores
KER1_Subsample_Z <- KER1_Subsample
KER1_Subsample_Z[, subsample_param] <- scale(KER1_Subsample[, subsample_param], center = TRUE, scale = TRUE)
KER1_Subsample_Z
write.csv(KER1_Subsample_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/KER1/KER1_Subsample_comp_Z.csv", row.names = FALSE)

KER1_Subsample_long <- KER1_Subsample %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
KER1_Subsample_long

KER1_Subsample_long_Z <- KER1_Subsample_Z %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
KER1_Subsample_long_Z

library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Plot vs Depth - all subsample data
p1_KER1_Subsample_depth <- ggplot(KER1_Subsample_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 1) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "KER1 Subsample data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p1_KER1_Subsample_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Subsample/Fig4.1b_KER1_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# KER3 subsample -------------------------------------------------------

KER3_Subsample <- ACE_Subsample%>% 
  filter(Site == "KER3")
KER3_Subsample

subsample_param <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density_g_cm3")

# Standardise and centre dataframe - Z-scores
KER3_Subsample_Z <- KER3_Subsample
KER3_Subsample_Z[, subsample_param] <- scale(KER3_Subsample[, subsample_param], center = TRUE, scale = TRUE)
KER3_Subsample_Z
write.csv(KER3_Subsample_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/KER3/KER3_Subsample_comp_Z.csv", row.names = FALSE)

KER3_Subsample_long <- KER3_Subsample %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
KER3_Subsample_long

KER3_Subsample_long_Z <- KER3_Subsample_Z %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
KER3_Subsample_long_Z

library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Plot vs Depth - all subsample data
p1_KER3_Subsample_depth <- ggplot(KER3_Subsample_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 1) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "KER3 Subsample data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p1_KER3_Subsample_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Subsample/Fig4.1b_KER3_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# PB1 subsample -------------------------------------------------------

PB1_Subsample <- ACE_Subsample%>% 
  filter(Site == "PB1")
PB1_Subsample

subsample_param <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density_g_cm3")

# Standardise and centre dataframe - Z-scores
PB1_Subsample_Z <- PB1_Subsample
PB1_Subsample_Z[, subsample_param] <- scale(PB1_Subsample[, subsample_param], center = TRUE, scale = TRUE)
PB1_Subsample_Z
write.csv(PB1_Subsample_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/PB1/PB1_Subsample_comp_Z.csv", row.names = FALSE)

PB1_Subsample_long <- PB1_Subsample %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
PB1_Subsample_long

PB1_Subsample_long_Z <- PB1_Subsample_Z %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_mean_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
PB1_Subsample_long_Z

library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Plot vs Depth - all subsample data
p1_PB1_Subsample_depth <- ggplot(PB1_Subsample_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "black", color = "black", size = 1) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "", y = "Depth (cm)") +
  labs(title = "PB1 Subsample data") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p1_PB1_Subsample_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Subsample/Fig4.1b_PB1_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")




# Matched summary plots - not much use but code might be useful later ---------------------------------------

# xrf matched
# Select site & order that plots appear in
Matched_reorder_log_inc <- ACE_matched_xrf_icp_log_inc_long_xrf %>% # define order plots appear in
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")) 
Matched_reorder_log_inc_Z <- ACE_matched_xrf_icp_log_inc_Z_long_xrf %>% 
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")) 

# Depth - matched log_inc
HER42PB_Matched_log_inc_depth <- ggplot(Matched_reorder_log_inc, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "white", color = "black", size = 2) +
  #scale_x_reverse() +
  scale_y_reverse() +
  ylim(410, 0) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Depth (cm)", y = "Ln(E/inc.)") +
  labs(title = "HER42PB Matched ITRAX log_inc") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_inc_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.1_Matched_log_inc_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Depth plot - matched log_inc Z-scores
HER42PB_Matched_log_inc_depth_Z <- ggplot(Matched_reorder_log_inc_Z, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
  geom_lineh(colour = "black") +
  #geom_point(shape = 21, fill = "white", color = "black", size = 2) +
  xlim(-2.5, 2.5) +
  #scale_x_reverse() +
  scale_y_reverse() +
  ylim(410, 0) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Depth (cm)", y = "Ln(E/inc.) Z-score") +
  labs(title = "HER42PB Matched ITRAX log_inc [Z-score]") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_inc_depth_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.2_Matched_log_inc_depth_Z.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# ICPMS matched plots

# Select site & order that plots appear in
Matched_reorder_log_icpms <- ACE_matched_xrf_icp_log_inc_long_icp %>% # define order plots appear in
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K_ICP", "Ca_ICP", "Ti_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP", "Mo_coh")) 
Matched_reorder_log_icpms_Z <- ACE_matched_xrf_icp_log_inc_Z_long_icp %>% 
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K_ICP", "Ca_ICP", "Ti_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP", "Mo_coh")) 

# Depth - matched log
HER42PB_Matched_log_icpms_depth <- ggplot(Matched_reorder_log_icpms, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
  geom_lineh(colour = "darkgreen") +
  geom_point(shape = 21, fill = "white", color = "darkgreen", size = 2) +
  #scale_x_reverse() +
  scale_y_reverse() +
  ylim(410, 0) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Depth (cm)", y = "log (ppm)") +
  labs(title = "HER42PB Matched log ICPMS") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_icpms_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.3_Matched_log_icpms_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Depth plot - matched log Z-scores + CONISS
HER42PB_Matched_log_icpms_depth_Z <- ggplot(Matched_reorder_log_icpms_Z, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
  geom_lineh(colour = "darkgreen") +
  geom_point(shape = 21, fill = "white", color = "darkgreen", size = 2) +
  xlim(-2.5, 2.5) +
  #scale_x_reverse() +
  scale_y_reverse() +
  ylim(410, 0) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Depth (cm)", y = "log (ppm) [Z-score]") +
  labs(title = "HER42PB Matched log ICPMS [Z-score]") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_icpms_depth_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.4_Matched_log_icpms_depth_Z.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Age plot - log_inc
HER42PB_Matched_log_inc_age <- ggplot(Matched_reorder_log_inc, aes(x = SH20_mean_age, y = value)) +
  geom_line(colour = "darkgrey") +
  #geom_point(shape =".") +
  #scale_y_reverse() +
  scale_x_reverse() +
  xlim(5000, -500) +
  facet_geochem_grid(vars(param)) +
  # scale_colour_manual(values = c("blue", "black"))
  labs(x = "Age (cal a BP)", y = "Ln(E/inc.)") +
  labs(title = "HER42PB ITRAX log_inc") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_inc_age
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.5_Matched_log_inc_age.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")

# Age plot - log_inc Z-scores
HER42PB_Matched_log_inc_age_Z <- ggplot(Matched_reorder_log_inc_Z, aes(x = SH20_mean_age, y = value)) +
  geom_line(colour = "darkgrey") +
  #geom_point() +
  #scale_y_reverse() +
  ylim(-2.5, 2.5) +
  scale_x_reverse() +
  xlim(5000, -500) +
  facet_geochem_grid(vars(param)) +
  # scale_colour_manual(values = c("blue", "black"))
  labs(x = "Age (cal a BP)", y = "Ln(E/inc.) Z-score") +
  labs(title = "HER42PB ITRAX log_inc [Z-score]") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_inc_age_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.6_Matched_log_inc_age_Z.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")


# Age plot - log icpms
HER42PB_Matched_log_icpms_age <- ggplot(Matched_reorder_log_icpms, aes(x = SH20_mean_age, y = value)) +
  geom_line(colour = "darkgreen") +
  geom_point(shape = 21, fill = "white", color = "darkgreen", size = 2) +
  #scale_y_reverse() +
  scale_x_reverse() +
  xlim(5000, -500) +
  facet_geochem_grid(vars(param)) +
  # scale_colour_manual(values = c("blue", "black"))
  labs(x = "Age (cal a BP)", y = "log (ppm)") +
  labs(title = "HER42PB log ICPMS") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_icpms_age
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.7_Matched_log_icpms_age.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")

# Age plot - log ICPMS Z-scores
HER42PB_Matched_log_icpms_age_Z <- ggplot(Matched_reorder_log_icpms_Z, aes(x = SH20_mean_age, y = value)) +
  geom_line(colour = "darkgreen") +
  geom_point(shape = 21, fill = "white", color = "darkgreen", size = 2) +
  #scale_y_reverse() +
  ylim(-2.5, 2.5) +
  scale_x_reverse() +
  xlim(5000, -500) +
  facet_geochem_grid(vars(param)) +
  # scale_colour_manual(values = c("blue", "black"))
  labs(x = "Age (cal a BP)", y = "log (ppm) Z-score") +
  labs(title = "HER42PB log ICPMS [Z-score]") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
HER42PB_Matched_log_icpms_age_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Matched/HER42PB/Fig4.8_Matched_log_icpms_age_Z.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")


# ITRAX log_inc -----------------------------------------------------------

# Depth plots
library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Select site & order that plots appear in
ITRAX_reorder_log_inc <- ACE_xrf_log_inc_long %>% # define order plots appear in
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")) 
ITRAX_reorder_log_inc_Z <- ACE_xrf_log_inc_Z_long %>% 
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")) 

# Depth - log_inc
p1_ITRAX_log_inc_depth <- ggplot(ITRAX_reorder_log_inc, aes(x = value, y = depth)) +
  geom_lineh(colour = "lightgrey") +
  #geom_point(shape = ".") +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Depth (cm)", y = "Ln(E/inc.)") +
  labs(title = "HER42PB ITRAX log_inc") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p1_ITRAX_log_inc_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.1_ITRAX_log_inc_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Depth plot - log_inc Z-scores + CONISS
p2_ITRAX_log_inc_depth_Z <- ggplot(ITRAX_reorder_log_inc_Z, 
                                   aes(x = value,y = depth, xmin = -2, xmax = 2)) +
  geom_col_segsh(colour = "lightgrey") +
  #xlim(-2.5, 2.5) +
  ylim(0, 410) +
  #scale_x_reverse() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Depth (cm)", y = "Ln(E/inc.) Z-score") +
  labs(title = "HER42PB ITRAX log_inc [Z-score]") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p2_ITRAX_log_inc_depth_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.2_ITRAX__log_inc_depth_Z.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Add CONISS & save
coniss_ITRAX_log_inc_Z <- ITRAX_reorder_log_inc_Z %>%
  nested_data(qualifiers = c(SH20_mean_age, depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()
p3_ITRAX__log_inc_Z_depth_coniss <- p2_ITRAX_log_inc_depth_Z +
  layer_dendrogram(coniss_ITRAX_log_inc_Z, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss_ITRAX_log_inc_Z, aes(y = depth))
p3_ITRAX__log_inc_Z_depth_coniss
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.3_ITRAX__log_inc_depth_Z_coniss.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Age plot - log_inc
p4_ITRAX_log_inc_age <- ggplot(ITRAX_reorder_log_inc, aes(x = SH20_mean_age, y = value)) +
  geom_line(colour = "lightgrey") +
  #geom_point(shape =".") +
  #scale_y_reverse() +
  scale_x_reverse() +
  xlim(5000, -500) +
  facet_geochem_grid(vars(param)) +
  # scale_colour_manual(values = c("blue", "black"))
  labs(x = "Age (cal a BP)", y = "Ln(E/inc.)") +
  labs(title = "HER42PB ITRAX log_inc") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p4_ITRAX_log_inc_age
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.4_ITRAX_log_inc_age.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")

# Age plot - log_inc Z-scores + CONISS
p5_ITRAX_log_inc_age_Z <- ggplot(ITRAX_reorder_log_inc_Z, aes(x = SH20_mean_age, y = value)) +
  geom_line(colour = "lightgrey") +
  #geom_point() +
  #scale_y_reverse() +
  ylim(-2.5, 2.5) +
  scale_x_reverse() +
  xlim(5000, -500) +
  facet_geochem_grid(vars(param)) +
  # scale_colour_manual(values = c("blue", "black"))
  labs(x = "Age (cal a BP)", y = "Ln(E/inc.) Z-score") +
  labs(title = "HER42PB ITRAX log_inc [Z-score]") +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() # remove gridlines from palaeo theme for clarity
  )
p5_ITRAX_log_inc_age_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.5_ITRAX_log_inc_age_Z.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")

# Add CONISS & save
coniss_ITRAX_log_inc_age_Z <- ITRAX_reorder_log_inc_Z %>%
  nested_data(qualifiers = c(SH20_mean_age, depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()
p6_ITRAX_log_inc_age_Z_coniss <- p5_ITRAX_log_inc_age_Z +
  layer_dendrogram(coniss_ITRAX_log_inc_age_Z, aes(x = SH20_mean_age), param = "CONISS") +
  layer_zone_boundaries(coniss_ITRAX_log_inc_age_Z, aes(x = SH20_mean_age))
p6_ITRAX_log_inc_age_Z_coniss
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.6_ITRAX_log_inc_age_Z_coniss.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")


# ITRAX with smoothing---------------

library(tidypaleo)
theme_set(theme_bw(base_size=12))

# Pivot to long format

# HER42PB
H42PB_xrf <- ACE_xrf_log_inc %>%
  filter(Site == "HER42PB") %>% 
  select(all_of(acf_icp_Elements_key), MSE,  depth, SH20_mean_age, label) %>% 
  select(-c(Zn, Mn)) # remove elements don't want to plot
# plotting the running mean stratigraphically - and add onto XRF plot for each Site
# with ACF elements only
H42PB_xrf_smooth <- full_join(y = H42PB_xrf %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_key)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = ACE_xrf_log_inc %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
)

Fig4a <- H42PB_xrf_smooth %>% 
  select(all_of(acf_icp_Elements_key2), depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "SH20_mean_age", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c(all_of(acf_icp_Elements_key2)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth, ymin= 0, ymax = 410)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Ln (Element/inc)", y = "Depth (cm)") +
  ggtitle("HER42PB")
print(Fig4a)

library(rioja)
library(repr)
library(patchwork)

H42PB_xrf_long <- H42PB_xrf_smooth %>% 
  select(all_of(acf_icp_Elements_key), depth, label, type) %>%
  select(-c(Zn, Mn, Rb)) %>%  # remove elemenst don't want to plot
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c(all_of(acf_icp_Elements_key)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  )
# adding dengrogram of layers with CONISS constrained cluster analysis from Rioja package
HER42PB_coniss <- H42PB_xrf_long %>%
  nested_data(qualifiers = c(depth), key = elements, value = peakarea, trans = scale) %>%
  nested_chclust_coniss()

HER42PB_coniss_plot <- ggplot() +
  layer_dendrogram(coniss, aes(y = depth, ymin= 0, ymax = 410)) +
  layer_zone_boundaries(coniss, aes(y = depth)) +
  scale_y_reverse() +
  facet_geochem_gridh(vars("CONISS"),
                      units = c("CONISS" = "SS")) +
  xlab(expression (SSquares))  +
  #labs(x = NULL) +
  theme(text=element_text(size=16, face = "plain"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,0), "cm")
  )


# Summary acf element plot vs depth
Fig4b <- H42PB_xrf %>% 
  #filter(qc == TRUE) %>% 
  select(any_of(acf_icp_Elements_key), depth, label) %>%
  pivot_longer(any_of(acf_icp_Elements_key), names_to = "param", values_to = "element") %>%
  filter(param %in% acf_icp_Elements_key) %>%
  mutate(param = fct_relevel(param, acfElementsList_max)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) + #don't use for ITRAX - too many datapoints 
  #geom_lineh(size = 0.5) + #this will make a single black line plot
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 5)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Peak area [cps]", y = "Depth [cm]") +
  ggtitle("ACE Composite ITRAX dataset: cps (ACF-filtered elements) max ACF >0.5")
print(Fig3.9)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section3/Figures/Fig3.9_ACFmax_key_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Ti - uXRF and ICPMS 

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
HER42PB_xrf_icp_age <- ggplot(data = ACE_xrf_log_inc, aes(SH20_mean_age, Ti)) + 
  geom_line(colour = "#74C476", alpha = 1, linewidth = 0.5) +
  geom_point(data = ACE_xrf_icp , aes(SH20_mean_age, Ti_ICP), colour = "#006D2C", fill = "blue", size = 0.5) + 
  geom_line(data = ACE_xrf_icp , aes(SH20_mean_age, Ti_ICP), colour = "#006D2C", linewidth = 0.5) +
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,5000,500), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Ardley Lake (ARD): Ca", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black"))+
  #geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), 
  #           linewidth = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
HER42PB_age

YAN_age <- ggplot(data = YAN_Ln_Ti_norm.Z_nonmarine, aes(SH20_mean_age, Ca)) + 
  geom_line(colour = "#969696", alpha = 1, size = 0.5) +
  geom_point(data = db_YAN.Z, aes(SH20_mean_age, Ca_Ti), colour = "#252525", fill = "#252525", size = 0.5) + 
  geom_line(data = db_YAN.Z, aes(SH20_mean_age, Ca_Ti), colour = "#252525", size = 0.5) + 
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Yanou Lake (YAN): Ca", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black")) +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
YAN_age

# align axes exactly for both graphs
p1.1 <- ARD_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
p1.2 <- YAN_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
pp1 <- list(p1.1, p1.2)
plot_grid(plotlist = pp1, ncol = 1, nrow = 2, align = "v")
ggsave("Figures/Fig 1C_ITRAX-XRF_Age_final_CONISS_Ca.pdf",
       height = c(6), width = c(10), dpi = 600, units = "cm")





# Set plotting parameters & universal plot size for base R  --------------------

# cm graph plot size *from* cm *to* inches:
plotinch_x <- 2 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch_y <- 7 / cm(1) # -> 8 cm  is  3.149606 inches
# or use aspect ratio
aspect_ratio <- 2.5
options(repr.plot.width=plotinch_x, repr.plot.height=plotinch_y)