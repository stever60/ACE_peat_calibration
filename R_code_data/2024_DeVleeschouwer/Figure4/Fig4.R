# Figure 4

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(tidypaleo)
library(dplyr)
library(readr)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(cowplot) # for plotting
library(vegan)
library(rioja)
library(ellipse)  # for PCA and cluster
library(factoextra) # for PCA and cluster
library(reshape2)
library(GGally)
library(ggsci)
library(ggdendro)
library(dendextend)
library(dynamicTreeCut)
library(colorspace)
library(cluster)
library(magrittr) # for piping %>%
library(mgcv)
library(gtable)
library(repr)
library(bestNormalize)
library(sjmisc)
library(chemometrics)
library(compositions)
#colour palettes
library(ggsci) #for npg etc
library(wesanderson) 
library(viridis)        
library(RColorBrewer)

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


# Set plotting parameters & universal plot size for base R  --------------------
 
# cm graph plot size *from* cm *to* inches:
plotinch_x <- 2 / cm(1) # -> 8 cm  is  3.149606 inches
plotinch_y <- 7 / cm(1) # -> 8 cm  is  3.149606 inches
# or use aspect ratio
aspect_ratio <- 2.5
options(repr.plot.width=plotinch_x, repr.plot.height=plotinch_y)

# Define elements for plotting--------------------------------------------------

# XRF 
# elements defined by ITRAX acf and matched to Francois ICPMS element list
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
                      "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")
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
mscl_param  <- c("Den1_SAT", "MS1_SAT", "DCMS1_SAT", "Imp_SAT", "FP_SAT", "RES_SAT")

# Subsample
subsample_param <- c("Water_Content_pc", "Dry_mass",	"LOI550",	"C_content",	
               "Wet_density_g_cm3", "Dry_density_g_cm3",	"DMAR_g_cm2_yr")



# Import datasets, standardise & centre (Z-scores) -----------------------------
# cps ---------------------------------------------------------------------

ACE_xrf_cps <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_inc.csv") %>% 
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc)
ACE_xrf_cps
write.csv(ACE_xrf_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_xrf_cps.csv", row.names = FALSE)

# cps Z-scores
ACE_xrf_cps <- ACE_xrf_cps 
ACE_xrf_cps[, acf_icp_Elements_key] <- scale(ACE_xrf_cps[, acf_icp_Elements_key], center = TRUE, scale = TRUE)
ACE_xrf_cps
write.csv(ACE_xrf_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_xrf_cps_Z.csv", row.names = FALSE)

# cps - convert to long format for facet plotting
acf_icp_Elements_key3 <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh") #no Zn, Mn for plotting
ACE_xrf_cps_long <- ACE_xrf_cps %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_cps_long

ACE_xrf_cps_long <- ACE_xrf_cps %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_cps_long

# log_inc -----------------------------------------------------------------

ACE_xrf_log_inc <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_inc.csv") %>% 
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc) %>%
  mutate_at(vars(all_of(acf_icp_Elements_key)), ## Replace zeros with half minimum value to allow linear modelling to work
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% # Recommended Bertrand et al. (2024) - retains dataframe structure
  mutate_at(vars(all_of(acf_icp_Elements_key)), 
            ~ (. == 0) * min(.[. != 0])/2 + .) %>%
  mutate(across(all_of(acf_icp_Elements_key), log)) %>% #log all xrf and icp data
  mutate(across(acf_icp_Elements_key, ~ ifelse(. <=-10, NA, .)))
ACE_xrf_log_inc
write.csv(ACE_xrf_log_inc,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_xrf_log_inc.csv", row.names = FALSE)

# log_inc Z-scores
ACE_xrf_log_inc_Z <- ACE_xrf_log_inc 
ACE_xrf_log_inc_Z[, acf_icp_Elements_key] <- scale(ACE_xrf_log_inc[, acf_icp_Elements_key], center = TRUE, scale = TRUE)
ACE_xrf_log_inc_Z
write.csv(ACE_xrf_log_inc_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_xrf_log_inc_Z.csv", row.names = FALSE)

# log_inc - convert to long format for facet plotting
acf_icp_Elements_key3 <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh") #no Zn, Mn for plotting
ACE_xrf_log_inc_long <- ACE_xrf_log_inc %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_log_inc_long

ACE_xrf_log_inc_Z_long <- ACE_xrf_log_inc_Z %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_log_inc_Z_long

# clr ---------------------------------------------------------------------

ACE_xrf_clr <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_ITRAX_qc_acf_clr.csv") %>% 
  select(Location:MSE, all_of(acf_icp_Elements_key), Total_scatter, inc_coh, coh_inc)
ACE_xrf_clr
write.csv(ACE_xrf_cps,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_xrf_clr.csv", row.names = FALSE)

# clr Z-scores - should be the same as already centred
ACE_xrf_clr_Z <- ACE_xrf_clr
ACE_xrf_clr_Z[, acf_icp_Elements_key] <- scale(ACE_xrf_clr[, acf_icp_Elements_key], center = TRUE, scale = TRUE)
ACE_xrf_clr_Z
write.csv(ACE_xrf_clr_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_xrf_clr_Z.csv", row.names = FALSE)

# clr - convert to long format for facet plotting
acf_icp_Elements_key3 <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh") #no Zn, Mn for plotting
ACE_xrf_clr_long <- ACE_xrf_clr %>% 
  select(c(all_of(acf_icp_Elements_key3), Site, depth, SH20_age)) %>%
  pivot_longer(c(`acf_icp_Elements_key3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_xrf_clr_long

# Matched XRF & ICPMS log_inc dataset used in Figure 2  ----------------------

ACE_matched_xrf_icp_log_inc <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_matched_xrf_icp_log_inc.csv") 
is.na(ACE_matched_xrf_icp_log_inc)<-sapply(ACE_matched_xrf_icp_log_inc, is.infinite) # replace any infinite values with NA
ACE_matched_xrf_icp_log_inc
write.csv(ACE_matched_xrf_icp_log_inc,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_matched_xrf_icp_log_inc.csv", row.names = FALSE)

# Standardise and centre dataframe - Z-scores
ACE_matched_xrf_icp_log_inc_Z <- ACE_matched_xrf_icp_log_inc 
ACE_matched_xrf_icp_log_inc_Z[, xrf_icp_elements] <- scale(ACE_matched_xrf_icp_log_inc[, xrf_icp_elements], center = TRUE, scale = TRUE)
ACE_matched_xrf_icp_log_inc_Z
write.csv(ACE_matched_xrf_icp_log_inc_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_matched_xrf_icp_log_inc_Z.csv", row.names = FALSE)

# matched - convert to long format for facet plotting
xrf_icp_elements3 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP",
                       "Fe", "Fe_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP", "Mo_coh")
xrf_icp_elements3_xrf <- c("K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")
xrf_icp_elements3_icp <- c("K_ICP", "Ca_ICP", "Ti_ICP","Fe_ICP", "Sr_ICP", "Zr_ICP", "Mo_coh")

ACE_matched_xrf_icp_log_inc_long <- ACE_matched_xrf_icp_log_inc %>% 
  select(c(all_of(xrf_icp_elements3), Site, depth, SH20_age)) %>%
  pivot_longer(c(`xrf_icp_elements3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_long

ACE_matched_xrf_icp_log_inc_Z_long <- ACE_matched_xrf_icp_log_inc_Z %>% 
  select(c(all_of(xrf_icp_elements3), Site, depth, SH20_age)) %>%
  pivot_longer(c(`xrf_icp_elements3`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_Z_long

ACE_matched_xrf_icp_log_inc_long_xrf <- ACE_matched_xrf_icp_log_inc %>% 
  select(c(all_of(xrf_icp_elements3_xrf), Site, depth, SH20_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_xrf`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_long_xrf

ACE_matched_xrf_icp_log_inc_Z_long_xrf <- ACE_matched_xrf_icp_log_inc_Z %>% 
  select(c(all_of(xrf_icp_elements3_xrf), Site, depth, SH20_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_xrf`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_Z_long_xrf

ACE_matched_xrf_icp_log_inc_long_icp <- ACE_matched_xrf_icp_log_inc %>% 
  select(c(all_of(xrf_icp_elements3_icp), Site, depth, SH20_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_icp`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_long_icp

ACE_matched_xrf_icp_log_inc_Z_long_icp <- ACE_matched_xrf_icp_log_inc_Z %>% 
  select(c(all_of(xrf_icp_elements3_icp), Site, depth, SH20_age)) %>%
  pivot_longer(c(`xrf_icp_elements3_icp`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_matched_xrf_icp_log_inc_Z_long_icp


# MSCL data ---------------------------------------------------------------

HER_MSCL <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/HER_MSCL_comp.csv")
HER_MSCL
# Standardise and centre dataframe - Z-scores
HER_MSCL_Z <- HER_MSCL
HER_MSCL_Z[, mscl_param] <- scale(HER_MSCL[, mscl_param], center = TRUE, scale = TRUE)
HER_MSCL_Z
write.csv(HER_MSCL_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER_MSCL_Comp_Z.csv", row.names = FALSE)

# Subsample data ----------------------------------------------------------

ACE_DEN <- read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Input/ACE_Den_comp.csv")
ACE_DEN
# Standardise and centre dataframe - Z-scores
ACE_DEN_Z <- ACE_DEN
ACE_DEN_Z[, subsample_param] <- scale(ACE_DEN[, subsample_param], center = TRUE, scale = TRUE)
ACE_DEN_Z
write.csv(ACE_DEN_Z,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_Den_comp_Z.csv", row.names = FALSE)

ACE_DEN_long <- ACE_DEN %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_DEN_long

ACE_DEN_long_Z <- ACE_DEN_Z %>% 
  select(c(all_of(subsample_param), Site, depth, SH20_age)) %>%
  pivot_longer(c(`subsample_param`), names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth)
ACE_DEN_long_Z


library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# Depth - log_inc
p1_ACE_DEN__depth <- ggplot(ACE_DEN_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "white", color = "black", size = 1) +
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
p1_ACE_DEN__depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.1_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")




# Plots -------------------------------------------------------------------
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
  nested_data(qualifiers = c(SH20_age, depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()
p3_ITRAX__log_inc_Z_depth_coniss <- p2_ITRAX_log_inc_depth_Z +
  layer_dendrogram(coniss_ITRAX_log_inc_Z, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss_ITRAX_log_inc_Z, aes(y = depth))
p3_ITRAX__log_inc_Z_depth_coniss
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.3_ITRAX__log_inc_depth_Z_coniss.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Age plot - log_inc
p4_ITRAX_log_inc_age <- ggplot(ITRAX_reorder_log_inc, aes(x = SH20_age, y = value)) +
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
p5_ITRAX_log_inc_age_Z <- ggplot(ITRAX_reorder_log_inc_Z, aes(x = SH20_age, y = value)) +
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
  nested_data(qualifiers = c(SH20_age, depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()
p6_ITRAX_log_inc_age_Z_coniss <- p5_ITRAX_log_inc_age_Z +
  layer_dendrogram(coniss_ITRAX_log_inc_age_Z, aes(x = SH20_age), param = "CONISS") +
  layer_zone_boundaries(coniss_ITRAX_log_inc_age_Z, aes(x = SH20_age))
p6_ITRAX_log_inc_age_Z_coniss
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.6_ITRAX_log_inc_age_Z_coniss.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")


# Matched log_inc xrf & log icpms  -----------------------------------------------------------

# Depth plots
library(tidypaleo) #https://cran.r-project.org/web/packages/tidypaleo/vignettes/strat_diagrams.html
theme_set(theme_paleo(12)) #theme_paleo

# xrf matched
# Select site & order that plots appear in
Matched_reorder_log_inc <- ACE_matched_xrf_icp_log_inc_long_xrf %>% # define order plots appear in
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")) 
Matched_reorder_log_inc_Z <- ACE_matched_xrf_icp_log_inc_Z_long_xrf %>% 
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K", "Ca", "Ti", "Fe", "Sr", "Zr", "Mo_coh")) 

# Depth - matched log_inc
p1_Matched_log_inc_depth <- ggplot(Matched_reorder_log_inc, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
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
p1_Matched_log_inc_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.1_Matched_log_inc_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Depth plot - matched log_inc Z-scores
p2_Matched_log_inc_depth_Z <- ggplot(Matched_reorder_log_inc_Z, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
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
p2_Matched_log_inc_depth_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.2_Matched__log_inc_depth_Z.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Age plot - log_inc
p4_Matched_log_inc_age <- ggplot(Matched_reorder_log_inc, aes(x = SH20_age, y = value)) +
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
p4_Matched_log_inc_age
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.4_Matched_log_inc_age.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")

# Age plot - log_inc Z-scores
p5_Matched_log_inc_age_Z <- ggplot(Matched_reorder_log_inc_Z, aes(x = SH20_age, y = value)) +
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
p5_Matched_log_inc_age_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.5_Matched_log_inc_age_Z.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")

# icpms matched

# Select site & order that plots appear in
Matched_reorder_log_icpms <- ACE_matched_xrf_icp_log_inc_long_icp %>% # define order plots appear in
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K_ICP", "Ca_ICP", "Ti_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP", "Mo_coh")) 
Matched_reorder_log_icpms_Z <- ACE_matched_xrf_icp_log_inc_Z_long_icp %>% 
  filter(Site =="HER42PB") %>% 
  mutate(param = fct_relevel(param,"K_ICP", "Ca_ICP", "Ti_ICP", "Fe_ICP", "Sr_ICP", "Zr_ICP", "Mo_coh")) 

# Depth - matched log
p1_Matched_log_icpms_depth <- ggplot(Matched_reorder_log_icpms, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
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
p1_Matched_log_icpms_depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.1_Matched_log_icpms_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Depth plot - matched log Z-scores + CONISS
p2_Matched_log_icpms_depth_Z <- ggplot(Matched_reorder_log_icpms_Z, aes(x = value, y = depth, ymin= 0, ymax = 410)) +
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
p2_Matched_log_icpms_depth_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.2_Matched_log_icpms_depth_Z.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# Age plot - log icpms
p4_Matched_log_icpms_age <- ggplot(Matched_reorder_log_icpms, aes(x = SH20_age, y = value)) +
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
p4_Matched_log_icpms_age
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.4_Matched_log_icpms_age.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")

# Age plot - log ICPMS Z-scores
p5_Matched_log_icpms_age_Z <- ggplot(Matched_reorder_log_icpms_Z, aes(x = SH20_age, y = value)) +
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
p5_Matched_log_icpms_age_Z
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.5_Matched_log_icpms_age_Z.pdf",
       height = c(36), width = c(24), dpi = 600, units = "cm")


p8 <- p2_ITRAX_log_inc_depth_Z + p2_Matched_log_inc_depth_Z
  
p8 



# Calibration  ----------------------------------------------------------

# Ti

# ICPMS data - change errors to Table 2 (previously +/-5%)
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
#calculate error in log inc ITRAX data
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
write.csv(ACE_ICP_ppm,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_matched_xrf_cps.csv", row.names = FALSE)

# Input ICPMS data for comparison / overlay
ACE_ICP_ppm_Ti <-read_csv("Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_matched_xrf_icp_cps.csv") %>%
  select(Site, depth, SH20_age, Ti, Ti_sd, Ti_ICP, Ti_ICP_sd) %>%
  filter(Site == "HER42PB")  
ACE_ICP_ppm_Ti

# Ti XRF-CS as ppm with RMSE errors
# Ti Ln equation: y = 11+0.67 x  where y = Ln(Ti_ICP) and x = Ln(Ti_XRF) & RMSE = 0.663
 
ACE_xrf_calib_Ti <- ACE_xrf_log_inc %>%
  select(Site, depth, SH20_age, Ti) %>%
  mutate(Ti_convert = 11+0.67*Ti) %>%
  mutate(Ti_ppm = exp(Ti_convert)) %>% 
  mutate(Ti_upper = Ti+0.663) %>% 
  mutate(Ti_convert_upper = 11+0.67*Ti_upper) %>%
  mutate(Ti_upper_RMSE = exp(Ti_convert_upper)) %>% 
  mutate(Ti_lower = Ti-0.663) %>% 
  mutate(Ti_convert_lower = 11+0.67*Ti_lower) %>%
  mutate(Ti_lower_RMSE = exp(Ti_convert_lower)) %>% 
  select(Site, depth, SH20_age, Ti, Ti_ppm, Ti_lower_RMSE, Ti_upper_RMSE)
ACE_xrf_calib_Ti
write.csv(ACE_xrf_calib_Ti,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/ACE_xrf_calib_Ti.csv", row.names = FALSE)

HER42PB_xrf_calib_Ti <- ACE_xrf_calib_Ti %>% 
  filter(Site == "HER42PB")
write.csv(HER42PB_xrf_calib_Ti,"Papers_R/2024_DeVleeschouwer/Figure4/Data/Output/HER42PB_xrf_calib_Ti.csv", row.names = FALSE)

# Depth plot
p12_Ti <- 
  ggplot() +
  geom_lineh(data = HER42PB_xrf_calib_Ti, aes(x=Ti_lower_RMSE, y=depth), color = "lightblue") +
  geom_lineh(data = HER42PB_xrf_calib_Ti, aes(x=Ti_upper_RMSE, y=depth), color = "lightblue") +
  geom_lineh(data = HER42PB_xrf_calib_Ti, aes(x=Ti_ppm, y=depth), color = "blue") +
  geom_lineh(data = ACE_ICP_ppm_Ti, aes(x=Ti_ICP, y=depth), color = "red") +
  geom_errorbarh(data = ACE_ICP_ppm_Ti, aes(x=Ti_ICP, y=depth, xmin=Ti_ICP-Ti_ICP_sd, xmax=Ti_ICP+Ti_ICP_sd), color = "red", height=0) +
  geom_point(data = ACE_ICP_ppm_Ti, aes(x=Ti_ICP, y=depth), fill = "white", color = "red", shape = 21, size = 2) +
  scale_y_reverse() +
  ylim(410, 0) +
  xlim(0,7000)
p12_Ti
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.12_ICPMS_Ti.pdf",
        height = c(24), width = c(8), dpi = 600, units = "cm")


# Figure 4 - plots 

# Depth -density plots
p1_ACE_DEN__depth <- ggplot(ACE_DEN_long, aes(x = value, y = depth)) +
  geom_lineh(colour = "black") +
  geom_point(shape = 21, fill = "white", color = "black", size = 1) +
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
p1_ACE_DEN__depth
ggsave("Papers_R/2024_DeVleeschouwer/Figure4/Plots/Fig4.1_Subsample_depth.pdf",
       height = c(24), width = c(36), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# EXTRA


# add units to the graph headers
#facet_geochem_gridh(
#  vars(param),
#  units = c("Shard_Counts" = "n", "Total_Shards" = "no. per g DM", "CONISS" = "SS")) +

# ITRAX with smoothing---------------

library(tidypaleo)
theme_set(theme_bw(base_size=12))

# Pivot to long format

# HER42PB
H42PB_xrf <- ACE_xrf_log_inc %>%
  filter(Site == "HER42PB") %>% 
  select(all_of(acf_icp_Elements_key), MSE,  depth, SH20_age, label) %>% 
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
  tidyr::pivot_longer(!c("depth", "SH20_age", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
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
HER42PB_xrf_icp_age <- ggplot(data = ACE_xrf_log_inc, aes(SH20_age, Ti)) + 
  geom_line(colour = "#74C476", alpha = 1, linewidth = 0.5) +
  geom_point(data = ACE_xrf_icp , aes(SH20_age, Ti_ICP), colour = "#006D2C", fill = "blue", size = 0.5) + 
  geom_line(data = ACE_xrf_icp , aes(SH20_age, Ti_ICP), colour = "#006D2C", linewidth = 0.5) +
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,5000,500), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Ardley Lake (ARD): Ca", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black"))+
  #geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), 
  #           linewidth = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
HER42PB_age

YAN_age <- ggplot(data = YAN_Ln_Ti_norm.Z_nonmarine, aes(SH20_age, Ca)) + 
  geom_line(colour = "#969696", alpha = 1, size = 0.5) +
  geom_point(data = db_YAN.Z, aes(SH20_age, Ca_Ti), colour = "#252525", fill = "#252525", size = 0.5) + 
  geom_line(data = db_YAN.Z, aes(SH20_age, Ca_Ti), colour = "#252525", size = 0.5) + 
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