# Figure 2 - cps plots

# Libraries ---------------------------------------------------------------

library(itraxR)
library(ggplot2)
library(ggpubr) # plotting
library(tidyverse) # all core tidyverse packages
library(tidypaleo) # Dewey Dunnington's ggplot extensions for palaeo-style plots
library(compositions)
library(scales)
library(ggsci) # colour palettes for publication
library(GGally) # for correlation and Prob density matrix plotting
library()
options(scipen = 999)

# Set up ------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()

# see "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/ACE/ACE_matching_cps.R" 
# for all matching, element correlation & linear modelling and plot code


# Make ACE matched dataset - import matched itrax-ICPMS datafiles from each site -------------------

ACE_matched_BI10 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/BI10_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)
ACE_matched_HER42PB <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/HER42PB_xrf_icp_matched_cps.csv", 
                              col_names = TRUE, skip = 0)
ACE_matched_KER1 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/KER1_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)
ACE_matched_KER3 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/KER3_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)
ACE_matched_PB1 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/PB1_xrf_icp_matched_cps.csv", 
                          col_names = TRUE, skip = 0)
ACE_matched_POB4 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/POB4_xrf_icp_matched_cps.csv", 
                           col_names = TRUE, skip = 0)

# Combine matched output from each site into ACE matched dataset
ACE_matched_cps <- bind_rows(ACE_matched_BI10, 
                                 ACE_matched_HER42PB, 
                                 ACE_matched_KER1, 
                                 ACE_matched_KER3, 
                                 ACE_matched_PB1, 
                                 ACE_matched_POB4) %>% 
  print()
write.csv(ACE_matched_cps,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/ACE_xrf_icp_matched_cps.csv", row.names = FALSE)

# Define elements to use ----------------------------------------------------------

# ICP elements as defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")

# Below were defined by autocorrelation (acf) analysis of ACE09 composite ITRAX dataframe analysis file - see other R files for more info

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

# key elements_reduced for more simplified plots 
xrf_icp_Elements_key_reduced <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")

# Import existing ACE cps matched file ----------------------------------------------

ACE_xrf_icp_matched_input <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/cps/ACE_xrf_icp_matched_cps.csv")
is.na(ACE_xrf_icp_matched)<-sapply(ACE_xrf_icp_matched, is.infinite) # replace any infinite values with NA

# ICPMS data - change errors supplied by FDV originally to CRM % values in Table 2 
ACE_xrf_icp_matched <- ACE_xrf_icp_matched_input %>%
  rename(K_ICP_sd_FDV2023 = K_ICP_sd) %>% 
  mutate(K_ICP_sd = K_ICP*0.18) %>%
  relocate(K_ICP_sd, .after = K_ICP) %>%
  rename(Ca_ICP_sd_FDV2023 = Ca_ICP_sd) %>% 
  mutate(Ca_ICP_sd = Ca_ICP*0.13) %>% 
  relocate(Ca_ICP_sd, .after = Ca_ICP) %>%
  rename(Ti_ICP_sd_FDV2023 = Ti_ICP_sd) %>% 
  mutate(Ti_ICP_sd = Ti_ICP*0.14) %>%
  relocate(Ti_ICP_sd, .after = Ti_ICP) %>%
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
  relocate(Sr_ICP_sd, .after = Sr_ICP)
write.csv(ACE_xrf_icp_matched,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/cps/ACE_matched_xrf_icp_cps.csv", row.names = FALSE)

# Define dataset to use for linear modelling
ACE_LM1 <- ACE_xrf_icp_matched %>%
  filter(!Site =="POB4") %>% #remove POB4 data
  print()
write.csv(ACE_LM1,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/cps/ACE_xrf_icp_matched_cps_noPOB4.csv", row.names = FALSE)

# Convert Site to use as a grouping variable
ACE_LM1$Site <- as.factor(ACE_LM1$Site)


# Correlation matrices --------------------------------------------------

# Fig 2a - ITRAX & ICP correlation matrix - key elements reduced
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_key_reduced], method = c("everything", "pearson"),
       size = 10, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 10)
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2a_Corr_matrix_key_reduced_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation & data plots -----------------------------------------------------------

# Define text for titles & labels for plotting
XRF_title <- " cps XRF-CS"
ICP_title <- " cps ICPMS"
correlation_title <- "cps Correlation"
method_title <- "cps"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Individual element plots -----------------------------------------------------------

# Only Ti shown below as shown in Fig 2 - full version of code is:
# "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/ACE/ACE_matching_cps.R"

# Ti ----------------------------------------------------------------------
element_title <- "Ti"

# Fig 2b - Linear regression model with 68% ellipses - all sites
theme_set(theme_classic(10))
Ti_corr1 <- ggscatter(ACE_LM1, x = "Ti", y = "Ti_ICP",
                      color = "Site", palette = palette_set,size = 2, alpha = 1, 
                      rug = TRUE, ellipse = TRUE, ellipse.level = 0.68, ellipse.alpha = 0.1,
                      add = "reg.line", conf.int = TRUE, add.params = list(color = "black", fill = "lightgrey")) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 0), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("ACE (OLS): ", element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=10, face="bold"))
Ti_corr1
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2b_Ti_Linear_reg_cps.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig 2c - Linear regression model - per site
theme_set(theme_classic(10))
Ti_corr_sites <- ggscatter(ACE_LM1, x = "Ti", y = "Ti_ICP", size = 1,
                      color = "Site", palette = palette_set,
                      facet.by = "Site", #scales = "free_x",
                      add = "reg.line", conf.int = TRUE, alpha = 0.5) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "none") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(strip.background = element_blank()) +
  stat_cor(aes(color = Site), method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.1, size = 3) + #label.sep = "\n"
  stat_regline_equation(aes(color = Site), label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.1, size = 3) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
Ti_corr_sites
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2c_Linear_reg_per_site_cps.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")




