# Figure 2 - clr plots

# Set up -----------------------------------------------------------------------

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
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis' )
lapply(packages, library, character.only=TRUE)
library()
options(scipen = 999)
# see "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/ACE/ACE_matching_cps.R" 
# for all matching, element correlation & linear modelling and plot code

# Make ACE matched dataset - import matched itrax-ICPMS datafiles from each site -------------------

ACE_matched_BI10 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/BI10_xrf_icp_matched_clr.csv", 
                             col_names = TRUE, skip = 0)
ACE_matched_HER42PB <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/HER42PB_xrf_icp_matched_clr.csv", 
                                col_names = TRUE, skip = 0)
ACE_matched_KER1 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/KER1_xrf_icp_matched_clr.csv", 
                             col_names = TRUE, skip = 0)
ACE_matched_KER3 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/KER3_xrf_icp_matched_clr.csv", 
                             col_names = TRUE, skip = 0)
ACE_matched_PB1 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/PB1_xrf_icp_matched_clr.csv", 
                            col_names = TRUE, skip = 0)
ACE_matched_POB4 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/POB4_xrf_icp_matched_clr.csv", 
                             col_names = TRUE, skip = 0)

# Combine matched output from each site into ACE matched dataset
ACE_matched_clr <- bind_rows(ACE_matched_BI10, 
                             ACE_matched_HER42PB, 
                             ACE_matched_KER1, 
                             ACE_matched_KER3, 
                             ACE_matched_PB1, 
                             ACE_matched_POB4) %>% 
  print()
write.csv(ACE_matched_clr,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/ACE_xrf_icp_matched_clr.csv", row.names = FALSE)


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

# key elements for Figure 2a
xrf_icp_Elements_key_Fig2 <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP", 
                               "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "dry_mass_pc")

# Import existing ACE cps matched file - alternate method for making clr dataset from existing ACE cps dataset ----------------------------------------------

# Matching code for each site and for are in "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/Sites/..."

# see "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/ACE/ACE_matching_clr.R" 
# for code that creates of matched ACE dataset from individual site matched data

ACE_xrf_icp_matched <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Input/clr/ACE_xrf_icp_matched_cps.csv")
is.na(ACE_xrf_icp_matched)<-sapply(ACE_xrf_icp_matched, is.infinite) # replace any infinite values with NA

# Replace zeros with half minimum value for each element 
ACE_all0 <- ACE_xrf_icp_matched %>% 
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

# Make clr file of while dataset after replacing zeros with half min value and using using composition package 
ACE_all_clr_itrax <- ACE_all0 %>% 
  select(all_of(acf_icp_Elements_min)) %>%
  select(-c(Mo_inc, Mo_coh)) %>% 
    clr() %>% # make centred log ratio dataframe
  as_tibble()

ACE_all_clr_itrax_sd <- ACE_all_clr_itrax %>% 
  mutate(across(where(is.numeric), ~ .x * 0.05)) %>% 
  rename_with(.fn = function(.x){paste0(.x,"_sd")})

ACE_all_clr_icp <- ACE_all0 %>% 
  select(all_of(icp_Elements_min)) %>%
  clr() %>% # make centred log ratio dataframe
  as_tibble()

ACE_all_clr_icp_sd <- ACE_all_clr_icp %>% 
  mutate(across(where(is.numeric), ~ .x * 0.05)) %>% 
  rename_with(.fn = function(.x){paste0(.x,"_sd")})

# Make other parameters / variables 
ACE_all_text <- ACE_all0 %>% 
  select(Location:MSE)
ACE_all_rest <- ACE_all0 %>% 
  select(Mo_inc, Mo_coh, Mo_inc_sd, Mo_coh_sd, coh_inc, coh_inc_sd,
         dry_mass_pc, dry_mass_err)

# Bind into single file & define dataset to use for linear modelling
ACE_LM1 <- bind_cols(ACE_all_text, ACE_all_clr_itrax, ACE_all_clr_itrax_sd, 
                     ACE_all_clr_icp, ACE_all_clr_icp_sd, ACE_all_rest) %>%
  filter(!Site =="POB4") %>% #remove POB4 data
  print()
write.csv(ACE_LM1,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/clr/ACE_xrf_icp_matched_noPOB4_clr.csv", row.names = FALSE)

# Convert Site to use as a grouping variable
ACE_LM1$Site <- as.factor(ACE_LM1$Site)

# Correlation matrices --------------------------------------------------

# Fig 2a - ITRAX & ICP correlation matrix - key elements reduced
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_key_Fig2], method = c("everything", "pearson"),
       size = 6, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 6)
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2a_Corr_matrix_key_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation & data plots -----------------------------------------------------------

# Define text for titles & labels for plotting
XRF_title <- " clr XRF-CS"
ICP_title <- " clr ICPMS"
correlation_title <- "clr Correlation"
method_title <- "clr"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Individual element plots -----------------------------------------------------------

# Only Ti shown below as shown in Fig 2 - full version of code is:
# "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/ACE/ACE_matching_clr.R"

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
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2b_Ti_Linear_reg_clr.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2c_Linear_reg_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

