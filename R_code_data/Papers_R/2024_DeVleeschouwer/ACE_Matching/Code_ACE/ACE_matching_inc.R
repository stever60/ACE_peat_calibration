# -------------------------------------------------------------------------

# ITRAX-ICPMS Data Matching, Correlation & Regression forinc normlaised dataset

# -------------------------------------------------------------------------

## Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library()
options(scipen = 999)
library(itraxR)
library(tidyverse) # all core tidyverse packages
library(tidypaleo) # Dewey Dunnington's ggplot extensions for palaeo-style plots
library(readr)
library(ggpubr) # plotting
library(GGally) # for correlation and Prob density matrix plotting
library(PeriodicTable)
library(errors)
library(chemometrics)
library(patchwork)
library(forecast) # Use autocorrelation function (acf) and plots to explore noise in a time-series
library(directlabels)
library(broom)
library(performance)
library(lmtest) # linear moedlling parameters
library(ggpmisc)
library(scales)
library(ggsci) # colour palettes for publication
library(cowplot) # graphing arrangement
library(Hmisc) # various stats tests

# Set up ------------------------------------------------------------------

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro 2013
# setwd("/Users/Steve/Dropbox/BAS/Data/R/")
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()

# Import matched itrax-ICPMS datafiles from each site -------------------

ACE_itrax_BI10 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/BI10/inc/BI10_xrf_icp_matched_inc.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_HER42PB <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/HER42PB/inc/HER42PB_xrf_icp_matched_inc.csv", 
                              col_names = TRUE, skip = 0)
ACE_itrax_KER1 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER1/inc/KER1_xrf_icp_matched_inc.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_KER3 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/inc/KER3_xrf_icp_matched_inc.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_PB1 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/PB1/inc/PB1_xrf_icp_matched_inc.csv", 
                          col_names = TRUE, skip = 0)
#ACE_itrax_POB4 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/POB4/inc/POB4_xrf_icp_matched_inc.csv", 
#                           col_names = TRUE, skip = 0)

# Combine matched output from each site
inc_cols <- c('Mo_inc', 'Mo_inc_sd')
ACE_xrf_icp_matched <- bind_rows(ACE_itrax_BI10, 
                                 ACE_itrax_HER42PB, 
                                 ACE_itrax_KER1, 
                                 ACE_itrax_KER3, 
                                 ACE_itrax_PB1,
                                 ACE_itrax_POB4) %>% 
  select(!all_of(inc_cols)) %>% 
  print()
write.csv(ACE_xrf_icp_matched,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_xrf_icp_matched.csv", row.names = FALSE)

# START HERE if ACE matched already written ----------------------------------------------------------

# Element lists

# Defined from ACF analysis of ACE09 composite ITRAX dataframe analysis file:
# https://www.dropbox.com/s/rhtlkp6uwp71ryc/ACE_itrax_Composite.R?dl=0

# ICP elements as defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")

# XRF elements defined by Francois & ITRAX acf
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_coh_sd")

# ICP elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

xrf_icp_Elements_min <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                          "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                          "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                          "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                          "Mo_coh", "coh_inc", "dry_mass_pc")

# define adjacent matched elements for correlation
xrf_icp_Elements_min <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                          "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                          "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                          "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                          "Mo_coh", "coh_inc", "dry_mass_pc")

# define as adjacent matched elements for log transformation
xrf_icp_Elements_min1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                           "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                           "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                           "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                           "Mo_coh")

# define as adjacent matched elements for log transformation
xrf_icp_Elements_min1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                           "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                           "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                           "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                           "Mo_coh")

# define as adjacent matched elements for PCA
xrf_icp_Elements_min_PCA <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                              "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                              "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                              "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                              "Mo_inc", "Mo_coh", "coh_inc", "dry_mass_pc")

xrf_icp_Elements_min_PCA_edited <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                                     "Mn", "Mn_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP",
                                     "coh_inc", "dry_mass_pc")

# key elements
xrf_icp_Elements_key <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP",
                          "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "dry_mass_pc")

# key elements_reduced
xrf_icp_Elements_key_reduced <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")


# Import existing ACE matched file ----------------------------------------------

ACE_xrf_icp_matched <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_xrf_icp_matched.csv")
is.na(ACE_xrf_icp_matched)<-sapply(ACE_xrf_icp_matched, is.infinite) # replace any infinite values with NA

# ACE matched - all sites
ACE_all <- ACE_xrf_icp_matched %>% 
  mutate(across(all_of(c(icp_Elements_min, icp_Elements_min_sd)), ~ ifelse(.x < 0, 0, .x))) %>%  # replace any ICPMS values <0 with zero
  mutate_at(vars(all_of(c(icp_Elements_min, acf_icp_Elements_min))), ## Replace zeros with half minimum value to allow linear modelling to work
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% # Recommended procedure from Bertrand et al. (submitted) - retains dataframe structure
  mutate_at(vars(all_of(c(icp_Elements_min_sd, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% 
  select(Location:MSE, all_of(acf_icp_Elements_min), coh_inc, all_of(acf_icp_Elements_min_sd), coh_inc_sd,
         all_of(icp_Elements_min), dry_mass_pc, all_of(icp_Elements_min_sd), dry_mass_err) %>% 
  print()
write.csv(ACE_all,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_all.csv", row.names = FALSE)
write.csv(ACE_all,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_xrf_icp_matched_inc.csv", row.names = FALSE)

# ACE matched - POB4 removed for regression analysis in cps runs
ACE_no_POB4 <- ACE_all %>% 
  filter(!Site =="POB4") %>% 
  print()
write.csv(ACE_no_POB4,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_no_POB4.csv", row.names = FALSE)

# ACE matched - POB4 & PB1 removed
ACE_no_POB4_PB1 <- ACE_all %>% 
  filter(!Site =="POB4") %>% 
  filter(!Site =="PB1") %>% 
  print()
write.csv(ACE_no_POB4_PB1,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_no_POB4_PB1.csv", row.names = FALSE)

ACE_no_PB1 <- ACE_all %>% # inc and element/inc only 
  filter(!Site =="PB1") %>% 
  print()
write.csv(ACE_no_PB1,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_no_PB1.csv", row.names = FALSE)

# Create summary stats table ----------------------------------------------

ACE_all_stats <- ACE_all %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_all_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_all_stats.csv", row.names = FALSE)

ACE_no_POB4_stats <- ACE_no_POB4 %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_no_POB4_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_no_POB4_stats.csv", row.names = FALSE)

ACE_no_POB4_PB1_stats <- ACE_no_POB4_PB1 %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_no_POB4_PB1_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_no_POB4_PB1_stats.csv", row.names = FALSE)

ACE_no_PB1_stats <- ACE_no_PB1 %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_no_PB1_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_no_PB1_stats.csv", row.names = FALSE)

# Correlation & OLS linear modelling  --------------------------------------------------

# Choose datasets to import - needs ACE_all, ACE_LM1, ACE_LM2 to run
ACE_all <- ACE_all

#ACE_LM1 <- ACE_all
ACE_LM1 <- ACE_no_POB4
#ACE_LM1 <- ACE_no_POB4_PB1
#ACE_LM1 <- ACE_no_PB1

#ACE_LM2 <- ACE_no_POB4
ACE_LM2 <- ACE_no_POB4_PB1
#ACE_LM2 <- ACE_no_PB1

# Convert Site as a grouping variable
ACE_all$Site <- as.factor(ACE_all$Site)
ACE_LM1$Site <- as.factor(ACE_LM1$Site)
ACE_LM2$Site <- as.factor(ACE_LM2$Site)

# Correlation matrices --------------------------------------------------

# ITRAX
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,acf_icp_Elements_min], method = c("everything", "pearson"), 
       size = 5, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 5)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ICP
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,icp_Elements_min], method = c("everything", "pearson"), 
       size = 5, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 5)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_ICP_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_min], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 3)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP key correlations
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_key], method = c("everything", "pearson"), 
       size = 7, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 7)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr_matrix_key.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP  key reduced correlations
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_key_reduced], method = c("everything", "pearson"), 
       size = 12, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 12)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr_matrix_key_reduced.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation & density matrices ----------------------------------------

# *** if the p-value is < 0.001, ** if the p-value is < 0.01, * if the p-value is < 0.05, . if the p-value is < 0.10

# ITRAX
theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = acf_icp_Elements_min, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE XRF-CS: Correlation-density plot - element/inc")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_Corr-den_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ICP
theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = icp_Elements_min, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE ICPMS: Correlation-density plot - element/inc")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_ICP_Corr-den_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP - key elements correlations
theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = xrf_icp_Elements_min, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE XRF-CS & ICPMS: Correlation-density plot - element/inc")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr-den_matrix_key.pdf", 
       height = c(15), width = c(15), dpi = 600, units = "cm")

theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE XRF-CS & ICPMS: Correlation-density plot - element/inc")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr-den_matrix_key.pdf", 
       height = c(15), width = c(15), dpi = 600, units = "cm")

# ITRAX & ICP - key elements correlations and regressions
# Function to plot each site as different colour usign jco scheme
palette_set <- "jco" # or "jco", or "npg", "uchicago" - set up colour scheme for site plots
cor_plot <- function(data, mapping, ...) {
  
  ggplot(data, mapping) + 
    geom_point(size = 0.7) +
    geom_smooth(formula = y~x, method = lm, color = "black") + 
    scale_color_jco()
}

# Run ggpairs - black regression line all - correlation all sites
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, 
        title="ACE XRF-CS & ICPMS: Correlation-density plot - element/inc",
        upper = list(
          mapping = aes(color=Site, palette = palette_set),
          continuous = wrap(ggally_cor, size = 3)
        ),
        lower=list(
          mapping = aes(color=Site, palette = palette_set, alpha = 0.5),
          continuous = cor_plot
        )
)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr-den_matrix_key_reg1.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs - black regression line all - correlation overall
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key,
        title="ACE XRF-CS & ICPMS: Correlation-density plot - element/inc",
        upper = list(
          continuous = wrap(ggally_cor, size = 5)
        ),
        lower=list(
          mapping = aes(color=Site, palette = palette_set, alpha = 0.5),
          continuous = cor_plot
        )
)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr-den_matrix_key_reg2.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs - regression line all sites - correlation all sites
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, upper = list(continuous = wrap("cor", size = 3)),
        #lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        lower = list(continuous = "smooth"),
        aes(color = Site,  # Color by group (cat. variable)
            alpha = 0.5),     # Transparency
        title="ACE XRF-CS & ICPMS: Correlation-density plot - element/inc")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr-den_matrix_key_reg3.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs - no regression - correlation all sites
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=1)),
        #lower = list(continuous = "smooth"),
        aes(color = Site,  # Color by group (cat. variable)
            alpha = 0.5),     # Transparency
        title="ACE XRF-CS & ICPMS: Correlation-density plot - element/inc")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_itrax_ICP_Corr-den_matrix_site_key_plots.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation & data plots -----------------------------------------------------------

# Define titles & labels for plotting
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "K Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Individual element plots -----------------------------------------------------------

# K -----------------------------------------------------------------------
element_title <- "K"
# Fig1 - Correlation plot with 68% CI ellipses 
theme_set(theme_classic(10))
K_corr1 <- ggscatter(ACE_LM1, x = "K", y = "K_ICP",
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
K_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density - all sites
# a) XRF-CS data
theme_set(theme_classic(10))
K_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=K, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
K_XRF_violin_boxplot <- K_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
K_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=K_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
K_ICP_violin_boxplot <- K_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
K_pmain <- ggplot(ACE_all, aes(x = K, y = K_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
K_xdens <- axis_canvas(K_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = K, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
K_ydens <- axis_canvas(K_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = K_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
K_p1 <- insert_xaxis_grob(K_pmain, K_xdens, grid::unit(.2, "null"), position = "top")
K_corr2 <- insert_yaxis_grob(K_p1, K_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(K_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
K_corr3 <- ggscatter(ACE_all, x = "K", y = "K_ICP", size = 1,
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
K_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(K_XRF_violin_boxplot, K_ICP_violin_boxplot, K_corr2, K_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests - K

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$K
y_all.reg <- ACE_all$K_ICP
x_all <- "K"
y_all <- "K_ICP"
x1.reg <- ACE_LM1$K
y1.reg <- ACE_LM1$K_ICP
x1 <- "K"
y1 <- "K_ICP"
x2.reg <- ACE_LM2$K
y2.reg <- ACE_LM2$K_ICP
x2 <- "K"
y2 <- "K_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, K, K_sd, K_ICP, K_ICP_sd, fit, lwr, upr)) %>% 
  print()
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
K_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 4, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 4, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
K_predict <- K_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
K_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
K_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
K_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
K_modf <- fortify(model_1)
K_res <- ggplot(K_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(K_predict,  K_qq.x1, K_res, K_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$K, ACE_all$K_ICP)
ctest_LM1 <- cor.test(ACE_LM1$K, ACE_LM1$K_ICP)
ctest_LM2 <- cor.test(ACE_LM2$K, ACE_LM2$K_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(K, K_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$K, ACE_all$K_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$K, ACE_LM1$K_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$K, ACE_LM2$K_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
# Test normality of one variable (univariate) at a time 
# If the p-value > 0.05 it implies that the distribution of the data 
# are not significantly different from normal distribution - i.e., can assume normality
# < 0.05 means data are not normally distributed
# Normal distribution is hard to achieve in  geochemical datasets!
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
# Use  to test for correlation among residuals
# H0 (null hypothesis) = no correlation among the residuals
# HA (alternative hypothesis) = residuals are autocorrelated
# Value of 2.0 indicates there is no autocorrelation detected in the sample. 
# 0-2 = positive autocorrelation; 2-4 = negative autocorrelation.
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/K/K_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Ca ----------------------------------------------------------------------
element_title <- "Ca"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Ca Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Ca_corr1 <- ggscatter(ACE_LM1, x = "Ca", y = "Ca_ICP",
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
Ca_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2 - Violin + box plots, scatterplot + density, Individual sites
# a) XRF-CS data
theme_set(theme_classic(10))
Ca_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Ca, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Ca_XRF_violin_boxplot <- Ca_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Ca_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Ca_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Ca_ICP_violin_boxplot <- Ca_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Ca_pmain <- ggplot(ACE_all, aes(x = Ca, y = Ca_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Ca_xdens <- axis_canvas(Ca_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Ca, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Ca_ydens <- axis_canvas(Ca_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Ca_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Ca_p1 <- insert_xaxis_grob(Ca_pmain, Ca_xdens, grid::unit(.2, "null"), position = "top")
Ca_corr2 <- insert_yaxis_grob(Ca_p1, Ca_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Ca_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Ca_corr3 <- ggscatter(ACE_all, x = "Ca", y = "Ca_ICP", size = 1,
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
Ca_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Ca_XRF_violin_boxplot, Ca_ICP_violin_boxplot, Ca_corr2, Ca_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - OLS Linear Model

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Ca
y_all.reg <- ACE_all$Ca_ICP
x_all <- "Ca"
y_all <- "Ca_ICP"
x1.reg <- ACE_LM1$Ca
y1.reg <- ACE_LM1$Ca_ICP
x1 <- "Ca"
y1 <- "Ca_ICP"
x2.reg <- ACE_LM2$Ca
y2.reg <- ACE_LM2$Ca_ICP
x2 <- "Ca"
y2 <- "Ca_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Ca, Ca_sd, Ca_ICP, Ca_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Ca_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Ca_predict <- Ca_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ca_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Ca_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ca_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Ca_modf <- fortify(model_1)
Ca_res <- ggplot(Ca_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Ca_predict,  Ca_qq.x1, Ca_res, Ca_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Ca, ACE_all$Ca_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Ca, ACE_LM1$Ca_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Ca, ACE_LM2$Ca_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Ca, Ca_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Ca, ACE_all$Ca_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Ca, ACE_LM1$Ca_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Ca, ACE_LM2$Ca_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ca/Ca_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Ti ----------------------------------------------------------------------
element_title <- "Ti"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Ti Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Ti_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Ti, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Ti_XRF_violin_boxplot <- Ti_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Ti_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Ti_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Ti_ICP_violin_boxplot <- Ti_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Ti_pmain <- ggplot(ACE_all, aes(x = Ti, y = Ti_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Ti_xdens <- axis_canvas(Ti_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Ti, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Ti_ydens <- axis_canvas(Ti_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Ti_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Ti_p1 <- insert_xaxis_grob(Ti_pmain, Ti_xdens, grid::unit(.2, "null"), position = "top")
Ti_corr2 <- insert_yaxis_grob(Ti_p1, Ti_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Ti_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Ti_corr3 <- ggscatter(ACE_all, x = "Ti", y = "Ti_ICP", size = 1,
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
Ti_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Ti_XRF_violin_boxplot, Ti_ICP_violin_boxplot, Ti_corr2, Ti_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests - Ti

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Ti
y_all.reg <- ACE_all$Ti_ICP
x_all <- "Ti"
y_all <- "Ti_ICP"
x1.reg <- ACE_LM1$Ti
y1.reg <- ACE_LM1$Ti_ICP
x1 <- "Ti"
y1 <- "Ti_ICP"
x2.reg <- ACE_LM2$Ti
y2.reg <- ACE_LM2$Ti_ICP
x2 <- "Ti"
y2 <- "Ti_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Ti, Ti_sd, Ti_ICP, Ti_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Ti_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Ti_predict <- Ti_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ti_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Ti_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ti_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Ti_modf <- fortify(model_1)
Ti_res <- ggplot(Ti_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Ti_predict,  Ti_qq.x1, Ti_res, Ti_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Ti, ACE_all$Ti_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Ti, ACE_LM1$Ti_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Ti, ACE_LM2$Ti_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Ti, Ti_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Ti, ACE_all$Ti_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Ti, ACE_LM1$Ti_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Ti, ACE_LM2$Ti_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ti/Ti_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Mn ----------------------------------------------------------------------
element_title <- "Mn"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Mn Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Mn_corr1 <- ggscatter(ACE_LM1, x = "Mn", y = "Mn_ICP",
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
Mn_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Mn_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Mn, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Mn_XRF_violin_boxplot <- Mn_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Mn_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Mn_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Mn_ICP_violin_boxplot <- Mn_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Mn_pmain <- ggplot(ACE_all, aes(x = Mn, y = Mn_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Mn_xdens <- axis_canvas(Mn_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Mn, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Mn_ydens <- axis_canvas(Mn_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Mn_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Mn_p1 <- insert_xaxis_grob(Mn_pmain, Mn_xdens, grid::unit(.2, "null"), position = "top")
Mn_corr2 <- insert_yaxis_grob(Mn_p1, Mn_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Mn_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Mn_corr3 <- ggscatter(ACE_all, x = "Mn", y = "Mn_ICP", size = 1,
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
Mn_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Mn_XRF_violin_boxplot, Mn_ICP_violin_boxplot, Mn_corr2, Mn_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests - Mn

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Mn
y_all.reg <- ACE_all$Mn_ICP
x_all <- "Mn"
y_all <- "Mn_ICP"
x1.reg <- ACE_LM1$Mn
y1.reg <- ACE_LM1$Mn_ICP
x1 <- "Mn"
y1 <- "Mn_ICP"
x2.reg <- ACE_LM2$Mn
y2.reg <- ACE_LM2$Mn_ICP
x2 <- "Mn"
y2 <- "Mn_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Mn, Mn_sd, Mn_ICP, Mn_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Mn_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Mn_predict <- Mn_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Mn_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Mn_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Mn_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Mn_modf <- fortify(model_1)
Mn_res <- ggplot(Mn_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Mn_predict,  Mn_qq.x1, Mn_res, Mn_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Mn, ACE_all$Mn_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Mn, ACE_LM1$Mn_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Mn, ACE_LM2$Mn_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Mn, Mn_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Mn, ACE_all$Mn_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Mn, ACE_LM1$Mn_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Mn, ACE_LM2$Mn_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Mn/Mn_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Fe ----------------------------------------------------------------------
element_title <- "Fe"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Fe Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Fe_corr1 <- ggscatter(ACE_LM1, x = "Fe", y = "Fe_ICP",
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
Fe_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Fe_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Fe, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Fe_XRF_violin_boxplot <- Fe_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Fe_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Fe_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = , face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Fe_ICP_violin_boxplot <- Fe_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Fe_pmain <- ggplot(ACE_all, aes(x = Fe, y = Fe_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Fe_xdens <- axis_canvas(Fe_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Fe, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Fe_ydens <- axis_canvas(Fe_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Fe_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Fe_p1 <- insert_xaxis_grob(Fe_pmain, Fe_xdens, grid::unit(.2, "null"), position = "top")
Fe_corr2 <- insert_yaxis_grob(Fe_p1, Fe_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Fe_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Fe_corr3 <- ggscatter(ACE_all, x = "Fe", y = "Fe_ICP", size = 1,
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
Fe_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Fe_XRF_violin_boxplot, Fe_ICP_violin_boxplot, Fe_corr2, Fe_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests - Fe

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Fe
y_all.reg <- ACE_all$Fe_ICP
x_all <- "Fe"
y_all <- "Fe_ICP"

x1.reg <- ACE_LM1$Fe
y1.reg <- ACE_LM1$Fe_ICP
x1 <- "Fe"
y1 <- "Fe_ICP"

x2.reg <- ACE_LM2$Fe
y2.reg <- ACE_LM2$Fe_ICP
x2 <- "Fe"
y2 <- "Fe_ICP"

# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)

# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)

# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")

# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Fe, Fe_sd, Fe_ICP, Fe_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Fe_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Fe_predict <- Fe_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Fe_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Fe_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Fe_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Fe_modf <- fortify(model_1)
Fe_res <- ggplot(Fe_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Fe_predict,  Fe_qq.x1, Fe_res, Fe_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Fe, ACE_all$Fe_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Fe, ACE_LM1$Fe_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Fe, ACE_LM2$Fe_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Fe, Fe_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Fe, ACE_all$Fe_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Fe, ACE_LM1$Fe_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Fe, ACE_LM2$Fe_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Fe/Fe_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Co ----------------------------------------------------------------------
element_title <- "Co"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Co Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Co_corr1 <- ggscatter(ACE_LM1, x = "Co", y = "Co_ICP",
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
Co_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Co_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Co, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Co_XRF_violin_boxplot <- Co_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Co_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Co_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Co_ICP_violin_boxplot <- Co_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Co_pmain <- ggplot(ACE_all, aes(x = Co, y = Co_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Co_xdens <- axis_canvas(Co_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Co, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Co_ydens <- axis_canvas(Co_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Co_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Co_p1 <- insert_xaxis_grob(Co_pmain, Co_xdens, grid::unit(.2, "null"), position = "top")
Co_corr2 <- insert_yaxis_grob(Co_p1, Co_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Co_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Co_corr3 <- ggscatter(ACE_all, x = "Co", y = "Co_ICP", size = 1,
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
Co_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Co_XRF_violin_boxplot, Co_ICP_violin_boxplot, Co_corr2, Co_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Co
y_all.reg <- ACE_all$Co_ICP
x_all <- "Co"
y_all <- "Co_ICP"
x1.reg <- ACE_LM1$Co
y1.reg <- ACE_LM1$Co_ICP
x1 <- "Co"
y1 <- "Co_ICP"
x2.reg <- ACE_LM2$Co
y2.reg <- ACE_LM2$Co_ICP
x2 <- "Co"
y2 <- "Co_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Co, Co_sd, Co_ICP, Co_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Co_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Co_predict <- Co_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Co_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Co_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Co_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Co_modf <- fortify(model_1)
Co_res <- ggplot(Co_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Co_predict,  Co_qq.x1, Co_res, Co_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Co, ACE_all$Co_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Co, ACE_LM1$Co_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Co, ACE_LM2$Co_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Co, Co_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Co, ACE_all$Co_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Co, ACE_LM1$Co_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Co, ACE_LM2$Co_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Co/Co_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Ni ----------------------------------------------------------------------
element_title <- "Ni"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Ni Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Ni_corr1 <- ggscatter(ACE_LM1, x = "Ni", y = "Ni_ICP",
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
Ni_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Ni_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Ni, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Ni_XRF_violin_boxplot <- Ni_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Ni_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Ni_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Ni_ICP_violin_boxplot <- Ni_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Ni_pmain <- ggplot(ACE_all, aes(x = Ni, y = Ni_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Ni_xdens <- axis_canvas(Ni_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Ni, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Ni_ydens <- axis_canvas(Ni_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Ni_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Ni_p1 <- insert_xaxis_grob(Ni_pmain, Ni_xdens, grid::unit(.2, "null"), position = "top")
Ni_corr2 <- insert_yaxis_grob(Ni_p1, Ni_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Ni_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Ni_corr3 <- ggscatter(ACE_all, x = "Ni", y = "Ni_ICP", size = 1,
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
Ni_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Ni_XRF_violin_boxplot, Ni_ICP_violin_boxplot, Ni_corr2, Ni_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Ni
y_all.reg <- ACE_all$Ni_ICP
x_all <- "Ni"
y_all <- "Ni_ICP"
x1.reg <- ACE_LM1$Ni
y1.reg <- ACE_LM1$Ni_ICP
x1 <- "Ni"
y1 <- "Ni_ICP"
x2.reg <- ACE_LM2$Ni
y2.reg <- ACE_LM2$Ni_ICP
x2 <- "Ni"
y2 <- "Ni_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Ni, Ni_sd, Ni_ICP, Ni_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Ni_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Ni_predict <- Ni_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ni_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Ni_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ni_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Ni_modf <- fortify(model_1)
Ni_res <- ggplot(Ni_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Ni_predict,  Ni_qq.x1, Ni_res, Ni_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Ni, ACE_all$Ni_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Ni, ACE_LM1$Ni_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Ni, ACE_LM2$Ni_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Ni, Ni_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Ni, ACE_all$Ni_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Ni, ACE_LM1$Ni_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Ni, ACE_LM2$Ni_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Ni/Ni_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Cu ----------------------------------------------------------------------
element_title <- "Cu"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Cu Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Cu_corr1 <- ggscatter(ACE_LM1, x = "Cu", y = "Cu_ICP",
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
Cu_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Cu_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Cu, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Cu_XRF_violin_boxplot <- Cu_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Cu_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Cu_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Cu_ICP_violin_boxplot <- Cu_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Cu_pmain <- ggplot(ACE_all, aes(x = Cu, y = Cu_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Cu_xdens <- axis_canvas(Cu_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Cu, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Cu_ydens <- axis_canvas(Cu_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Cu_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Cu_p1 <- insert_xaxis_grob(Cu_pmain, Cu_xdens, grid::unit(.2, "null"), position = "top")
Cu_corr2 <- insert_yaxis_grob(Cu_p1, Cu_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Cu_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Cu_corr3 <- ggscatter(ACE_all, x = "Cu", y = "Cu_ICP", size = 1,
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
Cu_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Cu_XRF_violin_boxplot, Cu_ICP_violin_boxplot, Cu_corr2, Cu_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Cu
y_all.reg <- ACE_all$Cu_ICP
x_all <- "Cu"
y_all <- "Cu_ICP"
x1.reg <- ACE_LM1$Cu
y1.reg <- ACE_LM1$Cu_ICP
x1 <- "Cu"
y1 <- "Cu_ICP"
x2.reg <- ACE_LM2$Cu
y2.reg <- ACE_LM2$Cu_ICP
x2 <- "Cu"
y2 <- "Cu_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Cu, Cu_sd, Cu_ICP, Cu_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Cu_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Cu_predict <- Cu_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Cu_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Cu_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Cu_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Cu_modf <- fortify(model_1)
Cu_res <- ggplot(Cu_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Cu_predict,  Cu_qq.x1, Cu_res, Cu_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Cu, ACE_all$Cu_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Cu, ACE_LM1$Cu_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Cu, ACE_LM2$Cu_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Cu, Cu_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Cu, ACE_all$Cu_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Cu, ACE_LM1$Cu_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Cu, ACE_LM2$Cu_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test 
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Cu/Cu_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Zn ----------------------------------------------------------------------
element_title <- "Zn"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Zn Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Zn_corr1 <- ggscatter(ACE_LM1, x = "Zn", y = "Zn_ICP",
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
Zn_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Zn_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Zn, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Zn_XRF_violin_boxplot <- Zn_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Zn_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Zn_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Zn_ICP_violin_boxplot <- Zn_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Zn_pmain <- ggplot(ACE_all, aes(x = Zn, y = Zn_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Zn_xdens <- axis_canvas(Zn_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Zn, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Zn_ydens <- axis_canvas(Zn_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Zn_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Zn_p1 <- insert_xaxis_grob(Zn_pmain, Zn_xdens, grid::unit(.2, "null"), position = "top")
Zn_corr2 <- insert_yaxis_grob(Zn_p1, Zn_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Zn_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Zn_corr3 <- ggscatter(ACE_all, x = "Zn", y = "Zn_ICP", size = 1,
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
Zn_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Zn_XRF_violin_boxplot, Zn_ICP_violin_boxplot, Zn_corr2, Zn_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Zn
y_all.reg <- ACE_all$Zn_ICP
x_all <- "Zn"
y_all <- "Zn_ICP"
x1.reg <- ACE_LM1$Zn
y1.reg <- ACE_LM1$Zn_ICP
x1 <- "Zn"
y1 <- "Zn_ICP"
x2.reg <- ACE_LM2$Zn
y2.reg <- ACE_LM2$Zn_ICP
x2 <- "Zn"
y2 <- "Zn_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Zn, Zn_sd, Zn_ICP, Zn_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Zn_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Zn_predict <- Zn_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Zn_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Zn_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Zn_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Zn_modf <- fortify(model_1)
Zn_res <- ggplot(Zn_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Zn_predict,  Zn_qq.x1, Zn_res, Zn_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Zn, ACE_all$Zn_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Zn, ACE_LM1$Zn_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Zn, ACE_LM2$Zn_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Zn, Zn_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Zn, ACE_all$Zn_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Zn, ACE_LM1$Zn_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Zn, ACE_LM2$Zn_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zn/Zn_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)
# Rb ----------------------------------------------------------------------
element_title <- "Rb"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Rb Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Rb_corr1 <- ggscatter(ACE_LM1, x = "Rb", y = "Rb_ICP",
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
Rb_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Rb_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Rb, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Rb_XRF_violin_boxplot <- Rb_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Rb_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Rb_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Rb_ICP_violin_boxplot <- Rb_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Rb_pmain <- ggplot(ACE_all, aes(x = Rb, y = Rb_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Rb_xdens <- axis_canvas(Rb_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Rb, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Rb_ydens <- axis_canvas(Rb_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Rb_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Rb_p1 <- insert_xaxis_grob(Rb_pmain, Rb_xdens, grid::unit(.2, "null"), position = "top")
Rb_corr2 <- insert_yaxis_grob(Rb_p1, Rb_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Rb_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Rb_corr3 <- ggscatter(ACE_all, x = "Rb", y = "Rb_ICP", size = 1,
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
Rb_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Rb_XRF_violin_boxplot, Rb_ICP_violin_boxplot, Rb_corr2, Rb_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Rb
y_all.reg <- ACE_all$Rb_ICP
x_all <- "Rb"
y_all <- "Rb_ICP"
x1.reg <- ACE_LM1$Rb
y1.reg <- ACE_LM1$Rb_ICP
x1 <- "Rb"
y1 <- "Rb_ICP"
x2.reg <- ACE_LM2$Rb
y2.reg <- ACE_LM2$Rb_ICP
x2 <- "Rb"
y2 <- "Rb_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Rb, Rb_sd, Rb_ICP, Rb_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Rb_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Rb_predict <- Rb_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Rb_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Rb_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Rb_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Rb_modf <- fortify(model_1)
Rb_res <- ggplot(Rb_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Rb_predict,  Rb_qq.x1, Rb_res, Rb_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Rb, ACE_all$Rb_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Rb, ACE_LM1$Rb_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Rb, ACE_LM2$Rb_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Rb, Rb_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Rb, ACE_all$Rb_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Rb, ACE_LM1$Rb_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Rb, ACE_LM2$Rb_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Rb/Rb_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Sr ----------------------------------------------------------------------
element_title <- "Sr"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Sr Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Sr_corr1 <- ggscatter(ACE_LM1, x = "Sr", y = "Sr_ICP",
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
Sr_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Sr_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Sr, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Sr_XRF_violin_boxplot <- Sr_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Sr_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Sr_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Sr_ICP_violin_boxplot <- Sr_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Sr_pmain <- ggplot(ACE_all, aes(x = Sr, y = Sr_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Sr_xdens <- axis_canvas(Sr_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Sr, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Sr_ydens <- axis_canvas(Sr_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Sr_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Sr_p1 <- insert_xaxis_grob(Sr_pmain, Sr_xdens, grid::unit(.2, "null"), position = "top")
Sr_corr2 <- insert_yaxis_grob(Sr_p1, Sr_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Sr_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Sr_corr3 <- ggscatter(ACE_all, x = "Sr", y = "Sr_ICP", size = 1,
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
Sr_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Sr_XRF_violin_boxplot, Sr_ICP_violin_boxplot, Sr_corr2, Sr_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Sr
y_all.reg <- ACE_all$Sr_ICP
x_all <- "Sr"
y_all <- "Sr_ICP"
x1.reg <- ACE_LM1$Sr
y1.reg <- ACE_LM1$Sr_ICP
x1 <- "Sr"
y1 <- "Sr_ICP"
x2.reg <- ACE_LM2$Sr
y2.reg <- ACE_LM2$Sr_ICP
x2 <- "Sr"
y2 <- "Sr_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Sr, Sr_sd, Sr_ICP, Sr_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Sr_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Sr_predict <- Sr_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Sr_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Sr_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Sr_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Sr_modf <- fortify(model_1)
Sr_res <- ggplot(Sr_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Sr_predict,  Sr_qq.x1, Sr_res, Sr_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Sr, ACE_all$Sr_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Sr, ACE_LM1$Sr_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Sr, ACE_LM2$Sr_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Sr, Sr_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Sr, ACE_all$Sr_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Sr, ACE_LM1$Sr_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Sr, ACE_LM2$Sr_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Sr/Sr_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# Zr ----------------------------------------------------------------------
element_title <- "Zr"
XRF_title <- "XRF-CS / inc"
ICP_title <- "ICPMS (ppm)"
correlation_title <- "Zr Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
Zr_corr1 <- ggscatter(ACE_LM1, x = "Zr", y = "Zr_ICP",
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
Zr_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
Zr_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=Zr, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Zr_XRF_violin_boxplot <- Zr_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
Zr_ICP_violin <- ggplot(ACE_all, aes(x=Site, y=Zr_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Zr_ICP_violin_boxplot <- Zr_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
Zr_pmain <- ggplot(ACE_all, aes(x = Zr, y = Zr_ICP, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
Zr_xdens <- axis_canvas(Zr_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = Zr, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Zr_ydens <- axis_canvas(Zr_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = Zr_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Zr_p1 <- insert_xaxis_grob(Zr_pmain, Zr_xdens, grid::unit(.2, "null"), position = "top")
Zr_corr2 <- insert_yaxis_grob(Zr_p1, Zr_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Zr_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
Zr_corr3 <- ggscatter(ACE_all, x = "Zr", y = "Zr_ICP", size = 1,
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
Zr_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Zr_XRF_violin_boxplot, Zr_ICP_violin_boxplot, Zr_corr2, Zr_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$Zr
y_all.reg <- ACE_all$Zr_ICP
x_all <- "Zr"
y_all <- "Zr_ICP"
x1.reg <- ACE_LM1$Zr
y1.reg <- ACE_LM1$Zr_ICP
x1 <- "Zr"
y1 <- "Zr_ICP"
x2.reg <- ACE_LM2$Zr
y2.reg <- ACE_LM2$Zr_ICP
x2 <- "Zr"
y2 <- "Zr_ICP"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, Zr, Zr_sd, Zr_ICP, Zr_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
Zr_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
Zr_predict <- Zr_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Rb_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
Zr_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Zr_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
Zr_modf <- fortify(model_1)
Zr_res <- ggplot(Zr_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(Zr_predict,  Zr_qq.x1, Zr_res, Zr_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$Zr, ACE_all$Zr_ICP)
ctest_LM1 <- cor.test(ACE_LM1$Zr, ACE_LM1$Zr_ICP)
ctest_LM2 <- cor.test(ACE_LM2$Zr, ACE_LM2$Zr_ICP)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(Zr, Zr_ICP))
# Examine co-variance in all models
cov_all <- cov(ACE_all$Zr, ACE_all$Zr_ICP)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$Zr, ACE_LM1$Zr_ICP)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$Zr, ACE_LM2$Zr_ICP)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/Zr/Zr_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# DM ----------------------------------------------------------------------
element_title <- "DM"
XRF_title <- "coh_inc XRF-CS"
ICP_title <- "Dry mass (%)"
correlation_title <- "DM Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Fig1 - Correlation plot
theme_set(theme_classic(10))
DM_corr1 <- ggscatter(ACE_LM1, x = "coh_inc", y = "dry_mass_pc",
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
  labs(x = paste(XRF_title), 
       y = paste(ICP_title)) +
  #scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("ACE (OLS): ", element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=10, face="bold"))
DM_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Fig1_Correlation_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Figure 2 - Violin + box plots, scatterplot + density, Individual Site Correlations
# a) XRF-CS data
theme_set(theme_classic(10))
DM_XRF_violin <- ggplot(ACE_all, aes(x=Site, y=coh_inc, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
DM_XRF_violin_boxplot <- DM_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
# b) ICPMS data
theme_set(theme_classic(10))
dry_mass_pc_violin <- ggplot(ACE_all, aes(x=Site, y=dry_mass_pc, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
dry_mass_pc_violin_boxplot <- dry_mass_pc_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
# c) Scatterplot and density plot
theme_set(theme_classic(10))
DM_pmain <- ggplot(ACE_all, aes(x = coh_inc, y = dry_mass_pc, color = Site)) +
  geom_point() + 
  ggpubr::color_palette(palette_set) +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "none") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
# Marginal densities along x axis
DM_xdens <- axis_canvas(DM_pmain, axis = "x") +
  geom_density(data = ACE_all, aes(x = coh_inc, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
DM_ydens <- axis_canvas(DM_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_all, aes(x = dry_mass_pc, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
DM_p1 <- insert_xaxis_grob(DM_pmain, DM_xdens, grid::unit(.2, "null"), position = "top")
DM_corr2 <- insert_yaxis_grob(DM_p1, DM_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(DM_corr2)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Fig2c_Scatterplot_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# d) Correlation per site
theme_set(theme_classic(10))
DM_corr3 <- ggscatter(ACE_all, x = "coh_inc", y = "dry_mass_pc", size = 1,
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
  labs(x = paste(XRF_title), 
       y = paste(ICP_title)) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  ggtitle(paste("ACE: ", element_title, ICP_title ,"vs", XRF_title))
DM_corr3
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(DM_XRF_violin_boxplot, dry_mass_pc_violin_boxplot, DM_corr2, DM_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Fig2_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Set up x and y variables - Models all, LM1, LM2
x_all.reg <- ACE_all$coh_inc
y_all.reg <- ACE_all$dry_mass_pc
x_all <- "coh_inc"
y_all <- "dry_mass_pc"
x1.reg <- ACE_LM1$coh_inc
y1.reg <- ACE_LM1$dry_mass_pc
x1 <- "coh_inc"
y1 <- "dry_mass_pc"
x2.reg <- ACE_LM2$coh_inc
y2.reg <- ACE_LM2$dry_mass_pc
x2 <- "coh_inc"
y2 <- "dry_mass_pc"
# Build linear regression models
model_all <- lm(y_all.reg~x_all.reg, data = ACE_all)
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
model_2 <- lm(y2.reg~x2.reg, data = ACE_LM2)
# Choose model to use & summary statistics
model_1
summary(model_1)
confint(model_1)
# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:SH20_age, coh_inc, coh_inc_sd, dry_mass_pc, dry_mass_err, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Model_predict.csv", row.names = FALSE)

# Figure 3 - qq plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Figure 3a - OLS Linear model 1
DM_p1 <- ggplot(data_1, aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Site), size = 2, alpha = 1) +
  ggpubr::color_palette(palette_set) +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, 
           vjust = 2, hjust = -0.2, size = 3, color = "black") + #label.sep = "\n"
  stat_regline_equation(label.x = -Inf, label.y = Inf,
                        vjust = 4, hjust = -0.3, size = 3, color = "black") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, bycol = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "lightblue") +
  labs(x = paste(element_title, XRF_title), 
       y = paste(element_title, ICP_title)) 
# Add prediction intervals
DM_predict <- DM_p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle(paste("ACE 95% CI & pred: ", element_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
DM_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Fig3a_Correlation_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
#Figure 3c, d
theme_set(theme_classic(base_size=10))
DM_qq.x1 <- ggqqplot(data_1, x = x1) +
  ggtitle(paste("q-q plot: ", element_title, XRF_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
DM_qq.y1 <- ggqqplot(data_1, x = y1) +
  ggtitle(paste("q-q plot: ", element_title, ICP_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
# Figure 3b - Residuals vs fitted values plot
theme_set(theme_classic(base_size=10))
DM_modf <- fortify(model_1)
DM_res <- ggplot(DM_modf, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle(paste("Residual: ", element_title)) +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))
ggarrange(DM_predict,  DM_qq.x1, DM_res, DM_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Fig3_Correlation_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics summary

# Hypothesis test 
ctest_all <- cor.test(ACE_all$coh_inc, ACE_all$dry_mass_pc)
ctest_LM1 <- cor.test(ACE_LM1$coh_inc, ACE_LM1$dry_mass_pc)
ctest_LM2 <- cor.test(ACE_LM2$coh_inc, ACE_LM2$dry_mass_pc)
# Correlation stats by site
correlate_ALL <- ACE_all %>%
  group_by(Site) %>% 
  summarise(r = cor(coh_inc, dry_mass_pc))
# Examine co-variance in all models
cov_all <- cov(ACE_all$coh_inc, ACE_all$dry_mass_pc)
round(cov_all,2)
cov_LM1 <- cov(ACE_LM1$coh_inc, ACE_LM1$dry_mass_pc)
round(cov_LM1,2)
cov_LM2 <- cov(ACE_LM2$coh_inc, ACE_LM2$dry_mass_pc)
round(cov_LM2,2)
# Shapiro-Wilk test all models
shapiro.test(x_all.reg)
shapiro.test(y_all.reg)
shapiro.test(x1.reg)
shapiro.test(y1.reg)
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# Durbin-Watson Test
dwtest(model_all)
dwtest(model_1)
dwtest(model_2)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/DM/DM_Summary_stats.txt")
print(correlate_ALL)
print(ctest_all)
print(ctest_LM1)
print(ctest_LM2)
print(summary(model_all))
print(confint(model_all))
print(summary(model_1))
print(confint(model_1))
print(summary(model_2))
print(confint(model_2))
print(dwtest(model_all))
print(dwtest(model_1))
print(dwtest(model_2))
print(shapiro.test(x_all.reg))
print(shapiro.test(y_all.reg))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
print(shapiro.test(x2.reg))
print(shapiro.test(y2.reg))
sink(file = NULL)

# -------------------------------------------------------------------------
# ACE Correlation summary
ggarrange(K_corr1, Ca_corr1, Ti_corr1, Mn_corr1, Fe_corr1, Co_corr1,
          Ni_corr1, Cu_corr1, Zn_corr1, Rb_corr1, Sr_corr1, Zr_corr1,
          nrow = 3, ncol = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_Correlation_Summary.pdf", 
       height = c(22.5), width = c(30), dpi = 600, units = "cm")
# ACE Scatterplot
plot_grid(K_corr2, Ca_corr2, Ti_corr2, Mn_corr2, Fe_corr2, Co_corr2,
          Ni_corr2, Cu_corr2, Zn_corr2, Rb_corr2, Sr_corr2, DM_corr2,
          nrow = 3, ncol = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_Correlation_inc/ACE_Scatterplot_Summary.pdf", 
       height = c(45), width = c(60), dpi = 600, units = "cm")
# -------------------------------------------------------------------------

# Linear modelling -----------------------------------------------------------

# Load required packages
library(performance) # linear model testing and graphical outputs
library(ggpmisc) # this package is used for simple/quick for labelling (stat_poly_eq)
library(lmtest) # linear model testing 

# Linear models & performance assessment /stats -----------------------------------------------

# Set up labels
#ACE_dataset <- ACE_all
ACE_dataset <- ACE_no_POB4
#ACE_dataset <-  ACE_no_POB4_PB1
#ACE_dataset <- ACE_no_PB1
site_title <- "ACE"
#itrax_dataset <- " (cps)"
itrax_dataset <- " / inc"
#itrax_dataset <- " [Ln(E/inc)]"
#itrax_dataset <- " (element/inc sum)"
#itrax_dataset <- " (clr)"

# Convert Site as a grouping variable
ACE_dataset$Site <- as.factor(ACE_dataset$Site)

# ACE_K -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_K_lm <- lm(K_ICP ~ K, data = ACE_dataset)
summary(ACE_K_lm)
glance(ACE_K_lm)
model_performance(ACE_K_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_OLS_summary.txt")
summary(ACE_K_lm)
glance(ACE_K_lm)
model_performance(ACE_K_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_lm) # Performance package summary check for heteroscedasticity
icc(ACE_K_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_K_wlm <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/(K_sd)^2)
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_OLS_wt_summary.txt")
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_lm) # Performance package summary check for heteroscedasticity
icc(ACE_K_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_K_model <- lm(K_ICP ~ K, data = ACE_dataset) # define model
ACE_K_wt <- 1 / lm(abs(ACE_K_model$residuals) ~ ACE_K_model$fitted.values)$fitted.values^2 #define weights to use
ACE_K_wls <- lm(K_ICP ~ K, data = ACE_dataset, weights=ACE_K_wt) #perform weighted least squares regression
# Checks
summary(ACE_K_wls) # summary stats
glance(ACE_K_wls) # summary stats including AIC
model_performance(ACE_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_WLS_summary.txt")
summary(ACE_K_wls)
glance(ACE_K_wls)
model_performance(ACE_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wls) # Performance package summary check for heteroscedasticity
icc(ACE_K_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_K_model_wt <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/K_sd^2) # define model
ACE_K_wt_wt <- 1 / lm(abs(ACE_K_model_wt$residuals) ~ ACE_K_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_K_wls_wt <- lm(K_ICP ~ K, data = ACE_dataset, weights=ACE_K_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_K_wls_wt) # summary stats
glance(ACE_K_wls_wt) # summary stats including AIC
model_performance(ACE_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_WLS_wt_summary.txt")
summary(ACE_K_wls_wt)
glance(ACE_K_wls_wt)
model_performance(ACE_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_K_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "K"
theme_set(theme_classic(10))
ACE_K <- ggplot(ACE_dataset, aes(x = K, y = K_ICP)) +
  geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/K_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_K_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_K_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_K
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_K_lm_p <- function(ACE_K_lm) {
  f <- summary(ACE_K_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_lm_p(ACE_K_lm)

ACE_K_lm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_lm_p(ACE_K_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_K_wlm_p <- function(ACE_K_wlm) {
  f <- summary(ACE_K_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wlm_p(ACE_K_wlm)

ACE_K_wlm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/(K_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wlm_p(ACE_K_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_K_wls_p <- function(ACE_K_wls) {
  f <- summary(ACE_K_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wls_p(ACE_K_wls)

ACE_K_wls_eqn <- function(df){
  m <- ACE_K_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wls_p(ACE_K_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_K_wls_wt_p <- function(ACE_K_wls_wt) {
  f <- summary(ACE_K_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wls_wt_p(ACE_K_wls_wt)

ACE_K_wls_wt_eqn <- function(df){
  m <- ACE_K_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wls_wt_p(ACE_K_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_K_final <- ACE_K + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_K_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_K_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_K_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_K_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_K_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/K/K_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# ACE_Ca -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ca_lm <- lm(Ca_ICP ~ Ca, data = ACE_dataset)
summary(ACE_Ca_lm)
glance(ACE_Ca_lm)
model_performance(ACE_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_OLS_summary.txt")
summary(ACE_Ca_lm)
glance(ACE_Ca_lm)
model_performance(ACE_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Ca_wlm <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/(Ca_sd)^2)
summary(ACE_Ca_wlm)
glance(ACE_Ca_wlm)
model_performance(ACE_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_OLS_wt_summary.txt")
summary(ACE_Ca_wlm)
glance(ACE_Ca_wlm)
model_performance(ACE_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Ca_model <- lm(Ca_ICP ~ Ca, data = ACE_dataset) # define model
ACE_Ca_wt <- 1 / lm(abs(ACE_Ca_model$residuals) ~ ACE_Ca_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ca_wls <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ca_wls) # summary stats
glance(ACE_Ca_wls) # summary stats including AIC
model_performance(ACE_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_WLS_summary.txt")
summary(ACE_Ca_wls)
glance(ACE_Ca_wls)
model_performance(ACE_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Ca_model_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/Ca_sd^2) # define model
ACE_Ca_wt_wt <- 1 / lm(abs(ACE_Ca_model_wt$residuals) ~ ACE_Ca_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ca_wls_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ca_wls_wt) # summary stats
glance(ACE_Ca_wls_wt) # summary stats including AIC
model_performance(ACE_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_WLS_wt_summary.txt")
summary(ACE_Ca_wls_wt)
glance(ACE_Ca_wls_wt)
model_performance(ACE_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Ca"
theme_set(theme_classic(10))
ACE_Ca <- ggplot(ACE_dataset, aes(x = Ca, y = Ca_ICP)) +
  geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Ca_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Ca_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Ca_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Ca
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ca_lm_p <- function(ACE_Ca_lm) {
  f <- summary(ACE_Ca_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_lm_p(ACE_Ca_lm)

ACE_Ca_lm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_lm_p(ACE_Ca_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ca_wlm_p <- function(ACE_Ca_wlm) {
  f <- summary(ACE_Ca_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wlm_p(ACE_Ca_wlm)

ACE_Ca_wlm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/(Ca_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wlm_p(ACE_Ca_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ca_wls_p <- function(ACE_Ca_wls) {
  f <- summary(ACE_Ca_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wls_p(ACE_Ca_wls)

ACE_Ca_wls_eqn <- function(df){
  m <- ACE_Ca_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wls_p(ACE_Ca_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ca_wls_wt_p <- function(ACE_Ca_wls_wt) {
  f <- summary(ACE_Ca_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wls_wt_p(ACE_Ca_wls_wt)

ACE_Ca_wls_wt_eqn <- function(df){
  m <- ACE_Ca_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wls_wt_p(ACE_Ca_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Ca_final <- ACE_Ca + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ca_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ca_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ca_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ca_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ca_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ca/Ca_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# ACE_Ti -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ti_lm <- lm(Ti_ICP ~ Ti, data = ACE_dataset)
summary(ACE_Ti_lm)
glance(ACE_Ti_lm)
model_performance(ACE_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_OLS_summary.txt")
summary(ACE_Ti_lm)
glance(ACE_Ti_lm)
model_performance(ACE_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Ti_wlm <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/(Ti_sd)^2)
summary(ACE_Ti_wlm)
glance(ACE_Ti_wlm)
model_performance(ACE_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_OLS_wt_summary.txt")
summary(ACE_Ti_wlm)
glance(ACE_Ti_wlm)
model_performance(ACE_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Ti_model <- lm(Ti_ICP ~ Ti, data = ACE_dataset) # define model
ACE_Ti_wt <- 1 / lm(abs(ACE_Ti_model$residuals) ~ ACE_Ti_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ti_wls <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ti_wls) # summary stats
glance(ACE_Ti_wls) # summary stats including AIC
model_performance(ACE_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_WLS_summary.txt")
summary(ACE_Ti_wls)
glance(ACE_Ti_wls)
model_performance(ACE_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Ti_model_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/Ti_sd^2) # define model
ACE_Ti_wt_wt <- 1 / lm(abs(ACE_Ti_model_wt$residuals) ~ ACE_Ti_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ti_wls_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ti_wls_wt) # summary stats
glance(ACE_Ti_wls_wt) # summary stats including AIC
model_performance(ACE_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_WLS_wt_summary.txt")
summary(ACE_Ti_wls_wt)
glance(ACE_Ti_wls_wt)
model_performance(ACE_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear models
element_title <- "Ti"
theme_set(theme_classic(10))
ACE_Ti <- ggplot(ACE_dataset, aes(x = Ti, y = Ti_ICP)) +
  geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Ti_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Ti_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Ti_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Ti
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ti_lm_p <- function(ACE_Ti_lm) {
  f <- summary(ACE_Ti_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_lm_p(ACE_Ti_lm)

ACE_Ti_lm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_lm_p(ACE_Ti_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ti_wlm_p <- function(ACE_Ti_wlm) {
  f <- summary(ACE_Ti_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wlm_p(ACE_Ti_wlm)

ACE_Ti_wlm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/(Ti_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wlm_p(ACE_Ti_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ti_wls_p <- function(ACE_Ti_wls) {
  f <- summary(ACE_Ti_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wls_p(ACE_Ti_wls)

ACE_Ti_wls_eqn <- function(df){
  m <- ACE_Ti_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wls_p(ACE_Ti_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ti_wls_wt_p <- function(ACE_Ti_wls_wt) {
  f <- summary(ACE_Ti_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wls_wt_p(ACE_Ti_wls_wt)

ACE_Ti_wls_wt_eqn <- function(df){
  m <- ACE_Ti_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wls_wt_p(ACE_Ti_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Ti_final <- ACE_Ti + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ti_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ti_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ti_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ti_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ti_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ti/Ti_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ACE_Mn -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Mn_lm <- lm(Mn_ICP ~ Mn, data = ACE_dataset)
summary(ACE_Mn_lm)
glance(ACE_Mn_lm)
model_performance(ACE_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_OLS_summary.txt")
summary(ACE_Mn_lm)
glance(ACE_Mn_lm)
model_performance(ACE_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Mn_wlm <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/(Mn_sd)^2)
summary(ACE_Mn_wlm)
glance(ACE_Mn_wlm)
model_performance(ACE_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_OLS_wt_summary.txt")
summary(ACE_Mn_wlm)
glance(ACE_Mn_wlm)
model_performance(ACE_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Mn_model <- lm(Mn_ICP ~ Mn, data = ACE_dataset) # define model
ACE_Mn_wt <- 1 / lm(abs(ACE_Mn_model$residuals) ~ ACE_Mn_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Mn_wls <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt) #perform weighted least squares regression
# Checks
summary(ACE_Mn_wls) # summary stats
glance(ACE_Mn_wls) # summary stats including AIC
model_performance(ACE_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_WLS_summary.txt")
summary(ACE_Mn_wls)
glance(ACE_Mn_wls)
model_performance(ACE_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Mn_model_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/Mn_sd^2) # define model
ACE_Mn_wt_wt <- 1 / lm(abs(ACE_Mn_model_wt$residuals) ~ ACE_Mn_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Mn_wls_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Mn_wls_wt) # summary stats
glance(ACE_Mn_wls_wt) # summary stats including AIC
model_performance(ACE_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_WLS_wt_summary.txt")
summary(ACE_Mn_wls_wt)
glance(ACE_Mn_wls_wt)
model_performance(ACE_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Mn"
theme_set(theme_classic(10))
ACE_Mn <- ggplot(ACE_dataset, aes(x = Mn, y = Mn_ICP)) +
  geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Mn_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Mn_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Mn_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Mn
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Mn_lm_p <- function(ACE_Mn_lm) {
  f <- summary(ACE_Mn_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_lm_p(ACE_Mn_lm)

ACE_Mn_lm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_lm_p(ACE_Mn_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Mn_wlm_p <- function(ACE_Mn_wlm) {
  f <- summary(ACE_Mn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wlm_p(ACE_Mn_wlm)

ACE_Mn_wlm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/(Mn_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wlm_p(ACE_Mn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Mn_wls_p <- function(ACE_Mn_wls) {
  f <- summary(ACE_Mn_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wls_p(ACE_Mn_wls)

ACE_Mn_wls_eqn <- function(df){
  m <- ACE_Mn_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wls_p(ACE_Mn_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Mn_wls_wt_p <- function(ACE_Mn_wls_wt) {
  f <- summary(ACE_Mn_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wls_wt_p(ACE_Mn_wls_wt)

ACE_Mn_wls_wt_eqn <- function(df){
  m <- ACE_Mn_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wls_wt_p(ACE_Mn_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Mn_final <- ACE_Mn + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Mn_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Mn_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Mn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Mn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Mn_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Mn/Mn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ACE_Fe -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Fe_lm <- lm(Fe_ICP ~ Fe, data = ACE_dataset)
summary(ACE_Fe_lm)
glance(ACE_Fe_lm)
model_performance(ACE_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_OLS_summary.txt")
summary(ACE_Fe_lm)
glance(ACE_Fe_lm)
model_performance(ACE_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Fe_wlm <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/(Fe_sd)^2)
summary(ACE_Fe_wlm)
glance(ACE_Fe_wlm)
model_performance(ACE_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_OLS_wt_summary.txt")
summary(ACE_Fe_wlm)
glance(ACE_Fe_wlm)
model_performance(ACE_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Fe_model <- lm(Fe_ICP ~ Fe, data = ACE_dataset) # define model
ACE_Fe_wt <- 1 / lm(abs(ACE_Fe_model$residuals) ~ ACE_Fe_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Fe_wls <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt) #perform weighted least squares regression
# Checks
summary(ACE_Fe_wls) # summary stats
glance(ACE_Fe_wls) # summary stats including AIC
model_performance(ACE_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_WLS_summary.txt")
summary(ACE_Fe_wls)
glance(ACE_Fe_wls)
model_performance(ACE_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Fe_model_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/Fe_sd^2) # define model
ACE_Fe_wt_wt <- 1 / lm(abs(ACE_Fe_model_wt$residuals) ~ ACE_Fe_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Fe_wls_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Fe_wls_wt) # summary stats
glance(ACE_Fe_wls_wt) # summary stats including AIC
model_performance(ACE_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_WLS_wt_summary.txt")
summary(ACE_Fe_wls_wt)
glance(ACE_Fe_wls_wt)
model_performance(ACE_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Fe"
theme_set(theme_classic(10))
ACE_Fe <- ggplot(ACE_dataset, aes(x = Fe, y = Fe_ICP)) +
  geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Fe_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Fe_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Fe_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Fe
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Fe_lm_p <- function(ACE_Fe_lm) {
  f <- summary(ACE_Fe_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_lm_p(ACE_Fe_lm)

ACE_Fe_lm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_lm_p(ACE_Fe_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Fe_wlm_p <- function(ACE_Fe_wlm) {
  f <- summary(ACE_Fe_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wlm_p(ACE_Fe_wlm)

ACE_Fe_wlm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/(Fe_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wlm_p(ACE_Fe_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Fe_wls_p <- function(ACE_Fe_wls) {
  f <- summary(ACE_Fe_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wls_p(ACE_Fe_wls)

ACE_Fe_wls_eqn <- function(df){
  m <- ACE_Fe_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wls_p(ACE_Fe_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Fe_wls_wt_p <- function(ACE_Fe_wls_wt) {
  f <- summary(ACE_Fe_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wls_wt_p(ACE_Fe_wls_wt)

ACE_Fe_wls_wt_eqn <- function(df){
  m <- ACE_Fe_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wls_wt_p(ACE_Fe_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Fe_final <- ACE_Fe + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Fe_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Fe_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Fe_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Fe_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Fe_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Fe/Fe_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ACE_Co -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Co_lm <- lm(Co_ICP ~ Co, data = ACE_dataset)
summary(ACE_Co_lm)
glance(ACE_Co_lm)
model_performance(ACE_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_OLS_summary.txt")
summary(ACE_Co_lm)
glance(ACE_Co_lm)
model_performance(ACE_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Co_wlm <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/(Co_sd)^2)
summary(ACE_Co_wlm)
glance(ACE_Co_wlm)
model_performance(ACE_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_OLS_wt_summary.txt")
summary(ACE_Co_wlm)
glance(ACE_Co_wlm)
model_performance(ACE_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Co_model <- lm(Co_ICP ~ Co, data = ACE_dataset) # define model
ACE_Co_wt <- 1 / lm(abs(ACE_Co_model$residuals) ~ ACE_Co_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Co_wls <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=ACE_Co_wt) #perform weighted least squares regression
# Checks
summary(ACE_Co_wls) # summary stats
glance(ACE_Co_wls) # summary stats including AIC
model_performance(ACE_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_WLS_summary.txt")
summary(ACE_Co_wls)
glance(ACE_Co_wls)
model_performance(ACE_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Co_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Co_model_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/Co_sd^2) # define model
ACE_Co_wt_wt <- 1 / lm(abs(ACE_Co_model_wt$residuals) ~ ACE_Co_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Co_wls_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=ACE_Co_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Co_wls_wt) # summary stats
glance(ACE_Co_wls_wt) # summary stats including AIC
model_performance(ACE_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_WLS_wt_summary.txt")
summary(ACE_Co_wls_wt)
glance(ACE_Co_wls_wt)
model_performance(ACE_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Co_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Co"
theme_set(theme_classic(10))
ACE_Co <- ggplot(ACE_dataset, aes(x = Co, y = Co_ICP)) +
  geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Co_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Co_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Co_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Co
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Co_lm_p <- function(ACE_Co_lm) {
  f <- summary(ACE_Co_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_lm_p(ACE_Co_lm)

ACE_Co_lm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_lm_p(ACE_Co_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Co_wlm_p <- function(ACE_Co_wlm) {
  f <- summary(ACE_Co_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wlm_p(ACE_Co_wlm)

ACE_Co_wlm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/(Co_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wlm_p(ACE_Co_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Co_wls_p <- function(ACE_Co_wls) {
  f <- summary(ACE_Co_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wls_p(ACE_Co_wls)

ACE_Co_wls_eqn <- function(df){
  m <- ACE_Co_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wls_p(ACE_Co_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Co_wls_wt_p <- function(ACE_Co_wls_wt) {
  f <- summary(ACE_Co_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wls_wt_p(ACE_Co_wls_wt)

ACE_Co_wls_wt_eqn <- function(df){
  m <- ACE_Co_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wls_wt_p(ACE_Co_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Co_final <- ACE_Co + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Co_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Co_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Co_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Co_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Co_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Co/Co_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ACE_Ni -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ni_lm <- lm(Ni_ICP ~ Ni, data = ACE_dataset)
summary(ACE_Ni_lm)
glance(ACE_Ni_lm)
model_performance(ACE_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_OLS_summary.txt")
summary(ACE_Ni_lm)
glance(ACE_Ni_lm)
model_performance(ACE_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Ni_wlm <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/(Ni_sd)^2)
summary(ACE_Ni_wlm)
glance(ACE_Ni_wlm)
model_performance(ACE_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_OLS_wt_summary.txt")
summary(ACE_Ni_wlm)
glance(ACE_Ni_wlm)
model_performance(ACE_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Ni_model <- lm(Ni_ICP ~ Ni, data = ACE_dataset) # define model
ACE_Ni_wt <- 1 / lm(abs(ACE_Ni_model$residuals) ~ ACE_Ni_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ni_wls <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ni_wls) # summary stats
glance(ACE_Ni_wls) # summary stats including AIC
model_performance(ACE_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_WLS_summary.txt")
summary(ACE_Ni_wls)
glance(ACE_Ni_wls)
model_performance(ACE_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_wls) # check for random effects - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Ni_model_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/Ni_sd^2) # define model
ACE_Ni_wt_wt <- 1 / lm(abs(ACE_Ni_model_wt$residuals) ~ ACE_Ni_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ni_wls_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ni_wls_wt) # summary stats
glance(ACE_Ni_wls_wt) # summary stats including AIC
model_performance(ACE_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_WLS_wt_summary.txt")
summary(ACE_Ni_wls_wt)
glance(ACE_Ni_wls_wt)
model_performance(ACE_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Ni"
theme_set(theme_classic(10))
ACE_Ni <- ggplot(ACE_dataset, aes(x = Ni, y = Ni_ICP)) +
  geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Ni_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Ni_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Ni_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Ni
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ni_lm_p <- function(ACE_Ni_lm) {
  f <- summary(ACE_Ni_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_lm_p(ACE_Ni_lm)

ACE_Ni_lm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_lm_p(ACE_Ni_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ni_wlm_p <- function(ACE_Ni_wlm) {
  f <- summary(ACE_Ni_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wlm_p(ACE_Ni_wlm)

ACE_Ni_wlm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/(Ni_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wlm_p(ACE_Ni_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ni_wls_p <- function(ACE_Ni_wls) {
  f <- summary(ACE_Ni_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wls_p(ACE_Ni_wls)

ACE_Ni_wls_eqn <- function(df){
  m <- ACE_Ni_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wls_p(ACE_Ni_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ni_wls_wt_p <- function(ACE_Ni_wls_wt) {
  f <- summary(ACE_Ni_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wls_wt_p(ACE_Ni_wls_wt)

ACE_Ni_wls_wt_eqn <- function(df){
  m <- ACE_Ni_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wls_wt_p(ACE_Ni_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Ni_final <- ACE_Ni + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ni_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ni_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ni_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ni_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ni_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Ni/Ni_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ACE_Cu -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Cu_lm <- lm(Cu_ICP ~ Cu, data = ACE_dataset)
summary(ACE_Cu_lm)
glance(ACE_Cu_lm)
model_performance(ACE_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_OLS_summary.txt")
summary(ACE_Cu_lm)
glance(ACE_Cu_lm)
model_performance(ACE_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_lm) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Cu_wlm <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/(Cu_sd)^2)
summary(ACE_Cu_wlm)
glance(ACE_Cu_wlm)
model_performance(ACE_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_OLS_wt_summary.txt")
summary(ACE_Cu_wlm)
glance(ACE_Cu_wlm)
model_performance(ACE_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_lm) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Cu_model <- lm(Cu_ICP ~ Cu, data = ACE_dataset) # define model
ACE_Cu_wt <- 1 / lm(abs(ACE_Cu_model$residuals) ~ ACE_Cu_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Cu_wls <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt) #perform weighted least squares regression
# Checks
summary(ACE_Cu_wls) # summary stats
glance(ACE_Cu_wls) # summary stats including AIC
model_performance(ACE_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_WLS_summary.txt")
summary(ACE_Cu_wls)
glance(ACE_Cu_wls)
model_performance(ACE_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_wls) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Cu_model_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/Cu_sd^2) # define model
ACE_Cu_wt_wt <- 1 / lm(abs(ACE_Cu_model_wt$residuals) ~ ACE_Cu_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Cu_wls_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Cu_wls_wt) # summary stats
glance(ACE_Cu_wls_wt) # summary stats including AIC
model_performance(ACE_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_WLS_wt_summary.txt")
summary(ACE_Cu_wls_wt)
glance(ACE_Cu_wls_wt)
model_performance(ACE_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_wls_wt) # check for random efCucts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Cu"
theme_set(theme_classic(10))
ACE_Cu <- ggplot(ACE_dataset, aes(x = Cu, y = Cu_ICP)) +
  geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Cu_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Cu_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Cu_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Cu
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Cu_lm_p <- function(ACE_Cu_lm) {
  f <- summary(ACE_Cu_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_lm_p(ACE_Cu_lm)

ACE_Cu_lm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_lm_p(ACE_Cu_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Cu_wlm_p <- function(ACE_Cu_wlm) {
  f <- summary(ACE_Cu_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wlm_p(ACE_Cu_wlm)

ACE_Cu_wlm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/(Cu_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wlm_p(ACE_Cu_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Cu_wls_p <- function(ACE_Cu_wls) {
  f <- summary(ACE_Cu_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wls_p(ACE_Cu_wls)

ACE_Cu_wls_eqn <- function(df){
  m <- ACE_Cu_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wls_p(ACE_Cu_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Cu_wls_wt_p <- function(ACE_Cu_wls_wt) {
  f <- summary(ACE_Cu_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wls_wt_p(ACE_Cu_wls_wt)

ACE_Cu_wls_wt_eqn <- function(df){
  m <- ACE_Cu_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wls_wt_p(ACE_Cu_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Cu_final <- ACE_Cu + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Cu_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Cu_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Cu_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Cu_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Cu_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Cu/Cu_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# ACE_Zn -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Zn_lm <- lm(Zn_ICP ~ Zn, data = ACE_dataset)
summary(ACE_Zn_lm)
glance(ACE_Zn_lm)
model_performance(ACE_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_OLS_summary.txt")
summary(ACE_Zn_lm)
glance(ACE_Zn_lm)
model_performance(ACE_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_lm) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Zn_wlm <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/(Zn_sd)^2)
summary(ACE_Zn_wlm)
glance(ACE_Zn_wlm)
model_performance(ACE_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_OLS_wt_summary.txt")
summary(ACE_Zn_wlm)
glance(ACE_Zn_wlm)
model_performance(ACE_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_lm) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Zn_model <- lm(Zn_ICP ~ Zn, data = ACE_dataset) # define model
ACE_Zn_wt <- 1 / lm(abs(ACE_Zn_model$residuals) ~ ACE_Zn_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Zn_wls <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zn_wls) # summary stats
glance(ACE_Zn_wls) # summary stats including AIC
model_performance(ACE_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_WLS_summary.txt")
summary(ACE_Zn_wls)
glance(ACE_Zn_wls)
model_performance(ACE_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_wls) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Zn_model_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/Zn_sd^2) # define model
ACE_Zn_wt_wt <- 1 / lm(abs(ACE_Zn_model_wt$residuals) ~ ACE_Zn_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Zn_wls_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zn_wls_wt) # summary stats
glance(ACE_Zn_wls_wt) # summary stats including AIC
model_performance(ACE_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_WLS_wt_summary.txt")
summary(ACE_Zn_wls_wt)
glance(ACE_Zn_wls_wt)
model_performance(ACE_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_wls_wt) # check for random efZncts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Zn"
theme_set(theme_classic(10))
ACE_Zn <- ggplot(ACE_dataset, aes(x = Zn, y = Zn_ICP)) +
  geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Zn_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Zn_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Zn_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Zn
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zn_lm_p <- function(ACE_Zn_lm) {
  f <- summary(ACE_Zn_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_lm_p(ACE_Zn_lm)

ACE_Zn_lm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_lm_p(ACE_Zn_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zn_wlm_p <- function(ACE_Zn_wlm) {
  f <- summary(ACE_Zn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wlm_p(ACE_Zn_wlm)

ACE_Zn_wlm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/(Zn_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wlm_p(ACE_Zn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zn_wls_p <- function(ACE_Zn_wls) {
  f <- summary(ACE_Zn_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wls_p(ACE_Zn_wls)

ACE_Zn_wls_eqn <- function(df){
  m <- ACE_Zn_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wls_p(ACE_Zn_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zn_wls_wt_p <- function(ACE_Zn_wls_wt) {
  f <- summary(ACE_Zn_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wls_wt_p(ACE_Zn_wls_wt)

ACE_Zn_wls_wt_eqn <- function(df){
  m <- ACE_Zn_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wls_wt_p(ACE_Zn_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Zn_final <- ACE_Zn + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Zn_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Zn_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Zn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Zn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Zn_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zn/Zn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# ACE_Rb -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Rb_lm <- lm(Rb_ICP ~ Rb, data = ACE_dataset)
summary(ACE_Rb_lm)
glance(ACE_Rb_lm)
model_performance(ACE_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_OLS_summary.txt")
summary(ACE_Rb_lm)
glance(ACE_Rb_lm)
model_performance(ACE_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_lm) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Rb_wlm <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/(Rb_sd)^2)
summary(ACE_Rb_wlm)
glance(ACE_Rb_wlm)
model_performance(ACE_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_OLS_wt_summary.txt")
summary(ACE_Rb_wlm)
glance(ACE_Rb_wlm)
model_performance(ACE_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_lm) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Rb_model <- lm(Rb_ICP ~ Rb, data = ACE_dataset) # define model
ACE_Rb_wt <- 1 / lm(abs(ACE_Rb_model$residuals) ~ ACE_Rb_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Rb_wls <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt) #perform weighted least squares regression
# Checks
summary(ACE_Rb_wls) # summary stats
glance(ACE_Rb_wls) # summary stats including AIC
model_performance(ACE_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_WLS_summary.txt")
summary(ACE_Rb_wls)
glance(ACE_Rb_wls)
model_performance(ACE_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_wls) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Rb_model_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/Rb_sd^2) # define model
ACE_Rb_wt_wt <- 1 / lm(abs(ACE_Rb_model_wt$residuals) ~ ACE_Rb_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Rb_wls_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Rb_wls_wt) # summary stats
glance(ACE_Rb_wls_wt) # summary stats including AIC
model_performance(ACE_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_WLS_wt_summary.txt")
summary(ACE_Rb_wls_wt)
glance(ACE_Rb_wls_wt)
model_performance(ACE_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_wls_wt) # check for random efRbcts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Rb"
theme_set(theme_classic(10))
ACE_Rb <- ggplot(ACE_dataset, aes(x = Rb, y = Rb_ICP)) +
  geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Rb_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Rb_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Rb_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Rb
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Rb_lm_p <- function(ACE_Rb_lm) {
  f <- summary(ACE_Rb_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_lm_p(ACE_Rb_lm)

ACE_Rb_lm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_lm_p(ACE_Rb_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Rb_wlm_p <- function(ACE_Rb_wlm) {
  f <- summary(ACE_Rb_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wlm_p(ACE_Rb_wlm)

ACE_Rb_wlm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/(Rb_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wlm_p(ACE_Rb_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Rb_wls_p <- function(ACE_Rb_wls) {
  f <- summary(ACE_Rb_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wls_p(ACE_Rb_wls)

ACE_Rb_wls_eqn <- function(df){
  m <- ACE_Rb_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wls_p(ACE_Rb_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Rb_wls_wt_p <- function(ACE_Rb_wls_wt) {
  f <- summary(ACE_Rb_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wls_wt_p(ACE_Rb_wls_wt)

ACE_Rb_wls_wt_eqn <- function(df){
  m <- ACE_Rb_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wls_wt_p(ACE_Rb_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Rb_final <- ACE_Rb + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Rb_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Rb_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Rb_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Rb_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Rb_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Rb/Rb_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# ACE_Sr -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Sr_lm <- lm(Sr_ICP ~ Sr, data = ACE_dataset)
summary(ACE_Sr_lm)
glance(ACE_Sr_lm)
model_performance(ACE_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_OLS_summary.txt")
summary(ACE_Sr_lm)
glance(ACE_Sr_lm)
model_performance(ACE_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_lm) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Sr_wlm <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/(Sr_sd)^2)
summary(ACE_Sr_wlm)
glance(ACE_Sr_wlm)
model_performance(ACE_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_OLS_wt_summary.txt")
summary(ACE_Sr_wlm)
glance(ACE_Sr_wlm)
model_performance(ACE_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_lm) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Sr_model <- lm(Sr_ICP ~ Sr, data = ACE_dataset) # define model
ACE_Sr_wt <- 1 / lm(abs(ACE_Sr_model$residuals) ~ ACE_Sr_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Sr_wls <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt) #perform weighted least squares regression
# Checks
summary(ACE_Sr_wls) # summary stats
glance(ACE_Sr_wls) # summary stats including AIC
model_performance(ACE_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_WLS_summary.txt")
summary(ACE_Sr_wls)
glance(ACE_Sr_wls)
model_performance(ACE_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_wls) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Sr_model_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/Sr_sd^2) # define model
ACE_Sr_wt_wt <- 1 / lm(abs(ACE_Sr_model_wt$residuals) ~ ACE_Sr_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Sr_wls_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Sr_wls_wt) # summary stats
glance(ACE_Sr_wls_wt) # summary stats including AIC
model_performance(ACE_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_WLS_wt_summary.txt")
summary(ACE_Sr_wls_wt)
glance(ACE_Sr_wls_wt)
model_performance(ACE_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_wls_wt) # check for random efSrcts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Sr"
theme_set(theme_classic(10))
ACE_Sr <- ggplot(ACE_dataset, aes(x = Sr, y = Sr_ICP)) +
  geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Sr_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Sr_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Sr_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Sr
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Sr_lm_p <- function(ACE_Sr_lm) {
  f <- summary(ACE_Sr_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_lm_p(ACE_Sr_lm)

ACE_Sr_lm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_lm_p(ACE_Sr_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Sr_wlm_p <- function(ACE_Sr_wlm) {
  f <- summary(ACE_Sr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wlm_p(ACE_Sr_wlm)

ACE_Sr_wlm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/(Sr_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wlm_p(ACE_Sr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Sr_wls_p <- function(ACE_Sr_wls) {
  f <- summary(ACE_Sr_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wls_p(ACE_Sr_wls)

ACE_Sr_wls_eqn <- function(df){
  m <- ACE_Sr_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wls_p(ACE_Sr_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Sr_wls_wt_p <- function(ACE_Sr_wls_wt) {
  f <- summary(ACE_Sr_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wls_wt_p(ACE_Sr_wls_wt)

ACE_Sr_wls_wt_eqn <- function(df){
  m <- ACE_Sr_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wls_wt_p(ACE_Sr_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Sr_final <- ACE_Sr + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Sr_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Sr_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Sr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Sr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Sr_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Sr/Sr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# ACE_Zr -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Zr_lm <- lm(Zr_ICP ~ Zr, data = ACE_dataset)
summary(ACE_Zr_lm)
glance(ACE_Zr_lm)
model_performance(ACE_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_OLS_summary.txt")
summary(ACE_Zr_lm)
glance(ACE_Zr_lm)
model_performance(ACE_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_lm) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_Zr_wlm <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/(Zr_sd)^2)
summary(ACE_Zr_wlm)
glance(ACE_Zr_wlm)
model_performance(ACE_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_OLS_wt_summary.txt")
summary(ACE_Zr_wlm)
glance(ACE_Zr_wlm)
model_performance(ACE_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_lm) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_Zr_model <- lm(Zr_ICP ~ Zr, data = ACE_dataset) # define model
ACE_Zr_wt <- 1 / lm(abs(ACE_Zr_model$residuals) ~ ACE_Zr_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Zr_wls <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zr_wls) # summary stats
glance(ACE_Zr_wls) # summary stats including AIC
model_performance(ACE_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_WLS_summary.txt")
summary(ACE_Zr_wls)
glance(ACE_Zr_wls)
model_performance(ACE_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_wls) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Zr_model_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/Zr_sd^2) # define model
ACE_Zr_wt_wt <- 1 / lm(abs(ACE_Zr_model_wt$residuals) ~ ACE_Zr_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Zr_wls_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zr_wls_wt) # summary stats
glance(ACE_Zr_wls_wt) # summary stats including AIC
model_performance(ACE_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_WLS_wt_summary.txt")
summary(ACE_Zr_wls_wt)
glance(ACE_Zr_wls_wt)
model_performance(ACE_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_wls_wt) # check for random efZrcts - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "Zr"
theme_set(theme_classic(10))
ACE_Zr <- ggplot(ACE_dataset, aes(x = Zr, y = Zr_ICP)) +
  geom_errorbar(aes(ymin=Zr_ICP-Zr_ICP_sd, ymax=Zr_ICP+Zr_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=Zr-Zr_sd, xmax=Zr+Zr_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/Zr_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_Zr_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_Zr_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_Zr
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zr_lm_p <- function(ACE_Zr_lm) {
  f <- summary(ACE_Zr_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_lm_p(ACE_Zr_lm)

ACE_Zr_lm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_lm_p(ACE_Zr_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zr_wlm_p <- function(ACE_Zr_wlm) {
  f <- summary(ACE_Zr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wlm_p(ACE_Zr_wlm)

ACE_Zr_wlm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/(Zr_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wlm_p(ACE_Zr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zr_wls_p <- function(ACE_Zr_wls) {
  f <- summary(ACE_Zr_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wls_p(ACE_Zr_wls)

ACE_Zr_wls_eqn <- function(df){
  m <- ACE_Zr_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wls_p(ACE_Zr_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zr_wls_wt_p <- function(ACE_Zr_wls_wt) {
  f <- summary(ACE_Zr_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wls_wt_p(ACE_Zr_wls_wt)

ACE_Zr_wls_wt_eqn <- function(df){
  m <- ACE_Zr_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wls_wt_p(ACE_Zr_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_Zr_final <- ACE_Zr + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Zr_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Zr_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Zr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Zr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_Zr_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/Zr/Zr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# ACE_DM -----------------------------------------------------------------------

# Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_DM_lm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset)
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_OLS_summary.txt")
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# Weighted OLS linear model & checks
ACE_DM_wlm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/(coh_inc_sd)^2)
summary(ACE_DM_wlm)
glance(ACE_DM_wlm)
model_performance(ACE_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_OLS_wt_summary.txt")
summary(ACE_DM_wlm)
glance(ACE_DM_wlm)
model_performance(ACE_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# Unweighted Linear Regression (WLS) model
ACE_DM_model <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset) # define model
ACE_DM_wt <- 1 / lm(abs(ACE_DM_model$residuals) ~ ACE_DM_model$fitted.values)$fitted.values^2 #define weights to use
ACE_DM_wls <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=ACE_DM_wt) #perform weighted least squares regression
# Checks
summary(ACE_DM_wls) # summary stats
glance(ACE_DM_wls) # summary stats including AIC
model_performance(ACE_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_WLS_summary.txt")
summary(ACE_DM_wls)
glance(ACE_DM_wls)
model_performance(ACE_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_wls) # Performance package summary check for heteroscedasticity
icc(ACE_DM_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_DM_model_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/coh_inc_sd^2) # define model
ACE_DM_wt_wt <- 1 / lm(abs(ACE_DM_model_wt$residuals) ~ ACE_DM_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_DM_wls_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=ACE_DM_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_DM_wls_wt) # summary stats
glance(ACE_DM_wls_wt) # summary stats including AIC
model_performance(ACE_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_WLS_wt_summary.txt")
summary(ACE_DM_wls_wt)
glance(ACE_DM_wls_wt)
model_performance(ACE_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_DM_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)

# Linear model summary elemental plots

element_title <- "coh_inc"
theme_set(theme_classic(10))
ACE_DM <- ggplot(ACE_dataset, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", colour = "lightblue") +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
              aes(weight = 1/coh_inc_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "dashed", aes(weight = ACE_DM_wt), colour="darkgrey") + # WLS regression
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              aes(weight = ACE_DM_wt_wt), colour="black") + # weighted WLS regression
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
  labs(x = paste(element_title, " [XRF-CS]", itrax_dataset), y = paste0(element_title, " [ICPMS] (ppm)")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset))
ACE_DM
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_DM_lm_p <- function(ACE_DM_lm) {
  f <- summary(ACE_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_lm_p(ACE_DM_lm)

ACE_DM_lm_eqn <- function(df){
  m <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_lm_p(ACE_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_DM_wlm_p <- function(ACE_DM_wlm) {
  f <- summary(ACE_DM_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_wlm_p(ACE_DM_wlm)

ACE_DM_wlm_eqn <- function(df){
  m <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/(coh_inc_sd)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_wlm_p(ACE_DM_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_DM_wls_p <- function(ACE_DM_wls) {
  f <- summary(ACE_DM_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_wls_p(ACE_DM_wls)

ACE_DM_wls_eqn <- function(df){
  m <- ACE_DM_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_wls_p(ACE_DM_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_DM_wls_wt_p <- function(ACE_DM_wls_wt) {
  f <- summary(ACE_DM_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_wls_wt_p(ACE_DM_wls_wt)

ACE_DM_wls_wt_eqn <- function(df){
  m <- ACE_DM_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_wls_wt_p(ACE_DM_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_DM_final <- ACE_DM + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_DM_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_DM_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_DM_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "darkgrey", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "lightblue", hjust = -0.15, vjust = 5, size = 3)
ACE_DM_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/DM/DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# -------------------------------------------------------------------------
# Summary 4x3 matrix plots of ITRAX-ACF matched elements

# OLS & WLS summary - unweighted stats on plot
ggarrange(ACE_K_final, ACE_Ca_final, ACE_Ti_final,
          ACE_Mn_final, ACE_Fe_final, ACE_Co_final,
          ACE_Ni_final, ACE_Cu_final, ACE_Zn_final,
          ACE_Rb_final, ACE_Sr_final, ACE_Zr_final,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/OLS_WLS_inc/ACE_OLS_WLS_Summary.pdf",
       height = c(50), width = c(50), dpi = 600, units = "cm")
# -------------------------------------------------------------------------




END


# OLS Linear models - forced through origin -------------------------------

# K unweighted ---------------------------
ACE_K_origin <- ggplot(ACE_all, aes(x = K, y = K_ICP)) +
  geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=.1, colour = "grey") +
  geom_smooth(method=glm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/K_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "K (XRF-CS) / cps", y = "K (ICPMS) / ppm" ) +
  ggtitle("ACE: K (ICP SD weighted = blue; unweighted = black)")
ACE_K_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/K_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# Mn unweighted ---------------------------
ACE_Mn_origin <- ggplot(ACE_all, aes(x = Mn, y = Mn_ICP)) +
  geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Mn_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Mn (XRF-CS) / cps", y = "Mn (ICPMS) / ppm" ) +
  ggtitle("ACE: Mn (ICP SD weighted = blue; unweighted = black)")
ACE_Mn_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Mn_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Ti 

# Ti unweighted ---------------------------
ACE_Ti_origin <- ggplot(ACE_all, aes(x = Ti, y = Ti_ICP)) +
  geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Ti_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Ti (XRF-CS) / cps", y = "Ti (ICPMS) / ppm" ) +
  ggtitle("ACE: Ti (ICP SD weighted = blue; unweighted = black)")
ACE_Ti_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Ti_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Mn

# Mn unweighted ---------------------------
ACE_Mn_origin <- ggplot(ACE_all, aes(x = Mn, y = Mn_ICP)) +
  geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Mn_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Mn (XRF-CS) / cps", y = "Mn (ICPMS) / ppm" ) +
  ggtitle("ACE: Mn (ICP SD weighted = blue; unweighted = black)")
ACE_Mn_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/cpsorigin//Mn_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Fe unweighted ---------------------------
ACE_Fe_origin <- ggplot(ACE_all, aes(x = Fe, y = Fe_ICP)) +
  geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Fe_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Fe (XRF-CS) / cps", y = "Fe (ICPMS) / ppm" ) +
  ggtitle("ACE: Fe (ICP SD weighted = blue; unweighted = black)")
ACE_Fe_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Fe_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Co

# Co unweighted ---------------------------
ACE_Co_origin <- ggplot(ACE_all, aes(x = Co, y = Co_ICP)) +
  geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Co_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Co (XRF-CS) / cps", y = "Co (ICPMS) / ppm" ) +
  ggtitle("ACE: Co (ICP SD weighted = blue; unweighted = black)")
ACE_Co_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Co_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Ni

# Ni unweighted ---------------------------
ACE_Ni_origin <- ggplot(ACE_all, aes(x = Ni, y = Ni_ICP)) +
  geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Ni_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Ni (XRF-CS) / cps", y = "Ni (ICPMS) / ppm" ) +
  ggtitle("ACE: Ni (ICP SD weighted = blue; unweighted = black)")
ACE_Ni_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Ni_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Cu

# Cu unweighted ---------------------------
ACE_Cu_origin <- ggplot(ACE_all, aes(x = Cu, y = Cu_ICP)) +
  geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Cu_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Cu (XRF-CS) / cps", y = "Cu (ICPMS) / ppm" ) +
  ggtitle("ACE: Cu (ICP SD weighted = blue; unweighted = black)")
ACE_Cu_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Cu_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Zn

# Zn unweighted ---------------------------
ACE_Zn_origin <- ggplot(ACE_all, aes(x = Zn, y = Zn_ICP)) +
  geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Zn_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Zn (XRF-CS) / cps", y = "Zn (ICPMS) / ppm" ) +
  ggtitle("ACE: Zn (ICP SD weighted = blue; unweighted = black)")
ACE_Zn_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Zn_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Rb

# Rb unweighted ---------------------------
ACE_Rb_origin <- ggplot(ACE_all, aes(x = Rb, y = Rb_ICP)) +
  geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Rb_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Rb (XRF-CS) / cps", y = "Rb (ICPMS) / ppm" ) +
  ggtitle("ACE: Rb (ICP SD weighted = blue; unweighted = black)")
ACE_Rb_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Rb_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Sr

# Sr unweighted ---------------------------
ACE_Sr_origin <- ggplot(ACE_all, aes(x = Sr, y = Sr_ICP)) +
  geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Sr_ICP_sd)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "Sr (XRF-CS) / cps", y = "Sr (ICPMS) / ppm" ) +
  ggtitle("ACE: Sr (ICP SD weighted = blue; unweighted = black)")
ACE_Sr_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Sr_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# DM

# DM unweighted ---------------------------
ACE_DM_origin <- ggplot(ACE_all, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, formula=y~x-1, aes(weight = 1/dry_mass_err)) +
  stat_poly_line(formula=y~x-1, colour = "black") + # unweighted line
  stat_poly_eq(formula=y~x-1, use_label(c("R2", "adj.R2", "p", "eq", "f")), colour = "black", cex = 4) + # unweighted stats displayed
  stat_poly_eq(use_label("n"), label.x = "right", label.y = "top", colour = "black") +
  geom_point() +
  scale_shape_manual(values = c(21)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8), 
        legend.justification = c(1,1),
        legend.position = c(1,1), 
        axis.text=element_text(size=10, colour = "black"), 
        axis.title=element_text(size=10, colour = "black"),
        title = element_text(size=10, colour = "black")) +
  labs(x = "coh/inc (XRF-CS) / cps", y = "DM (ICPMS) / ppm" ) +
  ggtitle("ACE: DM (ICP SD weighted = blue; unweighted = black)")
ACE_DM_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/DM_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")





# -------------------------------------------------------------------------
# Performance tests & stats - forced through origin - lm and wlm ------------------------

# K

#Performance: K unweighted lm stats forced through origin
ACE_K_lm_origin <- lm(K_ICP ~ K-1, data = ACE_all)
summary(ACE_K_lm_origin) # summary stats
check_model(ACE_K_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/K_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_K_lm_origin_mp <- model_performance(ACE_K_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_K_lm_origin_mp)
bptest(ACE_K_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_K_lm)

#Performance: K weighted lm forced through origin
ACE_K_wlm_origin <- lm(K_ICP ~ K-1, weight = 1/K_ICP_sd, data = ACE_all)
summary(ACE_K_wlm_origin) # summary stats

check_model(ACE_K_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/K_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_K_wlm_origin_mp <- model_performance(ACE_K_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_K_wlm_origin_mp)
bptest(ACE_K_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_K_wlm)

# Mn

#Performance: Mn unweighted lm stats forced through origin
ACE_Mn_lm_origin <- lm(Mn_ICP ~ Mn-1, data = ACE_all)
summary(ACE_Mn_lm_origin) # summary stats
check_model(ACE_Mn_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Mn_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Mn_lm_origin_mp <- model_performance(ACE_Mn_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Mn_lm_origin_mp)
bptest(ACE_Mn_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm)

#Performance: Mn weighted lm forced through origin
ACE_Mn_wlm_origin <- lm(Mn_ICP ~ Mn-1, weight = 1/Mn_ICP_sd, data = ACE_all)
summary(ACE_Mn_wlm_origin) # summary stats
check_model(ACE_Mn_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Mn_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Mn_wlm_origin_mp <- model_performance(ACE_Mn_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Mn_wlm_origin_mp)
bptest(ACE_Mn_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wlm)

# Ti

#Performance: Ti unweighted lm stats forced through origin
ACE_Ti_lm_origin <- lm(Ti_ICP ~ Ti-1, data = ACE_all)
summary(ACE_Ti_lm_origin) # summary stats
check_model(ACE_Ti_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Ti_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Ti_lm_origin_mp <- model_performance(ACE_Ti_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Ti_lm_origin_mp)
bptest(ACE_Ti_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_lm)

#Performance: Ti weighted lm forced through origin
ACE_Ti_wlm_origin <- lm(Ti_ICP ~ Ti-1, weight = 1/Ti_ICP_sd, data = ACE_all)
summary(ACE_Ti_wlm_origin) # summary stats
check_model(ACE_Ti_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Ti_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Ti_wlm_origin_mp <- model_performance(ACE_Ti_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Ti_wlm_origin_mp)
bptest(ACE_Ti_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_wlm)

# Mn

#Performance: Mn unweighted lm stats forced through origin
ACE_Mn_lm_origin <- lm(Mn_ICP ~ Mn-1, data = ACE_all)
summary(ACE_Mn_lm_origin) # summary stats
check_model(ACE_Mn_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Mn_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Mn_lm_origin_mp <- model_performance(ACE_Mn_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Mn_lm_origin_mp)
bptest(ACE_Mn_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm)

#Performance: Mn weighted lm forced through origin
ACE_Mn_wlm_origin <- lm(Mn_ICP ~ Mn-1, weight = 1/Mn_ICP_sd, data = ACE_all)
summary(ACE_Mn_wlm_origin) # summary stats
check_model(ACE_Mn_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Mn_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Mn_wlm_origin_mp <- model_performance(ACE_Mn_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Mn_wlm_origin_mp)
bptest(ACE_Mn_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wlm)

# Fe

#Performance: Fe unweighted lm stats forced through origin
ACE_Fe_lm_origin <- lm(Fe_ICP ~ Fe-1, data = ACE_all)
summary(ACE_Fe_lm_origin) # summary stats
check_model(ACE_Fe_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Fe_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Fe_lm_origin_mp <- model_performance(ACE_Fe_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Fe_lm_origin_mp)
bptest(ACE_Fe_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_lm)

#Performance: Fe weighted lm forced through origin
ACE_Fe_wlm_origin <- lm(Fe_ICP ~ Fe-1, weight = 1/Fe_ICP_sd, data = ACE_all)
summary(ACE_Fe_wlm_origin) # summary stats
check_model(ACE_Fe_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Fe_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Fe_wlm_origin_mp <- model_performance(ACE_Fe_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Fe_wlm_origin_mp)
bptest(ACE_Fe_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_wlm)

# Co 

#Performance: Co unweighted lm stats forced through origin
ACE_Co_lm_origin <- lm(Co_ICP ~ Co-1, data = ACE_all)
summary(ACE_Co_lm_origin) # summary stats
check_model(ACE_Co_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Co_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Co_lm_origin_mp <- model_performance(ACE_Co_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Co_lm_origin_mp)
bptest(ACE_Co_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Co_lm)

#Performance: Co weighted lm forced through origin
ACE_Co_wlm_origin <- lm(Co_ICP ~ Co-1, weight = 1/Co_ICP_sd, data = ACE_all)
summary(ACE_Co_wlm_origin) # summary stats
check_model(ACE_Co_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Co_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Co_wlm_origin_mp <- model_performance(ACE_Co_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Co_wlm_origin_mp)
bptest(ACE_Co_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Co_wlm)

# Ni

#Performance: Ni unweighted lm stats forced through origin
ACE_Ni_lm_origin <- lm(Ni_ICP ~ Ni-1, data = ACE_all)
summary(ACE_Ni_lm_origin) # summary stats
check_model(ACE_Ni_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Ni_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Ni_lm_origin_mp <- model_performance(ACE_Ni_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Ni_lm_origin_mp)
bptest(ACE_Ni_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_lm)

#Performance: Ni weighted lm forced through origin
ACE_Ni_wlm_origin <- lm(Ni_ICP ~ Ni-1, weight = 1/Ni_ICP_sd, data = ACE_all)
summary(ACE_Ni_wlm_origin) # summary stats
check_model(ACE_Ni_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Ni_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Ni_wlm_origin_mp <- model_performance(ACE_Ni_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Ni_wlm_origin_mp)
bptest(ACE_Ni_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_wlm)

# Cu

#Performance: Cu unweighted lm stats forced through origin
ACE_Cu_lm_origin <- lm(Cu_ICP ~ Cu-1, data = ACE_all)
summary(ACE_Cu_lm_origin) # summary stats
check_model(ACE_Cu_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Cu_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Cu_lm_origin_mp <- model_performance(ACE_Cu_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Cu_lm_origin_mp)
bptest(ACE_Cu_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_lm)

#Performance: Cu weighted lm forced through origin
ACE_Cu_wlm_origin <- lm(Cu_ICP ~ Cu-1, weight = 1/Cu_ICP_sd, data = ACE_all)
summary(ACE_Cu_wlm_origin) # summary stats
check_model(ACE_Cu_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Cu_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Cu_wlm_origin_mp <- model_performance(ACE_Cu_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Cu_wlm_origin_mp)
bptest(ACE_Cu_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_wlm)

# Zn

#Performance: Zn unweighted lm stats forced through origin
ACE_Zn_lm_origin <- lm(Zn_ICP ~ Zn-1, data = ACE_all)
summary(ACE_Zn_lm_origin) # summary stats
check_model(ACE_Zn_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Zn_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Zn_lm_origin_mp <- model_performance(ACE_Zn_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Zn_lm_origin_mp)
bptest(ACE_Zn_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_lm)

#Performance: Zn weighted lm forced through origin
ACE_Zn_wlm_origin <- lm(Zn_ICP ~ Zn-1, weight = 1/Zn_ICP_sd, data = ACE_all)
summary(ACE_Zn_wlm_origin) # summary stats
check_model(ACE_Zn_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Zn_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Zn_wlm_origin_mp <- model_performance(ACE_Zn_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Zn_wlm_origin_mp)
bptest(ACE_Zn_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_wlm)

# Rb

#Performance: Rb unweighted lm stats forced through origin
ACE_Rb_lm_origin <- lm(Rb_ICP ~ Rb-1, data = ACE_all)
summary(ACE_Rb_lm_origin) # summary stats
check_model(ACE_Rb_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Rb_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Rb_lm_origin_mp <- model_performance(ACE_Rb_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Rb_lm_origin_mp)
bptest(ACE_Rb_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_lm)

#Performance: Rb weighted lm forced through origin
ACE_Rb_wlm_origin <- lm(Rb_ICP ~ Rb-1, weight = 1/Rb_ICP_sd, data = ACE_all)
summary(ACE_Rb_wlm_origin) # summary stats
check_model(ACE_Rb_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Rb_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Rb_wlm_origin_mp <- model_performance(ACE_Rb_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Rb_wlm_origin_mp)
bptest(ACE_Rb_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_wlm)

# Sr

#Performance: Sr unweighted lm stats forced through origin
ACE_Sr_lm_origin <- lm(Sr_ICP ~ Sr-1, data = ACE_all)
summary(ACE_Sr_lm_origin) # summary stats
check_model(ACE_Sr_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Sr_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Sr_lm_origin_mp <- model_performance(ACE_Sr_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Sr_lm_origin_mp)
bptest(ACE_Sr_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_lm)

#Performance: Sr weighted lm forced through origin
ACE_Sr_wlm_origin <- lm(Sr_ICP ~ Sr-1, weight = 1/Sr_ICP_sd, data = ACE_all)
summary(ACE_Sr_wlm_origin) # summary stats
check_model(ACE_Sr_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/Sr_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_Sr_wlm_origin_mp <- model_performance(ACE_Sr_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_Sr_wlm_origin_mp)
bptest(ACE_Sr_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_wlm)

# DM

#Performance: DM unweighted lm stats forced through origin
ACE_DM_lm_origin <- lm(dry_mass_pc ~ coh_inc-1, data = ACE_all)
summary(ACE_DM_lm_origin) # summary stats
check_model(ACE_DM_lm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/DM_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_DM_lm_origin_mp <- model_performance(ACE_DM_lm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_DM_lm_origin_mp)
bptest(ACE_DM_lm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm)

#Performance: DM weighted lm forced through origin
ACE_DM_wlm_origin <- lm(dry_mass_pc ~ coh_inc-1, weight = 1/dry_mass_err, data = ACE_all)
summary(ACE_DM_wlm_origin) # summary stats
check_model(ACE_DM_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/origin/DM_origin_performance_wlm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_DM_wlm_origin_mp <- model_performance(ACE_DM_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_DM_wlm_origin_mp)
bptest(ACE_DM_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_DM_wlm)





# Summary 4x3 matrix plots of ITRAX-ACF matched variables --------
# OLS summary - forced through origin - unweighted stats on plot
ggarrange(ACE_K_origin, ACE_Mn_origin, ACE_Ti_origin, 
          ACE_Mn_origin, ACE_Fe_origin, ACE_Co_origin, 
          ACE_Ni_origin, ACE_Cu_origin, ACE_Zn_origin,
          ACE_Rb_origin, ACE_Sr_origin, ACE_DM_origin,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/inc/ACE_OLS_Summary_origin_lm_wlm.pdf", 
       height = c(50), width = c(50), dpi = 600, units = "cm")





