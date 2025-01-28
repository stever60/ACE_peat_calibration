
# -------------------------------------------------------------------------

# ITRAX-ICPMS Data Matching, Correlation & Regression for centred log ratio dataset

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

# -------------------------------------------------------------------------
# Matching - import matched itrax-ICPMS datafiles from each site -------------------

ACE_itrax_BI10 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/BI10/cps/BI10_xrf_icp_matched.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_HER42PB <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/HER42PB/cps/HER42PB_xrf_icp_matched.csv", 
                              col_names = TRUE, skip = 0)
ACE_itrax_KER1 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER1/cps/KER1_xrf_icp_matched.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_KER3 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/KER3/cps/KER3_xrf_icp_matched.csv", 
                           col_names = TRUE, skip = 0)
ACE_itrax_PB1 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/PB1/cps/PB1_xrf_icp_matched.csv", 
                          col_names = TRUE, skip = 0)
ACE_itrax_POB4 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/Sites/POB4/cps/POB4_xrf_icp_matched.csv", 
                           col_names = TRUE, skip = 0)

# Combine matched output from each site
#inc_cols <- c('Mo_inc', 'Mo_inc_sd')
ACE_xrf_icp_matched <- bind_rows(ACE_itrax_BI10, 
                                 ACE_itrax_HER42PB, 
                                 ACE_itrax_KER1, 
                                 ACE_itrax_KER3, 
                                 ACE_itrax_PB1,
                                 ACE_itrax_POB4) %>% 
  #select(-all_of(inc_cols)) %>% 
  print()

write.csv(ACE_xrf_icp_matched,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_xrf_icp_matched_cps.csv", row.names = FALSE)

# START HERE if ACE matched already written ----------------------------------------------------------

# Element lists

# Defined from ACF analysis of ACE09 composite ITRAX dataframe analysis file:
# https://www.dropbox.com/s/rhtlkp6uwp71ryc/ACE_itrax_Composite.R?dl=0

# ICP elements as defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")

# XRF elements defined by Francois & ITRAX acf
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")

# ICP elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

xrf_icp_Elements_min <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                          "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                          "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                          "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                          "Mo_inc", "Mo_coh", "coh_inc", "dry_mass_pc")

# define adjacent matched elements for correlation
xrf_icp_Elements_min <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                          "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                          "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                          "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                          "Mo_inc", "Mo_coh", "coh_inc", "dry_mass_pc")

# define as adjacent matched elements for clr transformation
xrf_icp_Elements_min1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                           "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                           "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                           "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")

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
                                     "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP",
                                     "coh_inc", "dry_mass_pc")

# key elements
xrf_icp_Elements_key <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP",
                          "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "dry_mass_pc")

# key elements_reduced
xrf_icp_Elements_key_reduced <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")


# Import existing ACE cps matched file ----------------------------------------------

ACE_xrf_icp_matched <-read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_xrf_icp_matched_cps.csv")
is.na(ACE_xrf_icp_matched)<-sapply(ACE_xrf_icp_matched, is.infinite) # replace any infinite values with NA

# ACE matched - all sites - import cps matched file to make clr from 
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

# clr - make clr file of while dataset after replacing zeros with half min value and using using composition package 
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

# Bind into single file 
ACE_all <- bind_cols(ACE_all_text, ACE_all_clr_itrax, ACE_all_clr_itrax_sd, 
                     ACE_all_clr_icp, ACE_all_clr_icp_sd, ACE_all_rest) %>%
  print()
write.csv(ACE_all,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_all.csv", row.names = FALSE)
write.csv(ACE_all,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_xrf_icp_matched_clr.csv", row.names = FALSE)

# ACE matched - POB4 removed for regression analysis in cps runs
ACE_no_POB4 <- ACE_all %>% 
  filter(!Site =="POB4") %>% 
  print()
write.csv(ACE_no_POB4,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_no_POB4.csv", row.names = FALSE)

# ACE matched - POB4 & PB1 removed
ACE_no_POB4_PB1 <- ACE_all %>% 
  filter(!Site =="POB4") %>% 
  filter(!Site =="PB1") %>% 
  print()
write.csv(ACE_no_POB4_PB1,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_no_POB4_PB1.csv", row.names = FALSE)

ACE_no_PB1 <- ACE_all %>% # inc and %cps only 
  filter(!Site =="PB1") %>% 
  print()
write.csv(ACE_no_PB1,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_no_PB1.csv", row.names = FALSE)

# Create summary stats table ----------------------------------------------

ACE_all_stats <- ACE_all %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_all_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_all_stats.csv", row.names = FALSE)

ACE_no_POB4_stats <- ACE_no_POB4 %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_no_POB4_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_no_POB4_stats.csv", row.names = FALSE)

ACE_no_POB4_PB1_stats <- ACE_no_POB4_PB1 %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_no_POB4_PB1_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_no_POB4_PB1_stats.csv", row.names = FALSE)

ACE_no_PB1_stats <- ACE_no_PB1 %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_no_PB1_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_no_PB1_stats.csv", row.names = FALSE)

# Correlation & OLS linear modelling set up --------------------------------------------------

# Import datasets - needs ACE_all, ACE_LM1, ACE_LM2 to run - need to define ACE_dataset as well
#ACE_all <- ACE_all
ACE_all <- ACE_no_POB4
#ACE_all <- ACE_no_POB4_PB1

#ACE_LM1 <- ACE_all
ACE_LM1 <- ACE_no_POB4
#ACE_LM1 <- ACE_no_POB4_PB1
#ACE_LM1 <- ACE_no_PB1

#ACE_LM2 <- ACE_no_POB4
ACE_LM2 <- ACE_no_POB4_PB1
#ACE_LM2 <- ACE_no_PB1

# Convert Site to use as a grouping variable
ACE_all$Site <- as.factor(ACE_all$Site)
ACE_LM1$Site <- as.factor(ACE_LM1$Site)
ACE_LM2$Site <- as.factor(ACE_LM2$Site)

# Correlation matrices --------------------------------------------------

# ITRAX
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,acf_icp_Elements_min], method = c("everything", "pearson"), 
       size = 5, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 5)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ICP
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,icp_Elements_min], method = c("everything", "pearson"), 
       size = 5, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 5)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_ICP_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP - key elements correlations
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_min], method = c("everything", "pearson"), # use xrf_icp_Elements_min for all matched elements
       size = 3, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 3)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP - key elements correlations
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_key], method = c("everything", "pearson"), # use xrf_icp_Elements_min for all matched elements
       size = 7, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 7)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr_matrix_key.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP - key elements reduced correlations
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_key_reduced], method = c("everything", "pearson"), # use xrf_icp_Elements_min for all matched elements
       size = 11, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 11)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr_matrix_key_reduced.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation & density matrices ----------------------------------------

# *** if the p-value is < 0.001, ** if the p-value is < 0.01, * if the p-value is < 0.05, . if the p-value is < 0.10

# ITRAX
theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = acf_icp_Elements_min, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE XRF-CS: Correlation-density plot - clr")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_Corr-den_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ICP
theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = icp_Elements_min, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE ICPMS: Correlation-density plot - clr")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_ICP_Corr-den_matrix.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ITRAX & ICP - key elements correlations
theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = xrf_icp_Elements_min, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE XRF-CS & ICPMS: Correlation-density plot - clr")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr-den_matrix.pdf", 
       height = c(15), width = c(15), dpi = 600, units = "cm")

theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE XRF-CS & ICPMS: Correlation-density plot - clr")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr-den_matrix_key.pdf", 
       height = c(15), width = c(15), dpi = 600, units = "cm")

# ITRAX & ICP - key elements correlations and regressions
# Function to plot each site as different colour usign jco scheme
theme_set(theme_bw(base_size=8))
palette_set <- "jco" # or "jco", or "npg", "uchicago" - set up colour scheme for site plots
cor_plot <- function(data, mapping, ...) {
  
  ggplot(data, mapping) + 
    geom_point(size = 0.7) +
    geom_smooth(formula = y~x, method = lm, color = "black") + 
    scale_color_jco()
}

# Run ggpairs - black regression line all - correlation all sites
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, 
        title="ACE XRF-CS & ICPMS: Correlation-density plot for clr",
        upper = list(
          mapping = aes(color=Site, palette = palette_set),
          continuous = wrap(ggally_cor, size = 3)
        ),
        lower=list(
          mapping = aes(color=Site, palette = palette_set, alpha = 0.5),
          continuous = cor_plot
        )
)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr-den_matrix_key_reg1.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs - black regression line all - correlation overall - key
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key,
        title="ACE XRF-CS & ICPMS: Correlation-density plot - clr",
        upper = list(
          continuous = wrap(ggally_cor, size = 5)
        ),
        lower=list(
          mapping = aes(color=Site, palette = palette_set, alpha = 0.5),
          continuous = cor_plot
        )
)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr-den_matrix_key_reg2.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs - black regression line all - correlation overall - key reduced
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key_reduced,
        title="ACE XRF-CS & ICPMS: Correlation-density plot - clr",
        upper = list(
          continuous = wrap(ggally_cor, size = 5)
        ),
        lower=list(
          mapping = aes(color=Site, palette = palette_set, alpha = 0.5),
          continuous = cor_plot
        )
)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr-den_matrix_key_reduced_reg2.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs -  regression line all sites - correlation all sites
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, upper = list(continuous = wrap("cor", size = 3)),
        #lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        lower = list(continuous = "smooth"),
        aes(color = Site,  # Color by group (cat. variable)
            alpha = 0.5),     # Transparency
        title="ACE XRF-CS & ICPMS: Correlation-density plot for clr")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr-den_matrix_key_reg3.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs - no regression - correlation all sites
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=1)),
        #lower = list(continuous = "smooth"),
        aes(color = Site,  # Color by group (cat. variable)
            alpha = 0.5),     # Transparency
        title="ACE XRF-CS & ICPMS: Correlation-density plot for clr")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_itrax_ICP_Corr-den_matrix_site_key_plots.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation & data plots -----------------------------------------------------------

# Define titles & labels for plotting

XRF_title <- " clr XRF-CS"
ICP_title <- " clr ICPMS"
correlation_title <- "clr Correlation"
method_title <- "clr"
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
  ggtitle(paste("ACE (OLS): ", element_title, method_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=10, face="bold"))
K_corr1
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(K_XRF_violin_boxplot, K_ICP_violin_boxplot, K_corr2, K_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, K, K_sd, K_ICP, K_ICP_sd, fit, lwr, upr)) %>% 
  print()
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
K_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/K/K_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Ca_XRF_violin_boxplot, Ca_ICP_violin_boxplot, Ca_corr2, Ca_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Ca, Ca_sd, Ca_ICP, Ca_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ca_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ca/Ca_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Ti_XRF_violin_boxplot, Ti_ICP_violin_boxplot, Ti_corr2, Ti_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Ti, Ti_sd, Ti_ICP, Ti_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ti_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ti/Ti_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Mn_XRF_violin_boxplot, Mn_ICP_violin_boxplot, Mn_corr2, Mn_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Mn, Mn_sd, Mn_ICP, Mn_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Mn_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Mn/Mn_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Fe_XRF_violin_boxplot, Fe_ICP_violin_boxplot, Fe_corr2, Fe_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Fe, Fe_sd, Fe_ICP, Fe_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Fe_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Fe/Fe_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Co_XRF_violin_boxplot, Co_ICP_violin_boxplot, Co_corr2, Co_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Co, Co_sd, Co_ICP, Co_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Co_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Co/Co_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Ni_XRF_violin_boxplot, Ni_ICP_violin_boxplot, Ni_corr2, Ni_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Ni, Ni_sd, Ni_ICP, Ni_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Ni_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Ni/Ni_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Cu_XRF_violin_boxplot, Cu_ICP_violin_boxplot, Cu_corr2, Cu_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Cu, Cu_sd, Cu_ICP, Cu_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Cu_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Cu/Cu_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Zn_XRF_violin_boxplot, Zn_ICP_violin_boxplot, Zn_corr2, Zn_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Zn, Zn_sd, Zn_ICP, Zn_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Zn_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zn/Zn_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Rb_XRF_violin_boxplot, Rb_ICP_violin_boxplot, Rb_corr2, Rb_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Rb, Rb_sd, Rb_ICP, Rb_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Rb_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Rb/Rb_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Sr_XRF_violin_boxplot, Sr_ICP_violin_boxplot, Sr_corr2, Sr_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Sr, Sr_sd, Sr_ICP, Sr_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Sr_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Sr/Sr_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(Zr_XRF_violin_boxplot, Zr_ICP_violin_boxplot, Zr_corr2, Zr_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, Zr, Zr_sd, Zr_ICP, Zr_ICP_sd, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
Zr_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/Zr/Zr_Summary_stats.txt")
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Fig1_Correlation_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Fig2c_Scatterplot_all_sites.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Fig2d_Correlation_per_site.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# Figure 2: 4-part plot 
ggarrange(DM_XRF_violin_boxplot, dry_mass_pc_violin_boxplot, DM_corr2, DM_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Fig2_Correlation_Sites.pdf", 
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
  select(c(Location:midpoint, coh_inc, coh_inc_sd, dry_mass_pc, dry_mass_err, fit, lwr, upr))
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Model_predict.csv", row.names = FALSE)

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
  ggtitle(paste("ACE 95% CI & pred: ", element_title, method_title)) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=10), 
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
DM_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Fig3a_Correlation_predict.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Fig3_Correlation_Sites.pdf", 
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
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/DM/DM_Summary_stats.txt")
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


# ACE Correlation summary-------------------------------------------------------------------------
ggarrange(K_corr1, Ca_corr1, Ti_corr1, Mn_corr1, Fe_corr1, Co_corr1,
          Ni_corr1, Cu_corr1, Zn_corr1, Rb_corr1, Sr_corr1, Zr_corr1,
          nrow = 3, ncol = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_Correlation_Summary.pdf", 
       height = c(22.5), width = c(30), dpi = 600, units = "cm")
# ACE Scatterplot
plot_grid(K_corr2, Ca_corr2, Ti_corr2, Mn_corr2, Fe_corr2, Co_corr2,
          Ni_corr2, Cu_corr2, Zn_corr2, Rb_corr2, Sr_corr2, Zr_corr2,
          nrow = 3, ncol = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_Scatterplot_Summary.pdf", 
       height = c(45), width = c(60), dpi = 600, units = "cm")

# Multivariate analysis -------------------------------------------------------------------------

# itrax.R unconstrained cluster analysis - similar geochemical units without depth constrained clustering / zoning
library(compositions)


ACE_mv0 <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/cps/ACE_no_POB4.csv") %>% #import original cps data to convert to clr 
  select(Location:MSE, all_of(xrf_icp_Elements_min_PCA))
ACE_mv0

# Give each site unique depths for plotting sequentially and unique uid for next step in process
BI10_mv <- ACE_mv0 %>%
  filter(Site == "BI10") %>%
  mutate(site_depth = (midpoint/1000) + 1) %>% 
  mutate(uid = paste0(Site, "_", site_depth)) %>% 
  relocate(c("site_depth", "uid"), .after = "Section")
BI10_mv

HER42PB_mv <- ACE_mv0 %>%
  filter(Site == "HER42PB") %>%
  mutate(site_depth = (midpoint/1000) + 2) %>% 
  mutate(uid = paste0(Site, "_", site_depth)) %>% 
  relocate(c("site_depth", "uid"), .after = "Section")
HER42PB_mv

KER1_mv <- ACE_mv0 %>%
  filter(Site == "KER1") %>%
  mutate(site_depth = (midpoint/1000) + 3) %>% 
  mutate(uid = paste0(Site, "_", site_depth)) %>% 
  relocate(c("site_depth", "uid"), .after = "Section")
KER1_mv

KER3_mv <- ACE_mv0 %>%
  filter(Site == "KER3") %>%
  mutate(site_depth = (midpoint/1000) + 4) %>% 
  mutate(uid = paste0(Site, "_", site_depth)) %>% 
  relocate(c("site_depth", "uid"), .after = "Section")
KER3_mv

PB1_mv <- ACE_mv0 %>%
  filter(Site == "PB1") %>%
  mutate(site_depth = (midpoint/10000) + 5) %>% 
  mutate(uid = paste0(Site, "_", site_depth)) %>% 
  relocate(c("site_depth", "uid"), .after = "Section")
PB1_mv

# Bind into single file 
ACE_mv <- bind_rows(BI10_mv, HER42PB_mv, KER1_mv, 
                    KER3_mv, PB1_mv) %>%
  print()
tail(ACE_mv)
write.csv(ACE_mv,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_multivariate_cps.csv", row.names = FALSE)

TRUE %in% ACE_mv$uid %>% duplicated() # check nothing duplicated

# Transform closed sum data into  acomp from compositions package: "a vector of class "acomp" represents one closed 
# composition or a matrix of class "acomp" representing multiple closed compositions each in one row"

ACE_mv_acomp <- ACE_mv %>%
  #filter(qc == TRUE) %>%
  select(any_of(c(xrf_icp_Elements_min_PCA, "uid"))) %>% #
  column_to_rownames("uid") %>%
  mutate(across(everything(), function(x){ifelse(x == 0, -1, x)})) %>%
  acomp() 
head(ACE_mv_acomp)

# From itrax.R: "For these count data we deal with values that are below the limit of 
# detection by letting the compositions package know the detection limit. For these 
# count data, the detection limit is 1, and thus is coded as -1.
# This allows zeroreplace() to correctly deal with this problem."

ACE_mv_acomp_meta <- full_join(ACE_mv_acomp %>% 
                                 as.data.frame() %>%
                                 rownames_to_column("uid"),
                               ACE_mv %>%
                                 select(c(Location:MSE)), #select(-any_of(xrf_icp_Elements_key)),
                               by = "uid"
) %>%
  arrange(site_depth, Site) 


# From itrax.R: "Principal component analysis (PCA) is a common method for exploring multivariate data. Note the use of zeroreplace()
# this is because the princomp() method defined for the acomp class uses a centred-log-ratio (clr()) transformation that is intolerant to zero-values."

ACE_mv_acomp %>%
  zeroreplace() %>%
  princomp() %>%
  biplot(xlabs = rep("o",times = nrow(ACE_mv_acomp)))

# PCA - summary plot

# choose elements to plot
#MyElements <- xrf_icp_Elements_min_PCA
MyElements <- xrf_icp_Elements_min_PCA_edited

# plot grouped by site as depth not relevant
theme_set(theme_bw(base_size=12))
cluster_col5 = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF")
ACE_mv_acomp_meta %>%
  ordr::ordinate(., cols = any_of(MyElements),
                 model = ~ princomp(clr(.)),
                 augment = any_of(c("site_depth", "Site"))
  ) %>%
  
  ordr::ggbiplot(., sec.axes = "cols", scale.factor = 8) +
  ordr::geom_rows_point(aes(colour = Site, shape = Site), alpha = 1, size = 3) +
  ordr::geom_cols_vector() +
  ordr::geom_cols_text_radiate(aes(label = name)) +
  scale_colour_manual(values = cluster_col5) +
  theme(legend.text=element_text(size=12)) +
  theme(text=element_text(size=12, face = "plain"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="plain"),
        plot.margin = unit(c(1,1,1,1), "cm"))
#scale_color_continuous(type = "viridis", trans = "reverse")
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_Correlation_clr/ACE_Elements_min_PCA_key.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Unconstrained cluster analysis - by site and then depth -----------------

# Perform unconstrained cluster analysis - use broken stick value as maximum number of lith units?
ACE_mv_unclust <- left_join(ACE_mv, 
                            tibble(uid = ACE_mv_acomp %>%
                                     rownames(),
                                   group = dist(ACE_mv_acomp) %>%
                                     hclust(method = "ward.D2") %>%
                                     cutree(k = 5) %>% # use number of groups generated by 
                                     as.factor()
                            ),
                            by = "uid"
) %>%
  relocate(group, .after = "Section")

# plot against depth in each section  
ACE_mv_unclust_plot <- ggplot(ACE_mv_unclust) + 
  geom_tile(aes(x = site_depth, y = 1, fill = group)) +
  scale_y_reverse() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ACE_mv_unclust_plot
write.csv(ACE_mv_unclust,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/log_inc/ACE_multivariate_gorups.csv", row.names = FALSE)




# Linear modelling -----------------------------------------------------------

# Load required packages
library(performance) # linear model testing and graphical outputs
library(ggpmisc) # this package is used for simple/quick for labelling (stat_poly_eq)
library(lmtest) # linear model testing 
library(car) # Leverage

# Linear model titles & performance assessment /stats -----------------------------------------------

# Set up labels
#ACE_dataset <- ACE_all
ACE_dataset <- ACE_no_POB4
#ACE_dataset <-  ACE_no_POB4_PB1
#ACE_dataset <- ACE_no_PB1
site_title <- "ACE"
#itrax_dataset <- " (cps)"
#itrax_dataset <- " / inc"
#itrax_dataset <- " (%cps sum)"
itrax_dataset <- " clr"
#icp_dataset <- " [Ln]

# After first ACE_dataset runs: Set these as ACE_daatset to remove outliers according to Cook's distance results

# OLS - unweighted
#ACE_dataset <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_all_OLS_no_outliers.csv", 
#                        col_names = TRUE, skip = 0)

# OLS - weighted
#ACE_dataset <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_all_OLS_wt_no_outliers.csv", 
#                        col_names = TRUE, skip = 0)

# WLS - unweighted
#ACE_dataset <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_all_WLS_no_outliers.csv", 
#                        col_names = TRUE, skip = 0)

# WLS - weighted
#ACE_dataset <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_all_WLS_wt_no_outliers.csv", 
#                        col_names = TRUE, skip = 0)

# Convert Site as a grouping variable
ACE_dataset$Site <- as.factor(ACE_dataset$Site)

# -------------------------------------------------------------------------
# ACE_K -----------------------------------------------------------------------

# 1) OLS (Ordinary Least Squares) - linear model
ACE_K_lm <- lm(K_ICP ~ K, data = ACE_dataset)
summary(ACE_K_lm)
glance(ACE_K_lm)
model_performance(ACE_K_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_lm_hats <- as.data.frame(hatvalues(ACE_K_lm))
# Cooks distance - 2-3 x difference from mean 
ACE_K_lm_cooksD <- cooks.distance(ACE_K_lm)
ACE_K_lm_influential <- ACE_K_lm_cooksD[(ACE_K_lm_cooksD > (3 * mean(ACE_K_lm_cooksD, na.rm = TRUE)))]
ACE_K_lm_influential
ACE_K_lm_influential_names <- names(ACE_K_lm_influential)
ACE_K_lm_outliers <- ACE_dataset[ACE_K_lm_influential_names,] # outliers only using of index values
ACE_K_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_K_lm_outliers) # new dataset with outliers removed
write.csv(ACE_K_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model & checks
ACE_K_wlm <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/(K_sd)^2)
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_wlm_hats <- as.data.frame(hatvalues(ACE_K_wlm))
# Cooks distance - 2-3 x difference from mean 
ACE_K_wlm_cooksD <- cooks.distance(ACE_K_wlm)
ACE_K_wlm_influential <- ACE_K_wlm_cooksD[(ACE_K_wlm_cooksD > (3 * mean(ACE_K_wlm_cooksD, na.rm = TRUE)))]
ACE_K_wlm_influential
ACE_K_wlm_influential_names <- names(ACE_K_wlm_influential)
ACE_K_wlm_outliers <- ACE_dataset[ACE_K_wlm_influential_names,] # outliers only using of index values
ACE_K_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_K_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_K_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_OLS_wt_no_outliers.csv", row.names = FALSE)


# 3) WLS (Weighted Least Squares)
ACE_K_model <- lm(K_ICP ~ K, data = ACE_dataset) # define model
ACE_K_wt <- 1 / lm(abs(ACE_K_model$residuals) ~ ACE_K_model$fitted.values)$fitted.values^2 #define weights to use
ACE_K_wls <- lm(K_ICP ~ K, data = ACE_dataset, weights=ACE_K_wt) #perform weighted least squares regression
# Checks
summary(ACE_K_wls) # summary stats
glance(ACE_K_wls) # summary stats including AIC
model_performance(ACE_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_wls_hats <- as.data.frame(hatvalues(ACE_K_wls))
# Cooks distance - 2-3 x difference from mean 
ACE_K_wls_cooksD <- cooks.distance(ACE_K_wls)
ACE_K_wls_influential <- ACE_K_wls_cooksD[(ACE_K_wls_cooksD > (3 * mean(ACE_K_wls_cooksD, na.rm = TRUE)))]
ACE_K_wls_influential
ACE_K_wls_influential_names <- names(ACE_K_wls_influential)
ACE_K_wls_outliers <- ACE_dataset[ACE_K_wls_influential_names,] # outliers only using of index values
ACE_K_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_K_wls_outliers) # new dataset with outliers removed
write.csv(ACE_K_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_WLS_no_outliers.csv", row.names = FALSE)

# 4) WLS (ITRAX error weighted)
ACE_K_model_wt <- lm(K_ICP ~ K, data = ACE_dataset, weight = 1/K_sd^2) # define model
ACE_K_wt_wt <- 1 / lm(abs(ACE_K_model_wt$residuals) ~ ACE_K_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_K_wls_wt <- lm(K_ICP ~ K, data = ACE_dataset, weights=ACE_K_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_K_wls_wt) # summary stats
glance(ACE_K_wls_wt) # summary stats including AIC
model_performance(ACE_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_K_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_K_wls_wt_hats <- as.data.frame(hatvalues(ACE_K_wls_wt))
# Cooks distance - 2-3 x difference from mean 
ACE_K_wls_wt_cooksD <- cooks.distance(ACE_K_wls_wt)
ACE_K_wls_wt_influential <- ACE_K_wls_wt_cooksD[(ACE_K_wls_wt_cooksD > (3 * mean(ACE_K_wls_wt_cooksD, na.rm = TRUE)))]
ACE_K_wls_wt_influential
ACE_K_wls_wt_influential_names <- names(ACE_K_wls_wt_influential)
ACE_K_wls_wt_outliers <- ACE_dataset[ACE_K_wls_wt_influential_names,] # outliers only using of index values
ACE_K_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_K_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_K_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_WLS_wt_no_outliers.csv", row.names = FALSE)


# Write stats & plots to file

# 1) OLS
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_OLS_wt_summary.txt")
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wlm) # Performance package summary check for heteroscedasticity
icc(ACE_K_wlm) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_lm_lev_bar.pdf")
barplot(hatvalues(ACE_K_wlm), col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_lm_lev.pdf")
leveragePlots(ACE_K_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wlm)
dev.off()

# 2) OLS weighted
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_OLS_summary.txt")
summary(ACE_K_wlm)
glance(ACE_K_wlm)
model_performance(ACE_K_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wlm) # Performance package summary check for heteroscedasticity
icc(ACE_K_wlm) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_K_wlm), col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wlm_lev.pdf")
leveragePlots(ACE_K_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wlm)
dev.off()

# 3) WLS
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_WLS_summary.txt")
summary(ACE_K_wls)
glance(ACE_K_wls)
model_performance(ACE_K_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wls) # Performance package summary check for heteroscedasticity
icc(ACE_K_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wls_lev_bar.pdf")
barplot(hatvalues(ACE_K_wls), col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wls_lev.pdf")
leveragePlots(ACE_K_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wls)
dev.off()

# 4) WLS weighted
# Performance indicators
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_WLS_wt_summary.txt")
summary(ACE_K_wls_wt)
glance(ACE_K_wls_wt)
model_performance(ACE_K_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_K_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_K_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
# Leverage & Cooks distance - influence tests
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_K_wls_wt), col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wls_wt_lev.pdf")
leveragePlots(ACE_K_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/ACE_K_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_K_wls_wt)
dev.off()

# Linear model summary elemental plots

# K
element_title <- "K"
K_WLS_wt = ACE_K_wt
K_WLS_err_wt = ACE_K_wt_wt
theme_set(theme_classic(10))

ACE_K <- ggplot(ACE_dataset, aes(x = K, y = K_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=K-K_sd, xmax=K+K_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/K_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = K_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = K_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_K

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_K_wlm_p <- function(ACE_K_wlm) {
  f <- summary(ACE_K_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_K_wlm_p(ACE_K_wlm)

ACE_K_wlm_eqn <- function(df){
  m <- lm(K_ICP ~ K, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_K_wlm_p(ACE_K_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_K_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_K_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_K_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_K_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_K_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
K_x.reg <- ACE_dataset$K
K_y.reg <- ACE_dataset$K_ICP
# Build linear regression models
K_model1 <- lm(K_y.reg ~ K, data = ACE_dataset) # OSL unweighted
K_model2 <- lm(K_y.reg ~ K, data = ACE_dataset, weights=ACE_K_wt) # WLS unweighted
# Model summary statistics
K_model1
summary(K_model1)
confint(K_model1)
K_model2
summary(K_model2)
confint(K_model2)
# Create predicted values and upper/lower CI to check model is working
K_new1.y <- data.frame(K_x.reg = c(1, 2, 3))
predict(K_model1, K_newdata1 = K_new1.y)
predict(K_model1, K_newdata1 = K_new1.y, K_interval1 = "confidence")
K_new2.y <- data.frame(K_x.reg = c(1, 2, 3))
predict(K_model2, K_newdata2 = K_new2.y)
predict(K_model2, K_newdata2 = K_new2.y, K_interval2 = "confidence")
# Add prediction intervals to model data frame
K_pred.int1 <- predict(K_model1, interval = "prediction")
K_data_1_out <- bind_cols(ACE_dataset, K_pred.int1) %>%
  rename(K_fit_OLS = fit, K_lwr_OLS = lwr, K_upr_OLS = upr) %>% 
  select(c(Location:midpoint, K, K_sd, K_ICP, K_ICP_sd, K_fit_OLS, K_lwr_OLS, K_upr_OLS))
K_pred.int2 <- predict(K_model2, interval = "prediction")
K_data_2_out <- bind_cols(K_data_1_out, K_pred.int2) %>%
  rename(K_fit_WLS = fit, K_lwr_WLS = lwr, K_upr_WLS = upr) %>% 
  select(c(Location:midpoint, K, K_sd, K_ICP, K_ICP_sd, K_fit_OLS, K_lwr_OLS, K_upr_OLS, K_fit_WLS, K_lwr_WLS, K_upr_WLS)) %>% 
  print()
write.csv(K_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_K_predict <- ACE_K_final + 
  geom_line(data = K_data_1_out, aes(y = K_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = K_data_1_out, aes(y = K_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = K_data_2_out, aes(y = K_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = K_data_2_out, aes(y = K_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_K_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/K/K_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")



# ACE_Ca -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ca_lm <- lm(Ca_ICP ~ Ca, data = ACE_dataset)
summary(ACE_Ca_lm)
glance(ACE_Ca_lm)
model_performance(ACE_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Ca_wlm <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/(Ca_sd)^2)
summary(ACE_Ca_wlm)
glance(ACE_Ca_wlm)
model_performance(ACE_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Ca_model <- lm(Ca_ICP ~ Ca, data = ACE_dataset) # define model
ACE_Ca_wt <- 1 / lm(abs(ACE_Ca_model$residuals) ~ ACE_Ca_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ca_wls <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ca_wls) # summary stats
glance(ACE_Ca_wls) # summary stats including AIC
model_performance(ACE_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Ca_model_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weight = 1/Ca_sd^2) # define model
ACE_Ca_wt_wt <- 1 / lm(abs(ACE_Ca_model_wt$residuals) ~ ACE_Ca_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ca_wls_wt <- lm(Ca_ICP ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ca_wls_wt) # summary stats
glance(ACE_Ca_wls_wt) # summary stats including AIC
model_performance(ACE_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ca_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_lm_hats <- as.data.frame(hatvalues(ACE_Ca_lm))
ACE_Ca_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_lm_cooksD <- cooks.distance(ACE_Ca_lm)
ACE_Ca_lm_influential <- ACE_Ca_lm_cooksD[(ACE_Ca_lm_cooksD > (3 * mean(ACE_Ca_lm_cooksD, na.rm = TRUE)))]
ACE_Ca_lm_influential
ACE_Ca_lm_influential_names <- names(ACE_Ca_lm_influential)
ACE_Ca_lm_outliers <- ACE_dataset[ACE_Ca_lm_influential_names,] # outliers only using of index values
ACE_Ca_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Ca_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_wlm_hats <- as.data.frame(hatvalues(ACE_Ca_wlm))
ACE_Ca_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_wlm_cooksD <- cooks.distance(ACE_Ca_wlm)
ACE_Ca_wlm_influential <- ACE_Ca_wlm_cooksD[(ACE_Ca_wlm_cooksD > (3 * mean(ACE_Ca_wlm_cooksD, na.rm = TRUE)))]
ACE_Ca_wlm_influential
ACE_Ca_wlm_influential_names <- names(ACE_Ca_wlm_influential)
ACE_Ca_wlm_outliers <- ACE_dataset[ACE_Ca_wlm_influential_names,] # outliers only using of index values
ACE_Ca_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Ca_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_wls_hats <- as.data.frame(hatvalues(ACE_Ca_wls))
ACE_Ca_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_wls_cooksD <- cooks.distance(ACE_Ca_wls)
ACE_Ca_wls_influential <- ACE_Ca_wls_cooksD[(ACE_Ca_wls_cooksD > (3 * mean(ACE_Ca_wls_cooksD, na.rm = TRUE)))]
ACE_Ca_wls_influential
ACE_Ca_wls_influential_names <- names(ACE_Ca_wls_influential)
ACE_Ca_wls_outliers <- ACE_dataset[ACE_Ca_wls_influential_names,] # outliers only using of index values
ACE_Ca_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Ca_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ca_wls_wt_hats <- as.data.frame(hatvalues(ACE_Ca_wls_wt))
ACE_Ca_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ca_wls_wt_cooksD <- cooks.distance(ACE_Ca_wls_wt)
ACE_Ca_wls_wt_influential <- ACE_Ca_wls_wt_cooksD[(ACE_Ca_wls_wt_cooksD > (3 * mean(ACE_Ca_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Ca_wls_wt_influential
ACE_Ca_wls_wt_influential_names <- names(ACE_Ca_wls_wt_influential)
ACE_Ca_wls_wt_outliers <- ACE_dataset[ACE_Ca_wls_wt_influential_names,] # outliers only using of index values
ACE_Ca_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Ca_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Ca_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_OLS_summary.txt")
summary(ACE_Ca_lm)
glance(ACE_Ca_lm)
model_performance(ACE_Ca_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_lm_lev.pdf")
leveragePlots(ACE_Ca_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_OLS_wt_summary.txt")
summary(ACE_Ca_wlm)
glance(ACE_Ca_wlm)
model_performance(ACE_Ca_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wlm_lev.pdf")
leveragePlots(ACE_Ca_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_WLS_summary.txt")
summary(ACE_Ca_wls)
glance(ACE_Ca_wls)
model_performance(ACE_Ca_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wls_lev.pdf")
leveragePlots(ACE_Ca_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_WLS_wt_summary.txt")
summary(ACE_Ca_wls_wt)
glance(ACE_Ca_wls_wt)
model_performance(ACE_Ca_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ca_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ca_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ca_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Ca_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wls_wt_lev.pdf")
leveragePlots(ACE_Ca_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/ACE_Ca_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ca_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Ca"
Ca_WLS_wt = ACE_Ca_wt
Ca_WLS_err_wt = ACE_Ca_wt_wt
theme_set(theme_classic(10))

ACE_Ca <- ggplot(ACE_dataset, aes(x = Ca, y = Ca_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Ca_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Ca_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Ca_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Ca

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ca_wlm_p <- function(ACE_Ca_wlm) {
  f <- summary(ACE_Ca_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ca_wlm_p(ACE_Ca_wlm)

ACE_Ca_wlm_eqn <- function(df){
  m <- lm(Ca_ICP ~ Ca, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ca_wlm_p(ACE_Ca_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ca_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ca_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ca_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ca_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ca_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Ca_x.reg <- ACE_dataset$Ca
Ca_y.reg <- ACE_dataset$Ca_ICP
# Build linear regression models
Ca_model1 <- lm(Ca_y.reg ~ Ca, data = ACE_dataset) # OSL unweighted
Ca_model2 <- lm(Ca_y.reg ~ Ca, data = ACE_dataset, weights=ACE_Ca_wt) # WLS unweighted
# Model summary statistics
Ca_model1
summary(Ca_model1)
confint(Ca_model1)
Ca_model2
summary(Ca_model2)
confint(Ca_model2)
# Create predicted values and upper/lower CI to check model is working
Ca_new1.y <- data.frame(Ca_x.reg = c(1, 2, 3))
predict(Ca_model1, Ca_newdata1 = Ca_new1.y)
predict(Ca_model1, Ca_newdata1 = Ca_new1.y, Ca_interval1 = "confidence")
Ca_new2.y <- data.frame(Ca_x.reg = c(1, 2, 3))
predict(Ca_model2, Ca_newdata2 = Ca_new2.y)
predict(Ca_model2, Ca_newdata2 = Ca_new2.y, Ca_interval2 = "confidence")
# Add prediction intervals to model data frame
Ca_pred.int1 <- predict(Ca_model1, interval = "prediction")
Ca_data_1_out <- bind_cols(ACE_dataset, Ca_pred.int1) %>%
  rename(Ca_fit_OLS = fit, Ca_lwr_OLS = lwr, Ca_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Ca, Ca_sd, Ca_ICP, Ca_ICP_sd, Ca_fit_OLS, Ca_lwr_OLS, Ca_upr_OLS))
Ca_pred.int2 <- predict(Ca_model2, interval = "prediction")
Ca_data_2_out <- bind_cols(Ca_data_1_out, Ca_pred.int2) %>%
  rename(Ca_fit_WLS = fit, Ca_lwr_WLS = lwr, Ca_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Ca, Ca_sd, Ca_ICP, Ca_ICP_sd, Ca_fit_OLS, Ca_lwr_OLS, Ca_upr_OLS, Ca_fit_WLS, Ca_lwr_WLS, Ca_upr_WLS)) %>% 
  print()
write.csv(Ca_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Ca_predict <- ACE_Ca_final + 
  geom_line(data = Ca_data_1_out, aes(y = Ca_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Ca_data_1_out, aes(y = Ca_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Ca_data_2_out, aes(y = Ca_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Ca_data_2_out, aes(y = Ca_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Ca_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ca/Ca_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")


# ACE_Ti -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ti_lm <- lm(Ti_ICP ~ Ti, data = ACE_dataset)
summary(ACE_Ti_lm)
glance(ACE_Ti_lm)
model_performance(ACE_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Ti_wlm <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/(Ti_sd)^2)
summary(ACE_Ti_wlm)
glance(ACE_Ti_wlm)
model_performance(ACE_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Ti_model <- lm(Ti_ICP ~ Ti, data = ACE_dataset) # define model
ACE_Ti_wt <- 1 / lm(abs(ACE_Ti_model$residuals) ~ ACE_Ti_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ti_wls <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ti_wls) # summary stats
glance(ACE_Ti_wls) # summary stats including AIC
model_performance(ACE_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Ti_model_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weight = 1/Ti_sd^2) # define model
ACE_Ti_wt_wt <- 1 / lm(abs(ACE_Ti_model_wt$residuals) ~ ACE_Ti_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ti_wls_wt <- lm(Ti_ICP ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ti_wls_wt) # summary stats
glance(ACE_Ti_wls_wt) # summary stats including AIC
model_performance(ACE_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ti_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_lm_hats <- as.data.frame(hatvalues(ACE_Ti_lm))
ACE_Ti_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_lm_cooksD <- cooks.distance(ACE_Ti_lm)
ACE_Ti_lm_influential <- ACE_Ti_lm_cooksD[(ACE_Ti_lm_cooksD > (3 * mean(ACE_Ti_lm_cooksD, na.rm = TRUE)))]
ACE_Ti_lm_influential
ACE_Ti_lm_influential_names <- names(ACE_Ti_lm_influential)
ACE_Ti_lm_outliers <- ACE_dataset[ACE_Ti_lm_influential_names,] # outliers only using of index values
ACE_Ti_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Ti_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_wlm_hats <- as.data.frame(hatvalues(ACE_Ti_wlm))
ACE_Ti_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_wlm_cooksD <- cooks.distance(ACE_Ti_wlm)
ACE_Ti_wlm_influential <- ACE_Ti_wlm_cooksD[(ACE_Ti_wlm_cooksD > (3 * mean(ACE_Ti_wlm_cooksD, na.rm = TRUE)))]
ACE_Ti_wlm_influential
ACE_Ti_wlm_influential_names <- names(ACE_Ti_wlm_influential)
ACE_Ti_wlm_outliers <- ACE_dataset[ACE_Ti_wlm_influential_names,] # outliers only using of index values
ACE_Ti_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Ti_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_wls_hats <- as.data.frame(hatvalues(ACE_Ti_wls))
ACE_Ti_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_wls_cooksD <- cooks.distance(ACE_Ti_wls)
ACE_Ti_wls_influential <- ACE_Ti_wls_cooksD[(ACE_Ti_wls_cooksD > (3 * mean(ACE_Ti_wls_cooksD, na.rm = TRUE)))]
ACE_Ti_wls_influential
ACE_Ti_wls_influential_names <- names(ACE_Ti_wls_influential)
ACE_Ti_wls_outliers <- ACE_dataset[ACE_Ti_wls_influential_names,] # outliers only using of index values
ACE_Ti_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Ti_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ti_wls_wt_hats <- as.data.frame(hatvalues(ACE_Ti_wls_wt))
ACE_Ti_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ti_wls_wt_cooksD <- cooks.distance(ACE_Ti_wls_wt)
ACE_Ti_wls_wt_influential <- ACE_Ti_wls_wt_cooksD[(ACE_Ti_wls_wt_cooksD > (3 * mean(ACE_Ti_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Ti_wls_wt_influential
ACE_Ti_wls_wt_influential_names <- names(ACE_Ti_wls_wt_influential)
ACE_Ti_wls_wt_outliers <- ACE_dataset[ACE_Ti_wls_wt_influential_names,] # outliers only using of index values
ACE_Ti_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Ti_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Ti_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_OLS_summary.txt")
summary(ACE_Ti_lm)
glance(ACE_Ti_lm)
model_performance(ACE_Ti_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_lm_lev.pdf")
leveragePlots(ACE_Ti_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_OLS_wt_summary.txt")
summary(ACE_Ti_wlm)
glance(ACE_Ti_wlm)
model_performance(ACE_Ti_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wlm_lev.pdf")
leveragePlots(ACE_Ti_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_WLS_summary.txt")
summary(ACE_Ti_wls)
glance(ACE_Ti_wls)
model_performance(ACE_Ti_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wls_lev.pdf")
leveragePlots(ACE_Ti_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_WLS_wt_summary.txt")
summary(ACE_Ti_wls_wt)
glance(ACE_Ti_wls_wt)
model_performance(ACE_Ti_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ti_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ti_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ti_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Ti_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wls_wt_lev.pdf")
leveragePlots(ACE_Ti_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/ACE_Ti_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ti_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Ti"
Ti_WLS_wt = ACE_Ti_wt
Ti_WLS_err_wt = ACE_Ti_wt_wt
theme_set(theme_classic(10))

ACE_Ti <- ggplot(ACE_dataset, aes(x = Ti, y = Ti_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Ti-Ti_sd, xmax=Ti+Ti_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Ti_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Ti_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Ti_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Ti

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ti_wlm_p <- function(ACE_Ti_wlm) {
  f <- summary(ACE_Ti_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ti_wlm_p(ACE_Ti_wlm)

ACE_Ti_wlm_eqn <- function(df){
  m <- lm(Ti_ICP ~ Ti, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ti_wlm_p(ACE_Ti_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ti_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ti_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ti_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ti_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ti_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Ti_x.reg <- ACE_dataset$Ti
Ti_y.reg <- ACE_dataset$Ti_ICP
# Build linear regression models
Ti_model1 <- lm(Ti_y.reg ~ Ti, data = ACE_dataset) # OSL unweighted
Ti_model2 <- lm(Ti_y.reg ~ Ti, data = ACE_dataset, weights=ACE_Ti_wt) # WLS unweighted
# Model summary statistics
Ti_model1
summary(Ti_model1)
confint(Ti_model1)
Ti_model2
summary(Ti_model2)
confint(Ti_model2)
# Create predicted values and upper/lower CI to check model is working
Ti_new1.y <- data.frame(Ti_x.reg = c(1, 2, 3))
predict(Ti_model1, Ti_newdata1 = Ti_new1.y)
predict(Ti_model1, Ti_newdata1 = Ti_new1.y, Ti_interval1 = "confidence")
Ti_new2.y <- data.frame(Ti_x.reg = c(1, 2, 3))
predict(Ti_model2, Ti_newdata2 = Ti_new2.y)
predict(Ti_model2, Ti_newdata2 = Ti_new2.y, Ti_interval2 = "confidence")
# Add prediction intervals to model data frame
Ti_pred.int1 <- predict(Ti_model1, interval = "prediction")
Ti_data_1_out <- bind_cols(ACE_dataset, Ti_pred.int1) %>%
  rename(Ti_fit_OLS = fit, Ti_lwr_OLS = lwr, Ti_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Ti, Ti_sd, Ti_ICP, Ti_ICP_sd, Ti_fit_OLS, Ti_lwr_OLS, Ti_upr_OLS))
Ti_pred.int2 <- predict(Ti_model2, interval = "prediction")
Ti_data_2_out <- bind_cols(Ti_data_1_out, Ti_pred.int2) %>%
  rename(Ti_fit_WLS = fit, Ti_lwr_WLS = lwr, Ti_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Ti, Ti_sd, Ti_ICP, Ti_ICP_sd, Ti_fit_OLS, Ti_lwr_OLS, Ti_upr_OLS, Ti_fit_WLS, Ti_lwr_WLS, Ti_upr_WLS)) %>% 
  print()
write.csv(Ti_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Ti_predict <- ACE_Ti_final + 
  geom_line(data = Ti_data_1_out, aes(y = Ti_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Ti_data_1_out, aes(y = Ti_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Ti_data_2_out, aes(y = Ti_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Ti_data_2_out, aes(y = Ti_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Ti_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ti/Ti_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# ACE_Mn -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Mn_lm <- lm(Mn_ICP ~ Mn, data = ACE_dataset)
summary(ACE_Mn_lm)
glance(ACE_Mn_lm)
model_performance(ACE_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Mn_wlm <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/(Mn_sd)^2)
summary(ACE_Mn_wlm)
glance(ACE_Mn_wlm)
model_performance(ACE_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Mn_model <- lm(Mn_ICP ~ Mn, data = ACE_dataset) # define model
ACE_Mn_wt <- 1 / lm(abs(ACE_Mn_model$residuals) ~ ACE_Mn_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Mn_wls <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt) #perform weighted least squares regression
# Checks
summary(ACE_Mn_wls) # summary stats
glance(ACE_Mn_wls) # summary stats including AIC
model_performance(ACE_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Mn_model_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weight = 1/Mn_sd^2) # define model
ACE_Mn_wt_wt <- 1 / lm(abs(ACE_Mn_model_wt$residuals) ~ ACE_Mn_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Mn_wls_wt <- lm(Mn_ICP ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Mn_wls_wt) # summary stats
glance(ACE_Mn_wls_wt) # summary stats including AIC
model_performance(ACE_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Mn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_lm_hats <- as.data.frame(hatvalues(ACE_Mn_lm))
ACE_Mn_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_lm_cooksD <- cooks.distance(ACE_Mn_lm)
ACE_Mn_lm_influential <- ACE_Mn_lm_cooksD[(ACE_Mn_lm_cooksD > (3 * mean(ACE_Mn_lm_cooksD, na.rm = TRUE)))]
ACE_Mn_lm_influential
ACE_Mn_lm_influential_names <- names(ACE_Mn_lm_influential)
ACE_Mn_lm_outliers <- ACE_dataset[ACE_Mn_lm_influential_names,] # outliers only using of index values
ACE_Mn_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Mn_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_wlm_hats <- as.data.frame(hatvalues(ACE_Mn_wlm))
ACE_Mn_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_wlm_cooksD <- cooks.distance(ACE_Mn_wlm)
ACE_Mn_wlm_influential <- ACE_Mn_wlm_cooksD[(ACE_Mn_wlm_cooksD > (3 * mean(ACE_Mn_wlm_cooksD, na.rm = TRUE)))]
ACE_Mn_wlm_influential
ACE_Mn_wlm_influential_names <- names(ACE_Mn_wlm_influential)
ACE_Mn_wlm_outliers <- ACE_dataset[ACE_Mn_wlm_influential_names,] # outliers only using of index values
ACE_Mn_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Mn_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_wls_hats <- as.data.frame(hatvalues(ACE_Mn_wls))
ACE_Mn_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_wls_cooksD <- cooks.distance(ACE_Mn_wls)
ACE_Mn_wls_influential <- ACE_Mn_wls_cooksD[(ACE_Mn_wls_cooksD > (3 * mean(ACE_Mn_wls_cooksD, na.rm = TRUE)))]
ACE_Mn_wls_influential
ACE_Mn_wls_influential_names <- names(ACE_Mn_wls_influential)
ACE_Mn_wls_outliers <- ACE_dataset[ACE_Mn_wls_influential_names,] # outliers only using of index values
ACE_Mn_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Mn_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Mn_wls_wt_hats <- as.data.frame(hatvalues(ACE_Mn_wls_wt))
ACE_Mn_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Mn_wls_wt_cooksD <- cooks.distance(ACE_Mn_wls_wt)
ACE_Mn_wls_wt_influential <- ACE_Mn_wls_wt_cooksD[(ACE_Mn_wls_wt_cooksD > (3 * mean(ACE_Mn_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Mn_wls_wt_influential
ACE_Mn_wls_wt_influential_names <- names(ACE_Mn_wls_wt_influential)
ACE_Mn_wls_wt_outliers <- ACE_dataset[ACE_Mn_wls_wt_influential_names,] # outliers only using of index values
ACE_Mn_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Mn_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Mn_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_OLS_summary.txt")
summary(ACE_Mn_lm)
glance(ACE_Mn_lm)
model_performance(ACE_Mn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_lm_lev.pdf")
leveragePlots(ACE_Mn_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_OLS_wt_summary.txt")
summary(ACE_Mn_wlm)
glance(ACE_Mn_wlm)
model_performance(ACE_Mn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wlm_lev.pdf")
leveragePlots(ACE_Mn_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_WLS_summary.txt")
summary(ACE_Mn_wls)
glance(ACE_Mn_wls)
model_performance(ACE_Mn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wls_lev.pdf")
leveragePlots(ACE_Mn_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_WLS_wt_summary.txt")
summary(ACE_Mn_wls_wt)
glance(ACE_Mn_wls_wt)
model_performance(ACE_Mn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Mn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Mn_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Mn_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Mn_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wls_wt_lev.pdf")
leveragePlots(ACE_Mn_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/ACE_Mn_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Mn_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Mn"
Mn_WLS_wt = ACE_Mn_wt
Mn_WLS_err_wt = ACE_Mn_wt_wt
theme_set(theme_classic(10))

ACE_Mn <- ggplot(ACE_dataset, aes(x = Mn, y = Mn_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Mn-Mn_sd, xmax=Mn+Mn_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Mn_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Mn_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Mn_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Mn

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Mn_wlm_p <- function(ACE_Mn_wlm) {
  f <- summary(ACE_Mn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Mn_wlm_p(ACE_Mn_wlm)

ACE_Mn_wlm_eqn <- function(df){
  m <- lm(Mn_ICP ~ Mn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Mn_wlm_p(ACE_Mn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Mn_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Mn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Mn_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Mn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Mn_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Mn_x.reg <- ACE_dataset$Mn
Mn_y.reg <- ACE_dataset$Mn_ICP
# Build linear regression models
Mn_model1 <- lm(Mn_y.reg ~ Mn, data = ACE_dataset) # OSL unweighted
Mn_model2 <- lm(Mn_y.reg ~ Mn, data = ACE_dataset, weights=ACE_Mn_wt) # WLS unweighted
# Model summary statistics
Mn_model1
summary(Mn_model1)
confint(Mn_model1)
Mn_model2
summary(Mn_model2)
confint(Mn_model2)
# Create predicted values and upper/lower CI to check model is working
Mn_new1.y <- data.frame(Mn_x.reg = c(1, 2, 3))
predict(Mn_model1, Mn_newdata1 = Mn_new1.y)
predict(Mn_model1, Mn_newdata1 = Mn_new1.y, Mn_interval1 = "confidence")
Mn_new2.y <- data.frame(Mn_x.reg = c(1, 2, 3))
predict(Mn_model2, Mn_newdata2 = Mn_new2.y)
predict(Mn_model2, Mn_newdata2 = Mn_new2.y, Mn_interval2 = "confidence")
# Add prediction intervals to model data frame
Mn_pred.int1 <- predict(Mn_model1, interval = "prediction")
Mn_data_1_out <- bind_cols(ACE_dataset, Mn_pred.int1) %>%
  rename(Mn_fit_OLS = fit, Mn_lwr_OLS = lwr, Mn_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Mn, Mn_sd, Mn_ICP, Mn_ICP_sd, Mn_fit_OLS, Mn_lwr_OLS, Mn_upr_OLS))
Mn_pred.int2 <- predict(Mn_model2, interval = "prediction")
Mn_data_2_out <- bind_cols(Mn_data_1_out, Mn_pred.int2) %>%
  rename(Mn_fit_WLS = fit, Mn_lwr_WLS = lwr, Mn_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Mn, Mn_sd, Mn_ICP, Mn_ICP_sd, Mn_fit_OLS, Mn_lwr_OLS, Mn_upr_OLS, Mn_fit_WLS, Mn_lwr_WLS, Mn_upr_WLS)) %>% 
  print()
write.csv(Mn_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Mn_predict <- ACE_Mn_final + 
  geom_line(data = Mn_data_1_out, aes(y = Mn_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Mn_data_1_out, aes(y = Mn_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Mn_data_2_out, aes(y = Mn_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Mn_data_2_out, aes(y = Mn_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Mn_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Mn/Mn_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")


# ACE_Fe -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Fe_lm <- lm(Fe_ICP ~ Fe, data = ACE_dataset)
summary(ACE_Fe_lm)
glance(ACE_Fe_lm)
model_performance(ACE_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Fe_wlm <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/(Fe_sd)^2)
summary(ACE_Fe_wlm)
glance(ACE_Fe_wlm)
model_performance(ACE_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Fe_model <- lm(Fe_ICP ~ Fe, data = ACE_dataset) # define model
ACE_Fe_wt <- 1 / lm(abs(ACE_Fe_model$residuals) ~ ACE_Fe_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Fe_wls <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt) #perform weighted least squares regression
# Checks
summary(ACE_Fe_wls) # summary stats
glance(ACE_Fe_wls) # summary stats including AIC
model_performance(ACE_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Fe_model_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weight = 1/Fe_sd^2) # define model
ACE_Fe_wt_wt <- 1 / lm(abs(ACE_Fe_model_wt$residuals) ~ ACE_Fe_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Fe_wls_wt <- lm(Fe_ICP ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Fe_wls_wt) # summary stats
glance(ACE_Fe_wls_wt) # summary stats including AIC
model_performance(ACE_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Fe_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_lm_hats <- as.data.frame(hatvalues(ACE_Fe_lm))
ACE_Fe_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_lm_cooksD <- cooks.distance(ACE_Fe_lm)
ACE_Fe_lm_influential <- ACE_Fe_lm_cooksD[(ACE_Fe_lm_cooksD > (3 * mean(ACE_Fe_lm_cooksD, na.rm = TRUE)))]
ACE_Fe_lm_influential
ACE_Fe_lm_influential_names <- names(ACE_Fe_lm_influential)
ACE_Fe_lm_outliers <- ACE_dataset[ACE_Fe_lm_influential_names,] # outliers only using of index values
ACE_Fe_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Fe_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_wlm_hats <- as.data.frame(hatvalues(ACE_Fe_wlm))
ACE_Fe_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_wlm_cooksD <- cooks.distance(ACE_Fe_wlm)
ACE_Fe_wlm_influential <- ACE_Fe_wlm_cooksD[(ACE_Fe_wlm_cooksD > (3 * mean(ACE_Fe_wlm_cooksD, na.rm = TRUE)))]
ACE_Fe_wlm_influential
ACE_Fe_wlm_influential_names <- names(ACE_Fe_wlm_influential)
ACE_Fe_wlm_outliers <- ACE_dataset[ACE_Fe_wlm_influential_names,] # outliers only using of index values
ACE_Fe_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Fe_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_wls_hats <- as.data.frame(hatvalues(ACE_Fe_wls))
ACE_Fe_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_wls_cooksD <- cooks.distance(ACE_Fe_wls)
ACE_Fe_wls_influential <- ACE_Fe_wls_cooksD[(ACE_Fe_wls_cooksD > (3 * mean(ACE_Fe_wls_cooksD, na.rm = TRUE)))]
ACE_Fe_wls_influential
ACE_Fe_wls_influential_names <- names(ACE_Fe_wls_influential)
ACE_Fe_wls_outliers <- ACE_dataset[ACE_Fe_wls_influential_names,] # outliers only using of index values
ACE_Fe_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Fe_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Fe_wls_wt_hats <- as.data.frame(hatvalues(ACE_Fe_wls_wt))
ACE_Fe_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Fe_wls_wt_cooksD <- cooks.distance(ACE_Fe_wls_wt)
ACE_Fe_wls_wt_influential <- ACE_Fe_wls_wt_cooksD[(ACE_Fe_wls_wt_cooksD > (3 * mean(ACE_Fe_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Fe_wls_wt_influential
ACE_Fe_wls_wt_influential_names <- names(ACE_Fe_wls_wt_influential)
ACE_Fe_wls_wt_outliers <- ACE_dataset[ACE_Fe_wls_wt_influential_names,] # outliers only using of index values
ACE_Fe_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Fe_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Fe_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_OLS_summary.txt")
summary(ACE_Fe_lm)
glance(ACE_Fe_lm)
model_performance(ACE_Fe_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_lm_lev.pdf")
leveragePlots(ACE_Fe_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_OLS_wt_summary.txt")
summary(ACE_Fe_wlm)
glance(ACE_Fe_wlm)
model_performance(ACE_Fe_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wlm_lev.pdf")
leveragePlots(ACE_Fe_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_WLS_summary.txt")
summary(ACE_Fe_wls)
glance(ACE_Fe_wls)
model_performance(ACE_Fe_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wls_lev.pdf")
leveragePlots(ACE_Fe_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_WLS_wt_summary.txt")
summary(ACE_Fe_wls_wt)
glance(ACE_Fe_wls_wt)
model_performance(ACE_Fe_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Fe_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Fe_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Fe_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Fe_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wls_wt_lev.pdf")
leveragePlots(ACE_Fe_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/ACE_Fe_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Fe_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Fe"
Fe_WLS_wt = ACE_Fe_wt
Fe_WLS_err_wt = ACE_Fe_wt_wt
theme_set(theme_classic(10))

ACE_Fe <- ggplot(ACE_dataset, aes(x = Fe, y = Fe_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Fe-Fe_sd, xmax=Fe+Fe_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Fe_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Fe_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Fe_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Fe

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Fe_wlm_p <- function(ACE_Fe_wlm) {
  f <- summary(ACE_Fe_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Fe_wlm_p(ACE_Fe_wlm)

ACE_Fe_wlm_eqn <- function(df){
  m <- lm(Fe_ICP ~ Fe, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Fe_wlm_p(ACE_Fe_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Fe_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Fe_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Fe_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Fe_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Fe_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Fe_x.reg <- ACE_dataset$Fe
Fe_y.reg <- ACE_dataset$Fe_ICP
# Build linear regression models
Fe_model1 <- lm(Fe_y.reg ~ Fe, data = ACE_dataset) # OSL unweighted
Fe_model2 <- lm(Fe_y.reg ~ Fe, data = ACE_dataset, weights=ACE_Fe_wt) # WLS unweighted
# Model summary statistics
Fe_model1
summary(Fe_model1)
confint(Fe_model1)
Fe_model2
summary(Fe_model2)
confint(Fe_model2)
# Create predicted values and upper/lower CI to check model is working
Fe_new1.y <- data.frame(Fe_x.reg = c(1, 2, 3))
predict(Fe_model1, Fe_newdata1 = Fe_new1.y)
predict(Fe_model1, Fe_newdata1 = Fe_new1.y, Fe_interval1 = "confidence")
Fe_new2.y <- data.frame(Fe_x.reg = c(1, 2, 3))
predict(Fe_model2, Fe_newdata2 = Fe_new2.y)
predict(Fe_model2, Fe_newdata2 = Fe_new2.y, Fe_interval2 = "confidence")
# Add prediction intervals to model data frame
Fe_pred.int1 <- predict(Fe_model1, interval = "prediction")
Fe_data_1_out <- bind_cols(ACE_dataset, Fe_pred.int1) %>%
  rename(Fe_fit_OLS = fit, Fe_lwr_OLS = lwr, Fe_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Fe, Fe_sd, Fe_ICP, Fe_ICP_sd, Fe_fit_OLS, Fe_lwr_OLS, Fe_upr_OLS))
Fe_pred.int2 <- predict(Fe_model2, interval = "prediction")
Fe_data_2_out <- bind_cols(Fe_data_1_out, Fe_pred.int2) %>%
  rename(Fe_fit_WLS = fit, Fe_lwr_WLS = lwr, Fe_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Fe, Fe_sd, Fe_ICP, Fe_ICP_sd, Fe_fit_OLS, Fe_lwr_OLS, Fe_upr_OLS, Fe_fit_WLS, Fe_lwr_WLS, Fe_upr_WLS)) %>% 
  print()
write.csv(Fe_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Fe_predict <- ACE_Fe_final + 
  geom_line(data = Fe_data_1_out, aes(y = Fe_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Fe_data_1_out, aes(y = Fe_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Fe_data_2_out, aes(y = Fe_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Fe_data_2_out, aes(y = Fe_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Fe_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Fe/Fe_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Co -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Co_lm <- lm(Co_ICP ~ Co, data = ACE_dataset)
summary(ACE_Co_lm)
glance(ACE_Co_lm)
model_performance(ACE_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Co_wlm <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/(Co_sd)^2)
summary(ACE_Co_wlm)
glance(ACE_Co_wlm)
model_performance(ACE_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Co_model <- lm(Co_ICP ~ Co, data = ACE_dataset) # define model
ACE_Co_wt <- 1 / lm(abs(ACE_Co_model$residuals) ~ ACE_Co_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Co_wls <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=ACE_Co_wt) #perform weighted least squares regression
# Checks
summary(ACE_Co_wls) # summary stats
glance(ACE_Co_wls) # summary stats including AIC
model_performance(ACE_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Co_model_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weight = 1/Co_sd^2) # define model
ACE_Co_wt_wt <- 1 / lm(abs(ACE_Co_model_wt$residuals) ~ ACE_Co_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Co_wls_wt <- lm(Co_ICP ~ Co, data = ACE_dataset, weights=ACE_Co_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Co_wls_wt) # summary stats
glance(ACE_Co_wls_wt) # summary stats including AIC
model_performance(ACE_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Co_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_lm_hats <- as.data.frame(hatvalues(ACE_Co_lm))
ACE_Co_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_lm_cooksD <- cooks.distance(ACE_Co_lm)
ACE_Co_lm_influential <- ACE_Co_lm_cooksD[(ACE_Co_lm_cooksD > (3 * mean(ACE_Co_lm_cooksD, na.rm = TRUE)))]
ACE_Co_lm_influential
ACE_Co_lm_influential_names <- names(ACE_Co_lm_influential)
ACE_Co_lm_outliers <- ACE_dataset[ACE_Co_lm_influential_names,] # outliers only using of index values
ACE_Co_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Co_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_wlm_hats <- as.data.frame(hatvalues(ACE_Co_wlm))
ACE_Co_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_wlm_cooksD <- cooks.distance(ACE_Co_wlm)
ACE_Co_wlm_influential <- ACE_Co_wlm_cooksD[(ACE_Co_wlm_cooksD > (3 * mean(ACE_Co_wlm_cooksD, na.rm = TRUE)))]
ACE_Co_wlm_influential
ACE_Co_wlm_influential_names <- names(ACE_Co_wlm_influential)
ACE_Co_wlm_outliers <- ACE_dataset[ACE_Co_wlm_influential_names,] # outliers only using of index values
ACE_Co_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Co_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_wls_hats <- as.data.frame(hatvalues(ACE_Co_wls))
ACE_Co_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_wls_cooksD <- cooks.distance(ACE_Co_wls)
ACE_Co_wls_influential <- ACE_Co_wls_cooksD[(ACE_Co_wls_cooksD > (3 * mean(ACE_Co_wls_cooksD, na.rm = TRUE)))]
ACE_Co_wls_influential
ACE_Co_wls_influential_names <- names(ACE_Co_wls_influential)
ACE_Co_wls_outliers <- ACE_dataset[ACE_Co_wls_influential_names,] # outliers only using of index values
ACE_Co_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Co_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Co_wls_wt_hats <- as.data.frame(hatvalues(ACE_Co_wls_wt))
ACE_Co_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Co_wls_wt_cooksD <- cooks.distance(ACE_Co_wls_wt)
ACE_Co_wls_wt_influential <- ACE_Co_wls_wt_cooksD[(ACE_Co_wls_wt_cooksD > (3 * mean(ACE_Co_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Co_wls_wt_influential
ACE_Co_wls_wt_influential_names <- names(ACE_Co_wls_wt_influential)
ACE_Co_wls_wt_outliers <- ACE_dataset[ACE_Co_wls_wt_influential_names,] # outliers only using of index values
ACE_Co_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Co_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Co_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_OLS_summary.txt")
summary(ACE_Co_lm)
glance(ACE_Co_lm)
model_performance(ACE_Co_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Co_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_lm_lev.pdf")
leveragePlots(ACE_Co_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_OLS_wt_summary.txt")
summary(ACE_Co_wlm)
glance(ACE_Co_wlm)
model_performance(ACE_Co_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Co_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Co_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wlm_lev.pdf")
leveragePlots(ACE_Co_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_WLS_summary.txt")
summary(ACE_Co_wls)
glance(ACE_Co_wls)
model_performance(ACE_Co_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Co_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Co_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wls_lev.pdf")
leveragePlots(ACE_Co_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_WLS_wt_summary.txt")
summary(ACE_Co_wls_wt)
glance(ACE_Co_wls_wt)
model_performance(ACE_Co_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Co_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Co_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Co_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Co_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wls_wt_lev.pdf")
leveragePlots(ACE_Co_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/ACE_Co_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Co_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Co"
Co_WLS_wt = ACE_Co_wt
Co_WLS_err_wt = ACE_Co_wt_wt
theme_set(theme_classic(10))

ACE_Co <- ggplot(ACE_dataset, aes(x = Co, y = Co_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Co-Co_sd, xmax=Co+Co_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Co_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Co_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Co_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Co

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Co_wlm_p <- function(ACE_Co_wlm) {
  f <- summary(ACE_Co_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Co_wlm_p(ACE_Co_wlm)

ACE_Co_wlm_eqn <- function(df){
  m <- lm(Co_ICP ~ Co, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Co_wlm_p(ACE_Co_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Co_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Co_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Co_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Co_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Co_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Co_x.reg <- ACE_dataset$Co
Co_y.reg <- ACE_dataset$Co_ICP
# Build linear regression models
Co_model1 <- lm(Co_y.reg ~ Co, data = ACE_dataset) # OSL unweighted
Co_model2 <- lm(Co_y.reg ~ Co, data = ACE_dataset, weights=ACE_Co_wt) # WLS unweighted
# Model summary statistics
Co_model1
summary(Co_model1)
confint(Co_model1)
Co_model2
summary(Co_model2)
confint(Co_model2)
# Create predicted values and upper/lower CI to check model is working
Co_new1.y <- data.frame(Co_x.reg = c(1, 2, 3))
predict(Co_model1, Co_newdata1 = Co_new1.y)
predict(Co_model1, Co_newdata1 = Co_new1.y, Co_interval1 = "confidence")
Co_new2.y <- data.frame(Co_x.reg = c(1, 2, 3))
predict(Co_model2, Co_newdata2 = Co_new2.y)
predict(Co_model2, Co_newdata2 = Co_new2.y, Co_interval2 = "confidence")
# Add prediction intervals to model data frame
Co_pred.int1 <- predict(Co_model1, interval = "prediction")
Co_data_1_out <- bind_cols(ACE_dataset, Co_pred.int1) %>%
  rename(Co_fit_OLS = fit, Co_lwr_OLS = lwr, Co_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Co, Co_sd, Co_ICP, Co_ICP_sd, Co_fit_OLS, Co_lwr_OLS, Co_upr_OLS))
Co_pred.int2 <- predict(Co_model2, interval = "prediction")
Co_data_2_out <- bind_cols(Co_data_1_out, Co_pred.int2) %>%
  rename(Co_fit_WLS = fit, Co_lwr_WLS = lwr, Co_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Co, Co_sd, Co_ICP, Co_ICP_sd, Co_fit_OLS, Co_lwr_OLS, Co_upr_OLS, Co_fit_WLS, Co_lwr_WLS, Co_upr_WLS)) %>% 
  print()
write.csv(Co_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Co_predict <- ACE_Co_final + 
  geom_line(data = Co_data_1_out, aes(y = Co_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Co_data_1_out, aes(y = Co_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Co_data_2_out, aes(y = Co_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Co_data_2_out, aes(y = Co_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Co_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Co/Co_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Ni -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Ni_lm <- lm(Ni_ICP ~ Ni, data = ACE_dataset)
summary(ACE_Ni_lm)
glance(ACE_Ni_lm)
model_performance(ACE_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Ni_wlm <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/(Ni_sd)^2)
summary(ACE_Ni_wlm)
glance(ACE_Ni_wlm)
model_performance(ACE_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Ni_model <- lm(Ni_ICP ~ Ni, data = ACE_dataset) # define model
ACE_Ni_wt <- 1 / lm(abs(ACE_Ni_model$residuals) ~ ACE_Ni_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Ni_wls <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ni_wls) # summary stats
glance(ACE_Ni_wls) # summary stats including AIC
model_performance(ACE_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Ni_model_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weight = 1/Ni_sd^2) # define model
ACE_Ni_wt_wt <- 1 / lm(abs(ACE_Ni_model_wt$residuals) ~ ACE_Ni_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Ni_wls_wt <- lm(Ni_ICP ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Ni_wls_wt) # summary stats
glance(ACE_Ni_wls_wt) # summary stats including AIC
model_performance(ACE_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Ni_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_lm_hats <- as.data.frame(hatvalues(ACE_Ni_lm))
ACE_Ni_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_lm_cooksD <- cooks.distance(ACE_Ni_lm)
ACE_Ni_lm_influential <- ACE_Ni_lm_cooksD[(ACE_Ni_lm_cooksD > (3 * mean(ACE_Ni_lm_cooksD, na.rm = TRUE)))]
ACE_Ni_lm_influential
ACE_Ni_lm_influential_names <- names(ACE_Ni_lm_influential)
ACE_Ni_lm_outliers <- ACE_dataset[ACE_Ni_lm_influential_names,] # outliers only using of index values
ACE_Ni_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Ni_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_wlm_hats <- as.data.frame(hatvalues(ACE_Ni_wlm))
ACE_Ni_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_wlm_cooksD <- cooks.distance(ACE_Ni_wlm)
ACE_Ni_wlm_influential <- ACE_Ni_wlm_cooksD[(ACE_Ni_wlm_cooksD > (3 * mean(ACE_Ni_wlm_cooksD, na.rm = TRUE)))]
ACE_Ni_wlm_influential
ACE_Ni_wlm_influential_names <- names(ACE_Ni_wlm_influential)
ACE_Ni_wlm_outliers <- ACE_dataset[ACE_Ni_wlm_influential_names,] # outliers only using of index values
ACE_Ni_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Ni_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_wls_hats <- as.data.frame(hatvalues(ACE_Ni_wls))
ACE_Ni_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_wls_cooksD <- cooks.distance(ACE_Ni_wls)
ACE_Ni_wls_influential <- ACE_Ni_wls_cooksD[(ACE_Ni_wls_cooksD > (3 * mean(ACE_Ni_wls_cooksD, na.rm = TRUE)))]
ACE_Ni_wls_influential
ACE_Ni_wls_influential_names <- names(ACE_Ni_wls_influential)
ACE_Ni_wls_outliers <- ACE_dataset[ACE_Ni_wls_influential_names,] # outliers only using of index values
ACE_Ni_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Ni_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Ni_wls_wt_hats <- as.data.frame(hatvalues(ACE_Ni_wls_wt))
ACE_Ni_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Ni_wls_wt_cooksD <- cooks.distance(ACE_Ni_wls_wt)
ACE_Ni_wls_wt_influential <- ACE_Ni_wls_wt_cooksD[(ACE_Ni_wls_wt_cooksD > (3 * mean(ACE_Ni_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Ni_wls_wt_influential
ACE_Ni_wls_wt_influential_names <- names(ACE_Ni_wls_wt_influential)
ACE_Ni_wls_wt_outliers <- ACE_dataset[ACE_Ni_wls_wt_influential_names,] # outliers only using of index values
ACE_Ni_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Ni_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Ni_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_OLS_summary.txt")
summary(ACE_Ni_lm)
glance(ACE_Ni_lm)
model_performance(ACE_Ni_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_lm_lev.pdf")
leveragePlots(ACE_Ni_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_OLS_wt_summary.txt")
summary(ACE_Ni_wlm)
glance(ACE_Ni_wlm)
model_performance(ACE_Ni_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wlm_lev.pdf")
leveragePlots(ACE_Ni_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_WLS_summary.txt")
summary(ACE_Ni_wls)
glance(ACE_Ni_wls)
model_performance(ACE_Ni_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wls_lev.pdf")
leveragePlots(ACE_Ni_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_WLS_wt_summary.txt")
summary(ACE_Ni_wls_wt)
glance(ACE_Ni_wls_wt)
model_performance(ACE_Ni_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Ni_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Ni_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Ni_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Ni_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wls_wt_lev.pdf")
leveragePlots(ACE_Ni_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/ACE_Ni_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Ni_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Ni"
Ni_WLS_wt = ACE_Ni_wt
Ni_WLS_err_wt = ACE_Ni_wt_wt
theme_set(theme_classic(10))

ACE_Ni <- ggplot(ACE_dataset, aes(x = Ni, y = Ni_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Ni-Ni_sd, xmax=Ni+Ni_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Ni_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Ni_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Ni_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Ni

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Ni_wlm_p <- function(ACE_Ni_wlm) {
  f <- summary(ACE_Ni_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Ni_wlm_p(ACE_Ni_wlm)

ACE_Ni_wlm_eqn <- function(df){
  m <- lm(Ni_ICP ~ Ni, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Ni_wlm_p(ACE_Ni_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Ni_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Ni_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Ni_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Ni_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Ni_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Ni_x.reg <- ACE_dataset$Ni
Ni_y.reg <- ACE_dataset$Ni_ICP
# Build linear regression models
Ni_model1 <- lm(Ni_y.reg ~ Ni, data = ACE_dataset) # OSL unweighted
Ni_model2 <- lm(Ni_y.reg ~ Ni, data = ACE_dataset, weights=ACE_Ni_wt) # WLS unweighted
# Model summary statistics
Ni_model1
summary(Ni_model1)
confint(Ni_model1)
Ni_model2
summary(Ni_model2)
confint(Ni_model2)
# Create predicted values and upper/lower CI to check model is working
Ni_new1.y <- data.frame(Ni_x.reg = c(1, 2, 3))
predict(Ni_model1, Ni_newdata1 = Ni_new1.y)
predict(Ni_model1, Ni_newdata1 = Ni_new1.y, Ni_interval1 = "confidence")
Ni_new2.y <- data.frame(Ni_x.reg = c(1, 2, 3))
predict(Ni_model2, Ni_newdata2 = Ni_new2.y)
predict(Ni_model2, Ni_newdata2 = Ni_new2.y, Ni_interval2 = "confidence")
# Add prediction intervals to model data frame
Ni_pred.int1 <- predict(Ni_model1, interval = "prediction")
Ni_data_1_out <- bind_cols(ACE_dataset, Ni_pred.int1) %>%
  rename(Ni_fit_OLS = fit, Ni_lwr_OLS = lwr, Ni_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Ni, Ni_sd, Ni_ICP, Ni_ICP_sd, Ni_fit_OLS, Ni_lwr_OLS, Ni_upr_OLS))
Ni_pred.int2 <- predict(Ni_model2, interval = "prediction")
Ni_data_2_out <- bind_cols(Ni_data_1_out, Ni_pred.int2) %>%
  rename(Ni_fit_WLS = fit, Ni_lwr_WLS = lwr, Ni_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Ni, Ni_sd, Ni_ICP, Ni_ICP_sd, Ni_fit_OLS, Ni_lwr_OLS, Ni_upr_OLS, Ni_fit_WLS, Ni_lwr_WLS, Ni_upr_WLS)) %>% 
  print()
write.csv(Ni_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Ni_predict <- ACE_Ni_final + 
  geom_line(data = Ni_data_1_out, aes(y = Ni_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Ni_data_1_out, aes(y = Ni_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Ni_data_2_out, aes(y = Ni_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Ni_data_2_out, aes(y = Ni_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Ni_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Ni/Ni_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Cu -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Cu_lm <- lm(Cu_ICP ~ Cu, data = ACE_dataset)
summary(ACE_Cu_lm)
glance(ACE_Cu_lm)
model_performance(ACE_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Cu_wlm <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/(Cu_sd)^2)
summary(ACE_Cu_wlm)
glance(ACE_Cu_wlm)
model_performance(ACE_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Cu_model <- lm(Cu_ICP ~ Cu, data = ACE_dataset) # define model
ACE_Cu_wt <- 1 / lm(abs(ACE_Cu_model$residuals) ~ ACE_Cu_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Cu_wls <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt) #perform weighted least squares regression
# Checks
summary(ACE_Cu_wls) # summary stats
glance(ACE_Cu_wls) # summary stats including AIC
model_performance(ACE_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Cu_model_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weight = 1/Cu_sd^2) # define model
ACE_Cu_wt_wt <- 1 / lm(abs(ACE_Cu_model_wt$residuals) ~ ACE_Cu_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Cu_wls_wt <- lm(Cu_ICP ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Cu_wls_wt) # summary stats
glance(ACE_Cu_wls_wt) # summary stats including AIC
model_performance(ACE_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Cu_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_lm_hats <- as.data.frame(hatvalues(ACE_Cu_lm))
ACE_Cu_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_lm_cooksD <- cooks.distance(ACE_Cu_lm)
ACE_Cu_lm_influential <- ACE_Cu_lm_cooksD[(ACE_Cu_lm_cooksD > (3 * mean(ACE_Cu_lm_cooksD, na.rm = TRUE)))]
ACE_Cu_lm_influential
ACE_Cu_lm_influential_names <- names(ACE_Cu_lm_influential)
ACE_Cu_lm_outliers <- ACE_dataset[ACE_Cu_lm_influential_names,] # outliers only using of index values
ACE_Cu_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Cu_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_wlm_hats <- as.data.frame(hatvalues(ACE_Cu_wlm))
ACE_Cu_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_wlm_cooksD <- cooks.distance(ACE_Cu_wlm)
ACE_Cu_wlm_influential <- ACE_Cu_wlm_cooksD[(ACE_Cu_wlm_cooksD > (3 * mean(ACE_Cu_wlm_cooksD, na.rm = TRUE)))]
ACE_Cu_wlm_influential
ACE_Cu_wlm_influential_names <- names(ACE_Cu_wlm_influential)
ACE_Cu_wlm_outliers <- ACE_dataset[ACE_Cu_wlm_influential_names,] # outliers only using of index values
ACE_Cu_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Cu_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_wls_hats <- as.data.frame(hatvalues(ACE_Cu_wls))
ACE_Cu_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_wls_cooksD <- cooks.distance(ACE_Cu_wls)
ACE_Cu_wls_influential <- ACE_Cu_wls_cooksD[(ACE_Cu_wls_cooksD > (3 * mean(ACE_Cu_wls_cooksD, na.rm = TRUE)))]
ACE_Cu_wls_influential
ACE_Cu_wls_influential_names <- names(ACE_Cu_wls_influential)
ACE_Cu_wls_outliers <- ACE_dataset[ACE_Cu_wls_influential_names,] # outliers only using of index values
ACE_Cu_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Cu_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Cu_wls_wt_hats <- as.data.frame(hatvalues(ACE_Cu_wls_wt))
ACE_Cu_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Cu_wls_wt_cooksD <- cooks.distance(ACE_Cu_wls_wt)
ACE_Cu_wls_wt_influential <- ACE_Cu_wls_wt_cooksD[(ACE_Cu_wls_wt_cooksD > (3 * mean(ACE_Cu_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Cu_wls_wt_influential
ACE_Cu_wls_wt_influential_names <- names(ACE_Cu_wls_wt_influential)
ACE_Cu_wls_wt_outliers <- ACE_dataset[ACE_Cu_wls_wt_influential_names,] # outliers only using of index values
ACE_Cu_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Cu_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Cu_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_OLS_summary.txt")
summary(ACE_Cu_lm)
glance(ACE_Cu_lm)
model_performance(ACE_Cu_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_lm_lev.pdf")
leveragePlots(ACE_Cu_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_OLS_wt_summary.txt")
summary(ACE_Cu_wlm)
glance(ACE_Cu_wlm)
model_performance(ACE_Cu_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wlm_lev.pdf")
leveragePlots(ACE_Cu_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_WLS_summary.txt")
summary(ACE_Cu_wls)
glance(ACE_Cu_wls)
model_performance(ACE_Cu_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wls_lev.pdf")
leveragePlots(ACE_Cu_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_WLS_wt_summary.txt")
summary(ACE_Cu_wls_wt)
glance(ACE_Cu_wls_wt)
model_performance(ACE_Cu_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Cu_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Cu_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Cu_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Cu_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wls_wt_lev.pdf")
leveragePlots(ACE_Cu_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/ACE_Cu_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Cu_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Cu"
Cu_WLS_wt = ACE_Cu_wt
Cu_WLS_err_wt = ACE_Cu_wt_wt
theme_set(theme_classic(10))

ACE_Cu <- ggplot(ACE_dataset, aes(x = Cu, y = Cu_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Cu-Cu_sd, xmax=Cu+Cu_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Cu_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Cu_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Cu_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Cu

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Cu_wlm_p <- function(ACE_Cu_wlm) {
  f <- summary(ACE_Cu_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Cu_wlm_p(ACE_Cu_wlm)

ACE_Cu_wlm_eqn <- function(df){
  m <- lm(Cu_ICP ~ Cu, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Cu_wlm_p(ACE_Cu_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Cu_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Cu_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Cu_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Cu_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Cu_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Cu_x.reg <- ACE_dataset$Cu
Cu_y.reg <- ACE_dataset$Cu_ICP
# Build linear regression models
Cu_model1 <- lm(Cu_y.reg ~ Cu, data = ACE_dataset) # OSL unweighted
Cu_model2 <- lm(Cu_y.reg ~ Cu, data = ACE_dataset, weights=ACE_Cu_wt) # WLS unweighted
# Model summary statistics
Cu_model1
summary(Cu_model1)
confint(Cu_model1)
Cu_model2
summary(Cu_model2)
confint(Cu_model2)
# Create predicted values and upper/lower CI to check model is working
Cu_new1.y <- data.frame(Cu_x.reg = c(1, 2, 3))
predict(Cu_model1, Cu_newdata1 = Cu_new1.y)
predict(Cu_model1, Cu_newdata1 = Cu_new1.y, Cu_interval1 = "confidence")
Cu_new2.y <- data.frame(Cu_x.reg = c(1, 2, 3))
predict(Cu_model2, Cu_newdata2 = Cu_new2.y)
predict(Cu_model2, Cu_newdata2 = Cu_new2.y, Cu_interval2 = "confidence")
# Add prediction intervals to model data frame
Cu_pred.int1 <- predict(Cu_model1, interval = "prediction")
Cu_data_1_out <- bind_cols(ACE_dataset, Cu_pred.int1) %>%
  rename(Cu_fit_OLS = fit, Cu_lwr_OLS = lwr, Cu_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Cu, Cu_sd, Cu_ICP, Cu_ICP_sd, Cu_fit_OLS, Cu_lwr_OLS, Cu_upr_OLS))
Cu_pred.int2 <- predict(Cu_model2, interval = "prediction")
Cu_data_2_out <- bind_cols(Cu_data_1_out, Cu_pred.int2) %>%
  rename(Cu_fit_WLS = fit, Cu_lwr_WLS = lwr, Cu_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Cu, Cu_sd, Cu_ICP, Cu_ICP_sd, Cu_fit_OLS, Cu_lwr_OLS, Cu_upr_OLS, Cu_fit_WLS, Cu_lwr_WLS, Cu_upr_WLS)) %>% 
  print()
write.csv(Cu_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Cu_predict <- ACE_Cu_final + 
  geom_line(data = Cu_data_1_out, aes(y = Cu_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Cu_data_1_out, aes(y = Cu_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Cu_data_2_out, aes(y = Cu_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Cu_data_2_out, aes(y = Cu_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Cu_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Cu/Cu_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Zn -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Zn_lm <- lm(Zn_ICP ~ Zn, data = ACE_dataset)
summary(ACE_Zn_lm)
glance(ACE_Zn_lm)
model_performance(ACE_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Zn_wlm <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/(Zn_sd)^2)
summary(ACE_Zn_wlm)
glance(ACE_Zn_wlm)
model_performance(ACE_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Zn_model <- lm(Zn_ICP ~ Zn, data = ACE_dataset) # define model
ACE_Zn_wt <- 1 / lm(abs(ACE_Zn_model$residuals) ~ ACE_Zn_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Zn_wls <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zn_wls) # summary stats
glance(ACE_Zn_wls) # summary stats including AIC
model_performance(ACE_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Zn_model_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weight = 1/Zn_sd^2) # define model
ACE_Zn_wt_wt <- 1 / lm(abs(ACE_Zn_model_wt$residuals) ~ ACE_Zn_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Zn_wls_wt <- lm(Zn_ICP ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zn_wls_wt) # summary stats
glance(ACE_Zn_wls_wt) # summary stats including AIC
model_performance(ACE_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zn_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_lm_hats <- as.data.frame(hatvalues(ACE_Zn_lm))
ACE_Zn_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_lm_cooksD <- cooks.distance(ACE_Zn_lm)
ACE_Zn_lm_influential <- ACE_Zn_lm_cooksD[(ACE_Zn_lm_cooksD > (3 * mean(ACE_Zn_lm_cooksD, na.rm = TRUE)))]
ACE_Zn_lm_influential
ACE_Zn_lm_influential_names <- names(ACE_Zn_lm_influential)
ACE_Zn_lm_outliers <- ACE_dataset[ACE_Zn_lm_influential_names,] # outliers only using of index values
ACE_Zn_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Zn_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_wlm_hats <- as.data.frame(hatvalues(ACE_Zn_wlm))
ACE_Zn_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_wlm_cooksD <- cooks.distance(ACE_Zn_wlm)
ACE_Zn_wlm_influential <- ACE_Zn_wlm_cooksD[(ACE_Zn_wlm_cooksD > (3 * mean(ACE_Zn_wlm_cooksD, na.rm = TRUE)))]
ACE_Zn_wlm_influential
ACE_Zn_wlm_influential_names <- names(ACE_Zn_wlm_influential)
ACE_Zn_wlm_outliers <- ACE_dataset[ACE_Zn_wlm_influential_names,] # outliers only using of index values
ACE_Zn_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Zn_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_wls_hats <- as.data.frame(hatvalues(ACE_Zn_wls))
ACE_Zn_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_wls_cooksD <- cooks.distance(ACE_Zn_wls)
ACE_Zn_wls_influential <- ACE_Zn_wls_cooksD[(ACE_Zn_wls_cooksD > (3 * mean(ACE_Zn_wls_cooksD, na.rm = TRUE)))]
ACE_Zn_wls_influential
ACE_Zn_wls_influential_names <- names(ACE_Zn_wls_influential)
ACE_Zn_wls_outliers <- ACE_dataset[ACE_Zn_wls_influential_names,] # outliers only using of index values
ACE_Zn_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Zn_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zn_wls_wt_hats <- as.data.frame(hatvalues(ACE_Zn_wls_wt))
ACE_Zn_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zn_wls_wt_cooksD <- cooks.distance(ACE_Zn_wls_wt)
ACE_Zn_wls_wt_influential <- ACE_Zn_wls_wt_cooksD[(ACE_Zn_wls_wt_cooksD > (3 * mean(ACE_Zn_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Zn_wls_wt_influential
ACE_Zn_wls_wt_influential_names <- names(ACE_Zn_wls_wt_influential)
ACE_Zn_wls_wt_outliers <- ACE_dataset[ACE_Zn_wls_wt_influential_names,] # outliers only using of index values
ACE_Zn_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Zn_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Zn_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_OLS_summary.txt")
summary(ACE_Zn_lm)
glance(ACE_Zn_lm)
model_performance(ACE_Zn_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_lm_lev.pdf")
leveragePlots(ACE_Zn_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_OLS_wt_summary.txt")
summary(ACE_Zn_wlm)
glance(ACE_Zn_wlm)
model_performance(ACE_Zn_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wlm_lev.pdf")
leveragePlots(ACE_Zn_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_WLS_summary.txt")
summary(ACE_Zn_wls)
glance(ACE_Zn_wls)
model_performance(ACE_Zn_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wls_lev.pdf")
leveragePlots(ACE_Zn_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_WLS_wt_summary.txt")
summary(ACE_Zn_wls_wt)
glance(ACE_Zn_wls_wt)
model_performance(ACE_Zn_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zn_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zn_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Zn_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Zn_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wls_wt_lev.pdf")
leveragePlots(ACE_Zn_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/ACE_Zn_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zn_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Zn"
Zn_WLS_wt = ACE_Zn_wt
Zn_WLS_err_wt = ACE_Zn_wt_wt
theme_set(theme_classic(10))

ACE_Zn <- ggplot(ACE_dataset, aes(x = Zn, y = Zn_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Zn-Zn_sd, xmax=Zn+Zn_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Zn_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Zn_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Zn_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Zn

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zn_wlm_p <- function(ACE_Zn_wlm) {
  f <- summary(ACE_Zn_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zn_wlm_p(ACE_Zn_wlm)

ACE_Zn_wlm_eqn <- function(df){
  m <- lm(Zn_ICP ~ Zn, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zn_wlm_p(ACE_Zn_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Zn_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Zn_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Zn_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Zn_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Zn_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Zn_x.reg <- ACE_dataset$Zn
Zn_y.reg <- ACE_dataset$Zn_ICP
# Build linear regression models
Zn_model1 <- lm(Zn_y.reg ~ Zn, data = ACE_dataset) # OSL unweighted
Zn_model2 <- lm(Zn_y.reg ~ Zn, data = ACE_dataset, weights=ACE_Zn_wt) # WLS unweighted
# Model summary statistics
Zn_model1
summary(Zn_model1)
confint(Zn_model1)
Zn_model2
summary(Zn_model2)
confint(Zn_model2)
# Create predicted values and upper/lower CI to check model is working
Zn_new1.y <- data.frame(Zn_x.reg = c(1, 2, 3))
predict(Zn_model1, Zn_newdata1 = Zn_new1.y)
predict(Zn_model1, Zn_newdata1 = Zn_new1.y, Zn_interval1 = "confidence")
Zn_new2.y <- data.frame(Zn_x.reg = c(1, 2, 3))
predict(Zn_model2, Zn_newdata2 = Zn_new2.y)
predict(Zn_model2, Zn_newdata2 = Zn_new2.y, Zn_interval2 = "confidence")
# Add prediction intervals to model data frame
Zn_pred.int1 <- predict(Zn_model1, interval = "prediction")
Zn_data_1_out <- bind_cols(ACE_dataset, Zn_pred.int1) %>%
  rename(Zn_fit_OLS = fit, Zn_lwr_OLS = lwr, Zn_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Zn, Zn_sd, Zn_ICP, Zn_ICP_sd, Zn_fit_OLS, Zn_lwr_OLS, Zn_upr_OLS))
Zn_pred.int2 <- predict(Zn_model2, interval = "prediction")
Zn_data_2_out <- bind_cols(Zn_data_1_out, Zn_pred.int2) %>%
  rename(Zn_fit_WLS = fit, Zn_lwr_WLS = lwr, Zn_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Zn, Zn_sd, Zn_ICP, Zn_ICP_sd, Zn_fit_OLS, Zn_lwr_OLS, Zn_upr_OLS, Zn_fit_WLS, Zn_lwr_WLS, Zn_upr_WLS)) %>% 
  print()
write.csv(Zn_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Zn_predict <- ACE_Zn_final + 
  geom_line(data = Zn_data_1_out, aes(y = Zn_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Zn_data_1_out, aes(y = Zn_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Zn_data_2_out, aes(y = Zn_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Zn_data_2_out, aes(y = Zn_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Zn_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zn/Zn_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Rb -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Rb_lm <- lm(Rb_ICP ~ Rb, data = ACE_dataset)
summary(ACE_Rb_lm)
glance(ACE_Rb_lm)
model_performance(ACE_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Rb_wlm <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/(Rb_sd)^2)
summary(ACE_Rb_wlm)
glance(ACE_Rb_wlm)
model_performance(ACE_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Rb_model <- lm(Rb_ICP ~ Rb, data = ACE_dataset) # define model
ACE_Rb_wt <- 1 / lm(abs(ACE_Rb_model$residuals) ~ ACE_Rb_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Rb_wls <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt) #perform weighted least squares regression
# Checks
summary(ACE_Rb_wls) # summary stats
glance(ACE_Rb_wls) # summary stats including AIC
model_performance(ACE_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Rb_model_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weight = 1/Rb_sd^2) # define model
ACE_Rb_wt_wt <- 1 / lm(abs(ACE_Rb_model_wt$residuals) ~ ACE_Rb_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Rb_wls_wt <- lm(Rb_ICP ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Rb_wls_wt) # summary stats
glance(ACE_Rb_wls_wt) # summary stats including AIC
model_performance(ACE_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Rb_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_lm_hats <- as.data.frame(hatvalues(ACE_Rb_lm))
ACE_Rb_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_lm_cooksD <- cooks.distance(ACE_Rb_lm)
ACE_Rb_lm_influential <- ACE_Rb_lm_cooksD[(ACE_Rb_lm_cooksD > (3 * mean(ACE_Rb_lm_cooksD, na.rm = TRUE)))]
ACE_Rb_lm_influential
ACE_Rb_lm_influential_names <- names(ACE_Rb_lm_influential)
ACE_Rb_lm_outliers <- ACE_dataset[ACE_Rb_lm_influential_names,] # outliers only using of index values
ACE_Rb_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Rb_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_wlm_hats <- as.data.frame(hatvalues(ACE_Rb_wlm))
ACE_Rb_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_wlm_cooksD <- cooks.distance(ACE_Rb_wlm)
ACE_Rb_wlm_influential <- ACE_Rb_wlm_cooksD[(ACE_Rb_wlm_cooksD > (3 * mean(ACE_Rb_wlm_cooksD, na.rm = TRUE)))]
ACE_Rb_wlm_influential
ACE_Rb_wlm_influential_names <- names(ACE_Rb_wlm_influential)
ACE_Rb_wlm_outliers <- ACE_dataset[ACE_Rb_wlm_influential_names,] # outliers only using of index values
ACE_Rb_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Rb_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_wls_hats <- as.data.frame(hatvalues(ACE_Rb_wls))
ACE_Rb_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_wls_cooksD <- cooks.distance(ACE_Rb_wls)
ACE_Rb_wls_influential <- ACE_Rb_wls_cooksD[(ACE_Rb_wls_cooksD > (3 * mean(ACE_Rb_wls_cooksD, na.rm = TRUE)))]
ACE_Rb_wls_influential
ACE_Rb_wls_influential_names <- names(ACE_Rb_wls_influential)
ACE_Rb_wls_outliers <- ACE_dataset[ACE_Rb_wls_influential_names,] # outliers only using of index values
ACE_Rb_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Rb_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Rb_wls_wt_hats <- as.data.frame(hatvalues(ACE_Rb_wls_wt))
ACE_Rb_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Rb_wls_wt_cooksD <- cooks.distance(ACE_Rb_wls_wt)
ACE_Rb_wls_wt_influential <- ACE_Rb_wls_wt_cooksD[(ACE_Rb_wls_wt_cooksD > (3 * mean(ACE_Rb_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Rb_wls_wt_influential
ACE_Rb_wls_wt_influential_names <- names(ACE_Rb_wls_wt_influential)
ACE_Rb_wls_wt_outliers <- ACE_dataset[ACE_Rb_wls_wt_influential_names,] # outliers only using of index values
ACE_Rb_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Rb_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Rb_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_OLS_summary.txt")
summary(ACE_Rb_lm)
glance(ACE_Rb_lm)
model_performance(ACE_Rb_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_lm_lev.pdf")
leveragePlots(ACE_Rb_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_OLS_wt_summary.txt")
summary(ACE_Rb_wlm)
glance(ACE_Rb_wlm)
model_performance(ACE_Rb_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wlm_lev.pdf")
leveragePlots(ACE_Rb_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_WLS_summary.txt")
summary(ACE_Rb_wls)
glance(ACE_Rb_wls)
model_performance(ACE_Rb_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wls_lev.pdf")
leveragePlots(ACE_Rb_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_WLS_wt_summary.txt")
summary(ACE_Rb_wls_wt)
glance(ACE_Rb_wls_wt)
model_performance(ACE_Rb_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Rb_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Rb_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Rb_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Rb_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wls_wt_lev.pdf")
leveragePlots(ACE_Rb_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/ACE_Rb_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Rb_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Rb"
Rb_WLS_wt = ACE_Rb_wt
Rb_WLS_err_wt = ACE_Rb_wt_wt
theme_set(theme_classic(10))

ACE_Rb <- ggplot(ACE_dataset, aes(x = Rb, y = Rb_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Rb-Rb_sd, xmax=Rb+Rb_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Rb_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Rb_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Rb_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Rb

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Rb_wlm_p <- function(ACE_Rb_wlm) {
  f <- summary(ACE_Rb_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Rb_wlm_p(ACE_Rb_wlm)

ACE_Rb_wlm_eqn <- function(df){
  m <- lm(Rb_ICP ~ Rb, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Rb_wlm_p(ACE_Rb_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Rb_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Rb_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Rb_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Rb_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Rb_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Rb_x.reg <- ACE_dataset$Rb
Rb_y.reg <- ACE_dataset$Rb_ICP
# Build linear regression models
Rb_model1 <- lm(Rb_y.reg ~ Rb, data = ACE_dataset) # OSL unweighted
Rb_model2 <- lm(Rb_y.reg ~ Rb, data = ACE_dataset, weights=ACE_Rb_wt) # WLS unweighted
# Model summary statistics
Rb_model1
summary(Rb_model1)
confint(Rb_model1)
Rb_model2
summary(Rb_model2)
confint(Rb_model2)
# Create predicted values and upper/lower CI to check model is working
Rb_new1.y <- data.frame(Rb_x.reg = c(1, 2, 3))
predict(Rb_model1, Rb_newdata1 = Rb_new1.y)
predict(Rb_model1, Rb_newdata1 = Rb_new1.y, Rb_interval1 = "confidence")
Rb_new2.y <- data.frame(Rb_x.reg = c(1, 2, 3))
predict(Rb_model2, Rb_newdata2 = Rb_new2.y)
predict(Rb_model2, Rb_newdata2 = Rb_new2.y, Rb_interval2 = "confidence")
# Add prediction intervals to model data frame
Rb_pred.int1 <- predict(Rb_model1, interval = "prediction")
Rb_data_1_out <- bind_cols(ACE_dataset, Rb_pred.int1) %>%
  rename(Rb_fit_OLS = fit, Rb_lwr_OLS = lwr, Rb_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Rb, Rb_sd, Rb_ICP, Rb_ICP_sd, Rb_fit_OLS, Rb_lwr_OLS, Rb_upr_OLS))
Rb_pred.int2 <- predict(Rb_model2, interval = "prediction")
Rb_data_2_out <- bind_cols(Rb_data_1_out, Rb_pred.int2) %>%
  rename(Rb_fit_WLS = fit, Rb_lwr_WLS = lwr, Rb_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Rb, Rb_sd, Rb_ICP, Rb_ICP_sd, Rb_fit_OLS, Rb_lwr_OLS, Rb_upr_OLS, Rb_fit_WLS, Rb_lwr_WLS, Rb_upr_WLS)) %>% 
  print()
write.csv(Rb_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Rb_predict <- ACE_Rb_final + 
  geom_line(data = Rb_data_1_out, aes(y = Rb_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Rb_data_1_out, aes(y = Rb_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Rb_data_2_out, aes(y = Rb_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Rb_data_2_out, aes(y = Rb_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Rb_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Rb/Rb_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Sr -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Sr_lm <- lm(Sr_ICP ~ Sr, data = ACE_dataset)
summary(ACE_Sr_lm)
glance(ACE_Sr_lm)
model_performance(ACE_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Sr_wlm <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/(Sr_sd)^2)
summary(ACE_Sr_wlm)
glance(ACE_Sr_wlm)
model_performance(ACE_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Sr_model <- lm(Sr_ICP ~ Sr, data = ACE_dataset) # define model
ACE_Sr_wt <- 1 / lm(abs(ACE_Sr_model$residuals) ~ ACE_Sr_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Sr_wls <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt) #perform weighted least squares regression
# Checks
summary(ACE_Sr_wls) # summary stats
glance(ACE_Sr_wls) # summary stats including AIC
model_performance(ACE_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Sr_model_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weight = 1/Sr_sd^2) # define model
ACE_Sr_wt_wt <- 1 / lm(abs(ACE_Sr_model_wt$residuals) ~ ACE_Sr_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Sr_wls_wt <- lm(Sr_ICP ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Sr_wls_wt) # summary stats
glance(ACE_Sr_wls_wt) # summary stats including AIC
model_performance(ACE_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Sr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_lm_hats <- as.data.frame(hatvalues(ACE_Sr_lm))
ACE_Sr_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_lm_cooksD <- cooks.distance(ACE_Sr_lm)
ACE_Sr_lm_influential <- ACE_Sr_lm_cooksD[(ACE_Sr_lm_cooksD > (3 * mean(ACE_Sr_lm_cooksD, na.rm = TRUE)))]
ACE_Sr_lm_influential
ACE_Sr_lm_influential_names <- names(ACE_Sr_lm_influential)
ACE_Sr_lm_outliers <- ACE_dataset[ACE_Sr_lm_influential_names,] # outliers only using of index values
ACE_Sr_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Sr_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_wlm_hats <- as.data.frame(hatvalues(ACE_Sr_wlm))
ACE_Sr_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_wlm_cooksD <- cooks.distance(ACE_Sr_wlm)
ACE_Sr_wlm_influential <- ACE_Sr_wlm_cooksD[(ACE_Sr_wlm_cooksD > (3 * mean(ACE_Sr_wlm_cooksD, na.rm = TRUE)))]
ACE_Sr_wlm_influential
ACE_Sr_wlm_influential_names <- names(ACE_Sr_wlm_influential)
ACE_Sr_wlm_outliers <- ACE_dataset[ACE_Sr_wlm_influential_names,] # outliers only using of index values
ACE_Sr_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Sr_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_wls_hats <- as.data.frame(hatvalues(ACE_Sr_wls))
ACE_Sr_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_wls_cooksD <- cooks.distance(ACE_Sr_wls)
ACE_Sr_wls_influential <- ACE_Sr_wls_cooksD[(ACE_Sr_wls_cooksD > (3 * mean(ACE_Sr_wls_cooksD, na.rm = TRUE)))]
ACE_Sr_wls_influential
ACE_Sr_wls_influential_names <- names(ACE_Sr_wls_influential)
ACE_Sr_wls_outliers <- ACE_dataset[ACE_Sr_wls_influential_names,] # outliers only using of index values
ACE_Sr_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Sr_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Sr_wls_wt_hats <- as.data.frame(hatvalues(ACE_Sr_wls_wt))
ACE_Sr_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Sr_wls_wt_cooksD <- cooks.distance(ACE_Sr_wls_wt)
ACE_Sr_wls_wt_influential <- ACE_Sr_wls_wt_cooksD[(ACE_Sr_wls_wt_cooksD > (3 * mean(ACE_Sr_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Sr_wls_wt_influential
ACE_Sr_wls_wt_influential_names <- names(ACE_Sr_wls_wt_influential)
ACE_Sr_wls_wt_outliers <- ACE_dataset[ACE_Sr_wls_wt_influential_names,] # outliers only using of index values
ACE_Sr_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Sr_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Sr_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_OLS_summary.txt")
summary(ACE_Sr_lm)
glance(ACE_Sr_lm)
model_performance(ACE_Sr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_lm_lev.pdf")
leveragePlots(ACE_Sr_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_OLS_wt_summary.txt")
summary(ACE_Sr_wlm)
glance(ACE_Sr_wlm)
model_performance(ACE_Sr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wlm_lev.pdf")
leveragePlots(ACE_Sr_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_WLS_summary.txt")
summary(ACE_Sr_wls)
glance(ACE_Sr_wls)
model_performance(ACE_Sr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wls_lev.pdf")
leveragePlots(ACE_Sr_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_WLS_wt_summary.txt")
summary(ACE_Sr_wls_wt)
glance(ACE_Sr_wls_wt)
model_performance(ACE_Sr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Sr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Sr_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Sr_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Sr_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wls_wt_lev.pdf")
leveragePlots(ACE_Sr_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/ACE_Sr_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Sr_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Sr"
Sr_WLS_wt = ACE_Sr_wt
Sr_WLS_err_wt = ACE_Sr_wt_wt
theme_set(theme_classic(10))

ACE_Sr <- ggplot(ACE_dataset, aes(x = Sr, y = Sr_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Sr-Sr_sd, xmax=Sr+Sr_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Sr_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Sr_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Sr_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Sr

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Sr_wlm_p <- function(ACE_Sr_wlm) {
  f <- summary(ACE_Sr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Sr_wlm_p(ACE_Sr_wlm)

ACE_Sr_wlm_eqn <- function(df){
  m <- lm(Sr_ICP ~ Sr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Sr_wlm_p(ACE_Sr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Sr_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Sr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Sr_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Sr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Sr_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Sr_x.reg <- ACE_dataset$Sr
Sr_y.reg <- ACE_dataset$Sr_ICP
# Build linear regression models
Sr_model1 <- lm(Sr_y.reg ~ Sr, data = ACE_dataset) # OSL unweighted
Sr_model2 <- lm(Sr_y.reg ~ Sr, data = ACE_dataset, weights=ACE_Sr_wt) # WLS unweighted
# Model summary statistics
Sr_model1
summary(Sr_model1)
confint(Sr_model1)
Sr_model2
summary(Sr_model2)
confint(Sr_model2)
# Create predicted values and upper/lower CI to check model is working
Sr_new1.y <- data.frame(Sr_x.reg = c(1, 2, 3))
predict(Sr_model1, Sr_newdata1 = Sr_new1.y)
predict(Sr_model1, Sr_newdata1 = Sr_new1.y, Sr_interval1 = "confidence")
Sr_new2.y <- data.frame(Sr_x.reg = c(1, 2, 3))
predict(Sr_model2, Sr_newdata2 = Sr_new2.y)
predict(Sr_model2, Sr_newdata2 = Sr_new2.y, Sr_interval2 = "confidence")
# Add prediction intervals to model data frame
Sr_pred.int1 <- predict(Sr_model1, interval = "prediction")
Sr_data_1_out <- bind_cols(ACE_dataset, Sr_pred.int1) %>%
  rename(Sr_fit_OLS = fit, Sr_lwr_OLS = lwr, Sr_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Sr, Sr_sd, Sr_ICP, Sr_ICP_sd, Sr_fit_OLS, Sr_lwr_OLS, Sr_upr_OLS))
Sr_pred.int2 <- predict(Sr_model2, interval = "prediction")
Sr_data_2_out <- bind_cols(Sr_data_1_out, Sr_pred.int2) %>%
  rename(Sr_fit_WLS = fit, Sr_lwr_WLS = lwr, Sr_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Sr, Sr_sd, Sr_ICP, Sr_ICP_sd, Sr_fit_OLS, Sr_lwr_OLS, Sr_upr_OLS, Sr_fit_WLS, Sr_lwr_WLS, Sr_upr_WLS)) %>% 
  print()
write.csv(Sr_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Sr_predict <- ACE_Sr_final + 
  geom_line(data = Sr_data_1_out, aes(y = Sr_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Sr_data_1_out, aes(y = Sr_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Sr_data_2_out, aes(y = Sr_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Sr_data_2_out, aes(y = Sr_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Sr_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Sr/Sr_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_Zr -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_Zr_lm <- lm(Zr_ICP ~ Zr, data = ACE_dataset)
summary(ACE_Zr_lm)
glance(ACE_Zr_lm)
model_performance(ACE_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 2) Weighted OLS linear model & checks
ACE_Zr_wlm <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/(Zr_sd)^2)
summary(ACE_Zr_wlm)
glance(ACE_Zr_wlm)
model_performance(ACE_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 3) Unweighted Linear Regression (WLS) model
ACE_Zr_model <- lm(Zr_ICP ~ Zr, data = ACE_dataset) # define model
ACE_Zr_wt <- 1 / lm(abs(ACE_Zr_model$residuals) ~ ACE_Zr_model$fitted.values)$fitted.values^2 #define weights to use
ACE_Zr_wls <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zr_wls) # summary stats
glance(ACE_Zr_wls) # summary stats including AIC
model_performance(ACE_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_Zr_model_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weight = 1/Zr_sd^2) # define model
ACE_Zr_wt_wt <- 1 / lm(abs(ACE_Zr_model_wt$residuals) ~ ACE_Zr_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_Zr_wls_wt <- lm(Zr_ICP ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_Zr_wls_wt) # summary stats
glance(ACE_Zr_wls_wt) # summary stats including AIC
model_performance(ACE_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_Zr_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Leverage & Cooks distance

# 1) OLS (Ordinary Least Squares): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_lm_hats <- as.data.frame(hatvalues(ACE_Zr_lm))
ACE_Zr_lm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_lm_cooksD <- cooks.distance(ACE_Zr_lm)
ACE_Zr_lm_influential <- ACE_Zr_lm_cooksD[(ACE_Zr_lm_cooksD > (3 * mean(ACE_Zr_lm_cooksD, na.rm = TRUE)))]
ACE_Zr_lm_influential
ACE_Zr_lm_influential_names <- names(ACE_Zr_lm_influential)
ACE_Zr_lm_outliers <- ACE_dataset[ACE_Zr_lm_influential_names,] # outliers only using of index values
ACE_Zr_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_lm_outliers) # new dataset with outliers removed
write.csv(ACE_Zr_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_OLS_no_outliers.csv", row.names = FALSE)

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_wlm_hats <- as.data.frame(hatvalues(ACE_Zr_wlm))
ACE_Zr_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_wlm_cooksD <- cooks.distance(ACE_Zr_wlm)
ACE_Zr_wlm_influential <- ACE_Zr_wlm_cooksD[(ACE_Zr_wlm_cooksD > (3 * mean(ACE_Zr_wlm_cooksD, na.rm = TRUE)))]
ACE_Zr_wlm_influential
ACE_Zr_wlm_influential_names <- names(ACE_Zr_wlm_influential)
ACE_Zr_wlm_outliers <- ACE_dataset[ACE_Zr_wlm_influential_names,] # outliers only using of index values
ACE_Zr_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_Zr_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_OLS_wt_no_outliers.csv", row.names = FALSE)

# 3) Weighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_wls_hats <- as.data.frame(hatvalues(ACE_Zr_wls))
ACE_Zr_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_wls_cooksD <- cooks.distance(ACE_Zr_wls)
ACE_Zr_wls_influential <- ACE_Zr_wls_cooksD[(ACE_Zr_wls_cooksD > (3 * mean(ACE_Zr_wls_cooksD, na.rm = TRUE)))]
ACE_Zr_wls_influential
ACE_Zr_wls_influential_names <- names(ACE_Zr_wls_influential)
ACE_Zr_wls_outliers <- ACE_dataset[ACE_Zr_wls_influential_names,] # outliers only using of index values
ACE_Zr_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_wls_outliers) # new dataset with outliers removed
write.csv(ACE_Zr_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_WLS_no_outliers.csv", row.names = FALSE)

# 4) Error weighted eighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_Zr_wls_wt_hats <- as.data.frame(hatvalues(ACE_Zr_wls_wt))
ACE_Zr_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_Zr_wls_wt_cooksD <- cooks.distance(ACE_Zr_wls_wt)
ACE_Zr_wls_wt_influential <- ACE_Zr_wls_wt_cooksD[(ACE_Zr_wls_wt_cooksD > (3 * mean(ACE_Zr_wls_wt_cooksD, na.rm = TRUE)))]
ACE_Zr_wls_wt_influential
ACE_Zr_wls_wt_influential_names <- names(ACE_Zr_wls_wt_influential)
ACE_Zr_wls_wt_outliers <- ACE_dataset[ACE_Zr_wls_wt_influential_names,] # outliers only using of index values
ACE_Zr_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_Zr_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_Zr_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_WLS_wt_no_outliers.csv", row.names = FALSE)

# Write stats to file
# 1) OLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_OLS_summary.txt")
summary(ACE_Zr_lm)
glance(ACE_Zr_lm)
model_performance(ACE_Zr_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_lm_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_lm_lev.pdf")
leveragePlots(ACE_Zr_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_lm)
dev.off()
# 2) OLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_OLS_wt_summary.txt")
summary(ACE_Zr_wlm)
glance(ACE_Zr_wlm)
model_performance(ACE_Zr_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_lm) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_lm) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wlm_lev.pdf")
leveragePlots(ACE_Zr_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_wlm)
dev.off()
# 3) WLS - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_WLS_summary.txt")
summary(ACE_Zr_wls)
glance(ACE_Zr_wls)
model_performance(ACE_Zr_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_wls) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_wls) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wls_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wls_lev.pdf")
leveragePlots(ACE_Zr_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_wls)
dev.off()
# 4) WLS weighted - write to file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_WLS_wt_summary.txt")
summary(ACE_Zr_wls_wt)
glance(ACE_Zr_wls_wt)
model_performance(ACE_Zr_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_Zr_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_Zr_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_Zr_wls_wt) # check for random effects - returns NULL if none present
sink(file = NULL)
#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_Zr_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wls_wt_lev.pdf")
leveragePlots(ACE_Zr_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/ACE_Zr_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_Zr_wls_wt)
dev.off()

# Linear model plotting
element_title <- "Zr"
Zr_WLS_wt = ACE_Zr_wt
Zr_WLS_err_wt = ACE_Zr_wt_wt
theme_set(theme_classic(10))

ACE_Zr <- ggplot(ACE_dataset, aes(x = Zr, y = Zr_ICP)) + #ACE_dataset
  #geom_errorbar(aes(ymin=Zr_ICP-Zr_ICP_sd, ymax=Zr_ICP+Zr_ICP_sd), width=0, colour = "grey", alpha = 0.7) +
  #geom_errorbar(aes(xmin=Zr-Zr_sd, xmax=Zr+Zr_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            linetype = "dashed", aes(weight = 1/Zr_sd^2), colour = "cadetblue4") + # OLS error weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") + # OLS unweighted
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            linetype = "dashed", aes(weight = Zr_WLS_err_wt), colour="darkgrey") + # WLS weighted
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = Zr_WLS_wt), colour="black") + # WLS unweighted
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
  labs(x = paste(element_title,  itrax_dataset, " [XRF-CS]") , y = paste0(element_title, itrax_dataset, " [ICPMS]")) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
ACE_Zr

# Define p value, equation & R2 as a string to add to plots
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
# Define p value,  OLS equation & R2 as a string to add to plot
# Source: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
ACE_Zr_wlm_p <- function(ACE_Zr_wlm) {
  f <- summary(ACE_Zr_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_Zr_wlm_p(ACE_Zr_wlm)

ACE_Zr_wlm_eqn <- function(df){
  m <- lm(Zr_ICP ~ Zr, data = ACE_dataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_Zr_wlm_p(ACE_Zr_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_Zr_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "gray30", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_Zr_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 3, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_Zr_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "cadetblue4", hjust = -0.1, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_Zr_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_Zr_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Add prediction intervals
library(lmtest)
Zr_x.reg <- ACE_dataset$Zr
Zr_y.reg <- ACE_dataset$Zr_ICP
# Build linear regression models
Zr_model1 <- lm(Zr_y.reg ~ Zr, data = ACE_dataset) # OSL unweighted
Zr_model2 <- lm(Zr_y.reg ~ Zr, data = ACE_dataset, weights=ACE_Zr_wt) # WLS unweighted
# Model summary statistics
Zr_model1
summary(Zr_model1)
confint(Zr_model1)
Zr_model2
summary(Zr_model2)
confint(Zr_model2)
# Create predicted values and upper/lower CI to check model is working
Zr_new1.y <- data.frame(Zr_x.reg = c(1, 2, 3))
predict(Zr_model1, Zr_newdata1 = Zr_new1.y)
predict(Zr_model1, Zr_newdata1 = Zr_new1.y, Zr_interval1 = "confidence")
Zr_new2.y <- data.frame(Zr_x.reg = c(1, 2, 3))
predict(Zr_model2, Zr_newdata2 = Zr_new2.y)
predict(Zr_model2, Zr_newdata2 = Zr_new2.y, Zr_interval2 = "confidence")
# Add prediction intervals to model data frame
Zr_pred.int1 <- predict(Zr_model1, interval = "prediction")
Zr_data_1_out <- bind_cols(ACE_dataset, Zr_pred.int1) %>%
  rename(Zr_fit_OLS = fit, Zr_lwr_OLS = lwr, Zr_upr_OLS = upr) %>% 
  select(c(Location:midpoint, Zr, Zr_sd, Zr_ICP, Zr_ICP_sd, Zr_fit_OLS, Zr_lwr_OLS, Zr_upr_OLS))
Zr_pred.int2 <- predict(Zr_model2, interval = "prediction")
Zr_data_2_out <- bind_cols(Zr_data_1_out, Zr_pred.int2) %>%
  rename(Zr_fit_WLS = fit, Zr_lwr_WLS = lwr, Zr_upr_WLS = upr) %>% 
  select(c(Location:midpoint, Zr, Zr_sd, Zr_ICP, Zr_ICP_sd, Zr_fit_OLS, Zr_lwr_OLS, Zr_upr_OLS, Zr_fit_WLS, Zr_lwr_WLS, Zr_upr_WLS)) %>% 
  print()
write.csv(Zr_data_2_out,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_OLS_WLS_predict.csv", row.names = FALSE)

# Add OLS & WLS prediction intervals

ACE_Zr_predict <- ACE_Zr_final + 
  geom_line(data = Zr_data_1_out, aes(y = Zr_lwr_OLS), color = "blue", linetype = "dashed") +
  geom_line(data = Zr_data_1_out, aes(y = Zr_upr_OLS), color = "blue", linetype = "dashed")  +
  geom_line(data = Zr_data_2_out, aes(y = Zr_lwr_WLS), color = "black", linetype = "dashed") +
  geom_line(data = Zr_data_2_out, aes(y = Zr_upr_WLS), color = "black", linetype = "dashed")  +
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title, itrax_dataset, " [95% CI & PI]"))
ACE_Zr_predict
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/Zr/Zr_OLS_WLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")
# ACE_DM -----------------------------------------------------------------------

# 1) Unweighted OLS (Ordinary LEast Squares) - linear model & checks
ACE_DM_lm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset)
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_OLS_summary.txt")
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_DM_wlm <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/(coh_inc_sd)^2)
summary(ACE_DM_wlm)
glance(ACE_DM_wlm)
model_performance(ACE_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_OLS_wt_summary.txt")
summary(ACE_DM_wlm)
glance(ACE_DM_wlm)
model_performance(ACE_DM_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
ACE_DM_model <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset) # define model
ACE_DM_wt <- 1 / lm(abs(ACE_DM_model$residuals) ~ ACE_DM_model$fitted.values)$fitted.values^2 #define weights to use
ACE_DM_wls <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=ACE_DM_wt) #perform weighted least squares regression
# Checks
summary(ACE_DM_wls) # summary stats
glance(ACE_DM_wls) # summary stats including AIC
model_performance(ACE_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_WLS_summary.txt")
summary(ACE_DM_wls)
glance(ACE_DM_wls)
model_performance(ACE_DM_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_wls) # Performance package summary check for heteroscedasticity
icc(ACE_DM_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ITRAX error weighted - model
ACE_DM_model_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weight = 1/coh_inc_sd^2) # define model
ACE_DM_wt_wt <- 1 / lm(abs(ACE_DM_model_wt$residuals) ~ ACE_DM_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_DM_wls_wt <- lm(dry_mass_pc ~ coh_inc, data = ACE_dataset, weights=ACE_DM_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_DM_wls_wt) # summary stats
glance(ACE_DM_wls_wt) # summary stats including AIC
model_performance(ACE_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_WLS_wt_summary.txt")
summary(ACE_DM_wls_wt)
glance(ACE_DM_wls_wt)
model_performance(ACE_DM_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_DM_wls_wt) # check for random effects - returns NULL if none present
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
ACE_DM_lm_outliers <- ACE_dataset[ACE_DM_lm_influential_names,] # outliers only using of index values
ACE_DM_lm_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_lm_outliers) # new dataset with outliers removed
write.csv(ACE_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_DM_lm_lev.pdf")
leveragePlots(ACE_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_DM_wlm_hats <- as.data.frame(hatvalues(ACE_DM_wlm))
ACE_DM_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_DM_wlm_cooksD <- cooks.distance(ACE_DM_wlm)
ACE_DM_wlm_influential <- ACE_DM_wlm_cooksD[(ACE_DM_wlm_cooksD > (3 * mean(ACE_DM_wlm_cooksD, na.rm = TRUE)))]
ACE_DM_wlm_influential
ACE_DM_wlm_influential_names <- names(ACE_DM_wlm_influential)
ACE_DM_wlm_outliers <- ACE_dataset[ACE_DM_wlm_influential_names,] # outliers only using of index values
ACE_DM_wlm_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_wlm_outliers) # new dataset with outliers removed
write.csv(ACE_DM_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/OLS_ACE_DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_DM_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_DM_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_DM_wlm_lev.pdf")
leveragePlots(ACE_DM_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_DM_wls_hats <- as.data.frame(hatvalues(ACE_DM_wls))
ACE_DM_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_DM_wls_cooksD <- cooks.distance(ACE_DM_wls)
ACE_DM_wls_influential <- ACE_DM_wls_cooksD[(ACE_DM_wls_cooksD > (3 * mean(ACE_DM_wls_cooksD, na.rm = TRUE)))]
ACE_DM_wls_influential
ACE_DM_wls_influential_names <- names(ACE_DM_wls_influential)
ACE_DM_wls_outliers <- ACE_dataset[ACE_DM_wls_influential_names,] # outliers only using of index values
ACE_DM_wls_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_wls_outliers) # new dataset with outliers removed
write.csv(ACE_DM_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_wls_lev_bar.pdf")
barplot(hatvalues(ACE_DM_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_wls_lev.pdf")
leveragePlots(ACE_DM_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_DM_wls_wt_hats <- as.data.frame(hatvalues(ACE_DM_wls_wt))
ACE_DM_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_DM_wls_wt_cooksD <- cooks.distance(ACE_DM_wls_wt)
ACE_DM_wls_wt_influential <- ACE_DM_wls_wt_cooksD[(ACE_DM_wls_wt_cooksD > (3 * mean(ACE_DM_wls_wt_cooksD, na.rm = TRUE)))]
ACE_DM_wls_wt_influential
ACE_DM_wls_wt_influential_names <- names(ACE_DM_wls_wt_influential)
ACE_DM_wls_wt_outliers <- ACE_dataset[ACE_DM_wls_wt_influential_names,] # outliers only using of index values
ACE_DM_wls_wt_no_outliers <- ACE_dataset %>% anti_join(ACE_DM_wls_wt_outliers) # new dataset with outliers removed
write.csv(ACE_DM_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_DM_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_wls_wt_lev.pdf")
leveragePlots(ACE_DM_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/ACE_ACE_DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_wls_wt)
dev.off()

# Linear model summary elemental plots
element_title <- "coh_inc"
theme_set(theme_classic(10))
ACE_DM <- ggplot(ACE_dataset, aes(x = coh_inc, y = dry_mass_pc)) +
  geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/coh_inc_sd^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = ACE_DM_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = ACE_DM_wt_wt), colour="black") + # weighted WLS regression
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
  ggtitle(paste("Site: ", site_title, "; ", "Element: ", element_title))
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
  #geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_DM_wls_wt_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  #geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_DM_wlm_eqn(df)), x = -Inf, y = Inf),
  #          parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_DM_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_DM_final
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/DM/DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# Summary 4x3 matrix plots of ITRAX-ACF matched elements-------------------------------------------------------------------------

# OLS & WLS summary - unweighted stats on plot
ggarrange(ACE_K_final, ACE_Ca_final, ACE_Ti_final,
          ACE_Mn_final, ACE_Fe_final, ACE_Co_final,
          ACE_Ni_final, ACE_Cu_final, ACE_Zn_final,
          ACE_Rb_final, ACE_Sr_final, ACE_Zr_final,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/ACE_OLS_WLS_Summary.pdf",
       height = c(50), width = c(50), dpi = 600, units = "cm")

# OLS & WLS summary - unweighted stats on plot
ggarrange(ACE_K_predict, ACE_Ca_predict, ACE_Ti_predict,
          ACE_Mn_predict, ACE_Fe_predict, ACE_Co_predict,
          ACE_Ni_predict, ACE_Cu_predict, ACE_Zn_predict,
          ACE_Rb_predict, ACE_Sr_predict, ACE_Zr_predict,
          ncol = 3, nrow = 4, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/OLS_WLS_clr/ACE_OLS_WLS_Summary_predict.pdf",
       height = c(50), width = c(50), dpi = 600, units = "cm")#
# -------------------------------------------------------------------------






END

# -------------------------------------------------------------------------
# OLS Linear models - forced through origin -------------------------------

# K unweighted ---------------------------
ACE_K_origin <- ggplot(ACE_all, aes(x = K, y = K_ICP)) +
  #geom_errorbar(aes(ymin=K_ICP-K_ICP_sd, ymax=K_ICP+K_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/K_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")
# Ca unweighted ---------------------------
ACE_Ca_origin <- ggplot(ACE_all, aes(x = Ca, y = Ca_ICP)) +
  #geom_errorbar(aes(ymin=Ca_ICP-Ca_ICP_sd, ymax=Ca_ICP+Ca_ICP_sd), width=.1, colour = "grey") +
  geom_errorbar(aes(xmin=Ca-Ca_sd, xmax=Ca+Ca_sd), width=.1, colour = "grey") +
  geom_smooth(method=lm, formula=y~x-1, se=TRUE, fullrange=FALSE,  colour = "blue", aes(weight = 1/Ca_ICP_sd)) +
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
  labs(x = "Ca (XRF-CS) / cps", y = "Ca (ICPMS) / ppm" ) +
  ggtitle("ACE: Ca (ICP SD weighted = blue; unweighted = black)")
ACE_Ca_origin
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Ca_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Ti 

# Ti unweighted ---------------------------
ACE_Ti_origin <- ggplot(ACE_all, aes(x = Ti, y = Ti_ICP)) +
  #geom_errorbar(aes(ymin=Ti_ICP-Ti_ICP_sd, ymax=Ti_ICP+Ti_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Ti_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Mn

# Mn unweighted ---------------------------
ACE_Mn_origin <- ggplot(ACE_all, aes(x = Mn, y = Mn_ICP)) +
  #geom_errorbar(aes(ymin=Mn_ICP-Mn_ICP_sd, ymax=Mn_ICP+Mn_ICP_sd), width=.1, colour = "grey") +
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
  #geom_errorbar(aes(ymin=Fe_ICP-Fe_ICP_sd, ymax=Fe_ICP+Fe_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Fe_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Co

# Co unweighted ---------------------------
ACE_Co_origin <- ggplot(ACE_all, aes(x = Co, y = Co_ICP)) +
  #geom_errorbar(aes(ymin=Co_ICP-Co_ICP_sd, ymax=Co_ICP+Co_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Co_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Ni

# Ni unweighted ---------------------------
ACE_Ni_origin <- ggplot(ACE_all, aes(x = Ni, y = Ni_ICP)) +
  #geom_errorbar(aes(ymin=Ni_ICP-Ni_ICP_sd, ymax=Ni_ICP+Ni_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Ni_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Cu

# Cu unweighted ---------------------------
ACE_Cu_origin <- ggplot(ACE_all, aes(x = Cu, y = Cu_ICP)) +
  #geom_errorbar(aes(ymin=Cu_ICP-Cu_ICP_sd, ymax=Cu_ICP+Cu_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Cu_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Zn

# Zn unweighted ---------------------------
ACE_Zn_origin <- ggplot(ACE_all, aes(x = Zn, y = Zn_ICP)) +
  #geom_errorbar(aes(ymin=Zn_ICP-Zn_ICP_sd, ymax=Zn_ICP+Zn_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Zn_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Rb

# Rb unweighted ---------------------------
ACE_Rb_origin <- ggplot(ACE_all, aes(x = Rb, y = Rb_ICP)) +
  #geom_errorbar(aes(ymin=Rb_ICP-Rb_ICP_sd, ymax=Rb_ICP+Rb_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Rb_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Sr

# Sr unweighted ---------------------------
ACE_Sr_origin <- ggplot(ACE_all, aes(x = Sr, y = Sr_ICP)) +
  #geom_errorbar(aes(ymin=Sr_ICP-Sr_ICP_sd, ymax=Sr_ICP+Sr_ICP_sd), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Sr_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# DM

# DM unweighted ---------------------------
ACE_DM_origin <- ggplot(ACE_all, aes(x = coh_inc, y = dry_mass_pc)) +
  #geom_errorbar(aes(ymin=dry_mass_pc-dry_mass_err, ymax=dry_mass_pc+dry_mass_err), width=.1, colour = "grey") +
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/DM_origin_lm_wlm.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")





# -------------------------------------------------------------------------
# Performance tests & stats - forced through origin - lm and wlm ------------------------

# K

#Performance: K unweighted lm stats forced through origin
ACE_K_wlm_origin <- lm(K_ICP ~ K-1, data = ACE_all)
summary(ACE_K_wlm_origin) # summary stats
check_model(ACE_K_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/K_origin_performance_lm.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")
ACE_K_wlm_origin_mp <- model_performance(ACE_K_wlm_origin) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
display(ACE_K_wlm_origin_mp)
bptest(ACE_K_wlm_origin) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_K_wlm_origin) # Performance package summary check for heteroscedasticity
icc(ACE_K_wlm)

#Performance: K weighted lm forced through origin
ACE_K_wlm_origin <- lm(K_ICP ~ K-1, weight = 1/K_ICP_sd, data = ACE_all)
summary(ACE_K_wlm_origin) # summary stats

check_model(ACE_K_wlm_origin) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/K_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Mn_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Mn_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Ti_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Ti_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Mn_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Mn_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Fe_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Fe_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Co_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Co_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Ni_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Ni_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Cu_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Cu_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Zn_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Zn_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Rb_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Rb_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Sr_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/Sr_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/DM_origin_performance_lm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/origin/DM_origin_performance_wlm.pdf", 
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/clr/ACE_OLS_Summary_origin_lm_wlm.pdf", 
       height = c(50), width = c(50), dpi = 600, units = "cm")







# END