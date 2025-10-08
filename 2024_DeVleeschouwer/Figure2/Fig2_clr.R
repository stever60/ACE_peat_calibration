# Figure 2 and Figure S3 - Correlation - clr

# 8/9/25: depth matching now defined by Subsample-ICPMS matched composite dataset
# n = 268 - increased from n = 265 after reduction from n = 272 - corrections to composite depths in BI10 10/6/25 and then KER3 8/9/25
# Code adapted from "Papers_R/2024_DeVleeschouwer/Output/itrax_Composite/Matching_mean/ACE/ACE_matching_clr.R" 

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
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis', 'psych' )
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

# Define elements to use -------------------------------------------------------

# ICP elements of interest defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "Dry_mass")

# XRF-CS acf elements matched to Francois ICPMS elements
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")
acf_icp_Elements_only_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_inc_sd", "Mo_coh_sd")
acf_icp_Elements_only_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd")

# Scatter parameters
scatter_param <- c("Mo_inc", "Mo_coh", "inc_coh", "Ln_inc_coh", "coh_inc", "Ln_coh_inc", "Total_scatter",
                   "Mo_inc_sd", "Mo_coh_sd", "inc_coh_sd", "Ln_inc_coh_sd", "coh_inc_sd", "Ln_coh_inc_sd", "Total_scatter_sd")

# ICPMS elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# key elements to simplify plotting 
xrf_icp_Elements_key <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP",
                          "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "Dry_mass")

# key elements_reduced for more simplified plots 
xrf_icp_Elements_key_reduced <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")

# key elements for Figure 2a
xrf_icp_Elements_key_Fig2 <- c("Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP", "Fe", "Fe_ICP", 
                               "Sr", "Sr_ICP", "Zr", "Zr_ICP", "coh_inc", "Dry_mass")

# Subsample parameters
subsample_param <- c("Water_Content", "Dry_mass", "Dry_mass_err", 
                     "LOI550", "LOI550_err", "C_org_pc", "C_org_pc_err", 
                     "Wet_density","Dry_density",	"Dry_density_err", 
                     "DMAR", "DMAR_err", "DMAR_err_pc")

subsample_param1 <- c("Dry_mass","LOI550", "C_org_pc", "Dry_density")


# # RUN ONCE - Import existing ACE cps matched file  ----------------------------------------
# ACE_subsample_xrf_icp_matched_cps_fr_err <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_Subsample_ICPMS_XRF_matched.csv") %>% 
#   select(c(Location:DMAR_err_pc, all_of(icp_Elements_min), all_of(icp_Elements_min_sd))) %>% 
#   mutate(K_sd_fr = (K_sd/K)) %>% 
#   mutate(Ca_sd_fr = (Ca_sd/Ca)) %>%
#   mutate(Ti_sd_fr = (Ti_sd/Ti)) %>%
#   mutate(Mn_sd_fr = (Mn_sd/Mn)) %>%
#   mutate(Fe_sd_fr = (Fe_sd/Fe)) %>%
#   mutate(Co_sd_fr = (Co_sd/Co)) %>%
#   mutate(Ni_sd_fr = (Ni_sd/Ni)) %>%
#   mutate(Cu_sd_fr = (Cu_sd/Cu)) %>%
#   mutate(Zn_sd_fr = (Zn_sd/Zn)) %>%
#   mutate(Rb_sd_fr = (Rb_sd/Rb)) %>%
#   mutate(Sr_sd_fr = (Sr_sd/Sr)) %>%
#   mutate(Zr_sd_fr = (Zr_sd/Zr))
# is.na(ACE_subsample_xrf_icp_matched_cps_fr_err)<-sapply(ACE_subsample_xrf_icp_matched_cps_fr_err, is.infinite) # replace any infinite values with NA
# ACE_subsample_xrf_icp_matched_cps_fr_err
# write.csv(ACE_subsample_xrf_icp_matched_cps_fr_err,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_xrf_icp_matched_cps_fr_err.csv", row.names = FALSE)
# 
# # Replace zeros with half minimum value for each element 
# ACE_all0 <- ACE_subsample_xrf_icp_matched_cps_fr_err %>% 
#   mutate(across(all_of(c(icp_Elements_min, icp_Elements_min_sd)), ~ ifelse(.x < 0, 0, .x))) %>%  # replace any ICPMS values <0 with zero
#   #zeroreplace() %>% # replace zeros as outlined in the compositions package
#   mutate_at(vars(all_of(c(icp_Elements_min, acf_icp_Elements_min))), ## replace zeros with half minimum value to allow logs
#             ~ (. == 0) * min(.[. != 0])/2 + .) %>% # Recommended procedure from Bertrand et al. (submitted) - retains dataframe structure
#   mutate_at(vars(all_of(c(icp_Elements_min_sd, acf_icp_Elements_min_sd))), 
#             ~ (. == 0) * min(.[. != 0])/2 + .) %>% 
#   na.omit() %>% #remove rows with NAs
#   print()
# 
# # RUN ONCE - Create clr dataset -----------------------------------------------------------
# ACE_clr_itrax <- ACE_all0 %>% 
#   select(all_of(acf_icp_Elements_min)) %>%
#   select(-c(Mo_inc, Mo_coh)) %>% 
#     clr() %>% # make centred log ratio dataframe
#   as_tibble()
# 
# ACE_cps_fract_err <- ACE_all0 %>% 
#   select(K_sd_fr:Zr_sd_fr)
# ACE_clr_itrax_sd  <- bind_cols(ACE_clr_itrax, ACE_cps_fract_err) %>%
#   #mutate(across(where(is.numeric), ~ .x * 0.1)) %>% # assumes measurement errors of 0.1%.- not correct - each row has Stdev errors which si much greater
#   #rename_with(.fn = function(.x){paste0(.x,"_sd")})
# # calculate percentage errors for XRF from Stdev fractional error data above 
#   mutate(K_sd = K*K_sd_fr) %>%
#   relocate(K_sd, .after = K) %>%
#   mutate(Ca_sd = Ca*Ca_sd_fr) %>% 
#   relocate(Ca_sd, .after = Ca) %>%
#   mutate(Ti_sd = Ti*Ti_sd_fr) %>%
#   relocate(Ti_sd, .after = Ti) %>%
#   mutate(Mn_sd = Mn*Mn_sd_fr) %>% 
#   relocate(Mn_sd, .after = Mn) %>%
#   mutate(Fe_sd = Fe*Fe_sd_fr) %>%
#   relocate(Fe_sd, .after = Fe) %>%
#   mutate(Co_sd = Co*Co_sd_fr) %>%
#   relocate(Co_sd, .after = Co) %>%
#   mutate(Ni_sd = Ni*Ni_sd_fr) %>%
#   relocate(Ni_sd, .after = Ni) %>%
#   mutate(Cu_sd = Cu*Cu_sd_fr) %>%
#   relocate(Cu_sd, .after = Cu) %>%
#   mutate(Zn_sd = Zn*Zn_sd_fr) %>% 
#   relocate(Zn_sd, .after = Zn) %>%
#   mutate(Rb_sd = Rb*Rb_sd_fr) %>%
#   relocate(Rb_sd, .after = Rb) %>%
#   mutate(Sr_sd = Sr*Sr_sd_fr) %>% 
#   relocate(Sr_sd, .after = Sr) %>% 
#   mutate(Zr_sd = Zr*Zr_sd_fr) %>% 
#   relocate(Zr_sd, .after = Zr) %>%
#   mutate(across(all_of(acf_icp_Elements_only_min_sd), abs)) %>% # change negative SD error values to absolute (positive) values
#   select(!(all_of(acf_icp_Elements_only_min)))
#   
# ACE_clr_icp <- ACE_all0 %>% 
#   select(all_of(icp_Elements_min)) %>%
#   clr() %>% # make centred log ratio dataframe
#   as_tibble()
# 
# ACE_clr_icp_sd <- ACE_clr_icp %>% # % errors from Table 2 in paper
#   mutate(K_ICP_sd = K_ICP*0.18) %>%
#   relocate(K_ICP_sd, .after = K_ICP) %>%
#   mutate(Ca_ICP_sd = Ca_ICP*0.13) %>% 
#   relocate(Ca_ICP_sd, .after = Ca_ICP) %>%
#   mutate(Ti_ICP_sd = Ti_ICP*0.14) %>%
#   relocate(Ti_ICP_sd, .after = Ti_ICP) %>%
#   mutate(Mn_ICP_sd = Mn_ICP*0.12) %>% 
#   relocate(Mn_ICP_sd, .after = Mn_ICP) %>%
#   mutate(Fe_ICP_sd = Fe_ICP*0.18) %>%
#   relocate(Fe_ICP_sd, .after = Fe_ICP) %>%
#   mutate(Co_ICP_sd = Co_ICP*0.07) %>%
#   relocate(Co_ICP_sd, .after = Co_ICP) %>%
#   mutate(Ni_ICP_sd = Ni_ICP*0.54) %>%
#   relocate(Ni_ICP_sd, .after = Ni_ICP) %>%
#   mutate(Cu_ICP_sd = Cu_ICP*0.14) %>%
#   relocate(Cu_ICP_sd, .after = Cu_ICP) %>%
#   mutate(Zn_ICP_sd = Zn_ICP*0.09) %>% 
#   relocate(Zn_ICP_sd, .after = Zn_ICP) %>%
#   mutate(Rb_ICP_sd = Rb_ICP*0.16) %>%
#   relocate(Rb_ICP_sd, .after = Rb_ICP) %>%
#   mutate(Sr_ICP_sd = Sr_ICP*0.02) %>% 
#   relocate(Sr_ICP_sd, .after = Sr_ICP) %>% 
#   mutate(Zr_ICP_sd = Zr_ICP*0.05) %>% 
#   relocate(Zr_ICP_sd, .after = Zr_ICP) %>% 
#   select(K_ICP_sd, Ca_ICP_sd, Ti_ICP_sd, Mn_ICP_sd, Fe_ICP_sd, Co_ICP_sd,
#          Ni_ICP_sd, Cu_ICP_sd, Zn_ICP_sd, Rb_ICP_sd, Sr_ICP_sd, Zr_ICP_sd)
# 
# # Make other parameters / variables 
# ACE_text <- ACE_all0 %>% 
#   select(Location:MSE)
# ACE_rest <- ACE_all0 %>% 
#   select(all_of(scatter_param), all_of(subsample_param))
# 
# # Bind into single file & define dataset to use for linear modelling
# ACE_subsample_icp_xrf_matched_clr <- bind_cols(ACE_text, ACE_clr_itrax, ACE_clr_itrax_sd, 
#                      ACE_clr_icp, ACE_clr_icp_sd, ACE_rest) %>%
#   relocate(scatter_param, .after = Zr_sd) %>% 
#   print()
# write.csv(ACE_subsample_icp_xrf_matched_clr,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_clr.csv", row.names = FALSE)
# write.csv(ACE_subsample_icp_xrf_matched_clr,"Papers_R/2024_DeVleeschouwer/FigureS6/Data/Input/ACE_subsample_icp_xrf_matched_clr.csv", row.names = FALSE)
# 
# 
# 
# 
# 






# Import existing ACE log_inc matched file -----------------------------------------
ACE_subsample_icp_xrf_matched_clr <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_clr.csv")

# Define dataset to use for Linear modelling  -----------------------------
ACE_LM1 <- ACE_subsample_icp_xrf_matched_clr

# Produce stats for info -----------------------------------------------
ACE_LM1_stats <- ACE_LM1 %>%
  select(accrate: Zr_ICP_sd) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_LM1_stats,"Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_xrf_icp_matched_stats_clr.csv", row.names = FALSE)

# Convert Site to use as a grouping variable
ACE_LM1$Site <- as.factor(ACE_LM1$Site)

# Fig 2 & S3 - Correlation matrices & density plots-----------------------------

# Fig 2a - ITRAX & ICP correlation matrix - key elements reduced
theme_set(theme_bw(base_size=2))
ggcorr(ACE_LM1[,xrf_icp_Elements_key_Fig2], method = c("everything", "pearson"),
       size = 7, label = TRUE, label_alpha = FALSE, label_round=2, label_size= 7)
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2a_Corr_matrix_key_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Fig S3 - ITRAX & ICP - key elements correlation density plots
theme_set(theme_bw(base_size=8))
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key_Fig2, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="ACE XRF-CS & ICPMS: Correlation-density plot - cps")
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/FigS3.1_itrax_ICP_corr-den_matrix_clr.pdf", 
       height = c(15), width = c(15), dpi = 600, units = "cm")

# Set up function to plot each site as different colour using jco scheme
theme_set(theme_bw(base_size=8))
palette_set <- "jco" # or "jco", or "npg", "uchicago" - set up colour scheme for site plots
cor_plot <- function(data, mapping, ...) {
  
  ggplot(data, mapping) + 
    geom_point(size = 0.7) +
    geom_smooth(formula = y~x, method = lm, color = "black") + 
    scale_color_jco()
}

# Run ggpairs - black regression line all - correlation all sites
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key_Fig2, 
        title="ACE XRF-CS & ICPMS: Correlation-density plot - cps",
        upper = list(
          mapping = aes(color=Site, palette = palette_set),
          continuous = wrap(ggally_cor, size = 2)
        ),
        lower=list(
          mapping = aes(color=Site, palette = palette_set, alpha = 0.5),
          continuous = cor_plot
        )
)
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/FigS3.2_itrax_ICP_corr-den_sites_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Run ggpairs - black regression line as correlation overall - key
ggpairs(ACE_LM1, columns = xrf_icp_Elements_key_Fig2,
        title="ACE XRF-CS & ICPMS: Correlation-density plot - cps",
        upper = list(
          continuous = wrap(ggally_cor, size = 5)
        ),
        lower=list(
          mapping = aes(color=Site, palette = palette_set, alpha = 0.5),
          continuous = cor_plot
        )
)
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/FigS3.3_itrax_ICP_corr-den_sites_all_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Correlation & data plots -----------------------------------------------------------

# Define text for titles & labels for plotting
XRF_title <- " clr [XRF-CS]"
ICP_title <- " clr [ICPMS]"
correlation_title <- "clr Correlation"
method_title <- "clr"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Individual element plots -----------------------------------------------------------

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
ggsave("Papers_R/2024_DeVleeschouwer/Figure2/Plots/Fig2c_Linear_reg_per_site_clr.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

