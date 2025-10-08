# Figure S11 - PCA clr 

# Set up ------------------------------------------------------------------

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
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis',
              'ggrepel', 'itraxR', 'PeriodicTable', 'errors', 'forecast', 'broom',
              'directlabels', 'performance', 'lmtest', 'ggpmisc', 'cowplot', 'Hmisc')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)



# Figure S11 - Multivariate PCA analysis for PLS regression --------------------
# Define elements for PCA -----------------------------------------------
xrf_icp_Elements_min_PCA <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                              "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                              "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                              "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                              "Mo_inc", "Mo_coh", "inc_coh", "coh_inc", "C_org_pc", "Dry_mass")

xrf_icp_Elements_min_PCA_edited <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                                     "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Sr", "Sr_ICP", 
                                     "Zr", "Zr_ICP", "inc_coh", "coh_inc", "C_org_pc", "Dry_mass")

xrf_icp_Elements_min_PCA_edited_no_scatter <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                                     "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP")

# Import ACE Subsample-ICPMS_XRF-cps matched file from Figure 2 -----------------------------------------------------------------------
ACE_subsample_icp_xrf_matched_cps <-read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv")
is.na(ACE_subsample_icp_xrf_matched_cps)<-sapply(ACE_subsample_icp_xrf_matched_cps, is.infinite) # replace any infinite values with NA
ACE_subsample_icp_xrf_matched_cps

# clr is created in PCA analysis from cps data

# Select elements for PCA - all or key ones from correlation/regression analysis
ACE_mv0 <- ACE_subsample_icp_xrf_matched_cps %>%  
  select(Location:MSE, all_of(xrf_icp_Elements_min_PCA))
  #select(Location:MSE, all_of(xrf_icp_Elements_min_PCA)) # or 6 major elements (5 with good correlation) & Zr for dust fluxes - this step is set in PCA
ACE_mv0

# Create site-specific unique depths for plotting sequentially ---------------
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
write.csv(ACE_mv,"Papers_R/2024_DeVleeschouwer/FigureS11/Data/ACE_PCA_cps.csv", row.names = FALSE)

TRUE %in% ACE_mv$uid %>% duplicated() # check nothing duplicated


# Transform data ----------------------------------------------------------

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
# Different to hlaf min value replacement used in clr dataset for Figure 3 - but zeroreplace() is still a very small number so insignificant  
ACE_mv_acomp %>%
  zeroreplace() %>%
  princomp() %>%
  biplot(xlabs = rep("o",times = nrow(ACE_mv_acomp)))


# PCA - summary plots - model = ~ princomp(clr(.)) ----------------------


# Figure S11a - all matched elements & scatter ---------------------------------

MyElements <- xrf_icp_Elements_min_PCA
# MyElements <- xrf_icp_Elements_min_PCA_edited
# MyElements <- xrf_icp_Elements_min_PCA_edited_no_scatter

theme_set(theme_bw(base_size = 16))

# custom colours
cluster_col5 <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF")

ACE_mv_acomp_meta %>%
  ordr::ordinate(
    cols = any_of(MyElements),
    model = ~ princomp(clr(.)),
    augment = any_of(c("site_depth", "Site"))
  ) %>%
  ordr::ggbiplot(sec.axes = "cols", scale.factor = 8) +
  # filled circles, outline same as fill
  ordr::geom_rows_point(
    aes(fill = Site, colour = Site),
    shape = 21,
    size = 1,
    stroke = 1,
    alpha = 1
  ) +
  ordr::geom_cols_vector() +
  ordr::geom_cols_text_radiate(aes(label = name)) +
  scale_fill_manual(values = cluster_col5, name = "Site") +
  scale_colour_manual(values = cluster_col5, guide = "none") +  # same colors for outline
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(
    legend.text  = element_text(size = 16),
    text         = element_text(size = 16, face = "plain"),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16, face = "plain"),
    plot.margin  = unit(c(1, 1, 1, 1), "cm")
  )
ggsave("Papers_R/2024_DeVleeschouwer/FigureS11/Plots/Fig11a_PCA_key_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")



# Figure S11b - all major elements & Zr - use to inform PLS ------------------

#MyElements <- xrf_icp_Elements_min_PCA
MyElements <- xrf_icp_Elements_min_PCA_edited
#MyElements <- xrf_icp_Elements_min_PCA_edited_no_scatter

# Define titles & labels for plotting
XRF_title <- "/ inc [Ln XRF-CS]"
ICP_title <- " [Ln ICPMS]"
correlation_title <- "Log Correlation"
palette_set <- "jco" # or "npg", "uchicago"

theme_set(theme_bw(base_size = 16))

# custom colours
cluster_col5 <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF")

ACE_mv_acomp_meta %>%
  ordr::ordinate(
    cols = any_of(MyElements),
    model = ~ princomp(clr(.)),
    augment = any_of(c("site_depth", "Site"))
  ) %>%
  ordr::ggbiplot(sec.axes = "cols", scale.factor = 8) +
  # filled circles, outline same as fill
  ordr::geom_rows_point(
    aes(fill = Site, colour = Site),
    shape = 21,
    size = 1,
    stroke = 1,
    alpha = 1
  ) +
  ordr::geom_cols_vector() +
  ordr::geom_cols_text_radiate(aes(label = name)) +
  scale_fill_manual(values = cluster_col5, name = "Site") +
  scale_colour_manual(values = cluster_col5, guide = "none") +  # same colors for outline
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(
    legend.text  = element_text(size = 16),
    text         = element_text(size = 16, face = "plain"),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16, face = "plain"),
    plot.margin  = unit(c(1, 1, 1, 1), "cm")
  )
ggsave("Papers_R/2024_DeVleeschouwer/FigureS11/Plots/Fig11b_PCA_key_edited_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure S11c - all major elements & Zr without scatter ratios, dry mass and Corg - use to inform PLS ------------------

#MyElements <- xrf_icp_Elements_min_PCA
#MyElements <- xrf_icp_Elements_min_PCA_edited
MyElements <- xrf_icp_Elements_min_PCA_edited_no_scatter

# Define titles & labels for plotting
XRF_title <- "/ inc [Ln XRF-CS]"
ICP_title <- " [Ln ICPMS]"
correlation_title <- "Log Correlation"
palette_set <- "jco" # or "npg", "uchicago"

theme_set(theme_bw(base_size = 16))

# custom colours
cluster_col5 <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF")

ACE_mv_acomp_meta %>%
  ordr::ordinate(
    cols = any_of(MyElements),
    model = ~ princomp(clr(.)),
    augment = any_of(c("site_depth", "Site"))
  ) %>%
  ordr::ggbiplot(sec.axes = "cols", scale.factor = 8) +
  # filled circles, outline same as fill
  ordr::geom_rows_point(
    aes(fill = Site, colour = Site),
    shape = 21,
    size = 1,
    stroke = 1,
    alpha = 1
  ) +
  ordr::geom_cols_vector() +
  ordr::geom_cols_text_radiate(aes(label = name)) +
  scale_fill_manual(values = cluster_col5, name = "Site") +
  scale_colour_manual(values = cluster_col5, guide = "none") +  # same colors for outline
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(
    legend.text  = element_text(size = 16),
    text         = element_text(size = 16, face = "plain"),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16, face = "plain"),
    plot.margin  = unit(c(1, 1, 1, 1), "cm")
  )
ggsave("Papers_R/2024_DeVleeschouwer/FigureS11/Plots/Fig11c_PCA_key_edited_no_scatter_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")




# # END ------------------------------------------------------------------------


