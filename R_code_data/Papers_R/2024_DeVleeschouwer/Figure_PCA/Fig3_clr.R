# Figure 3 - clr PCA plots

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

# Figure 3 - Multivariate analysis -------------------------------------------------------------------------

# Define as adjacent matched elements for PCA
xrf_icp_Elements_min_PCA <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                              "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Co", "Co_ICP",
                              "Ni", "Ni_ICP", "Cu", "Cu_ICP", "Zn", "Zn_ICP",
                              "Rb", "Rb_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP", 
                              "Mo_inc", "Mo_coh", "coh_inc", "dry_mass_pc")

xrf_icp_Elements_min_PCA_edited <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", 
                                     "Mn", "Mn_ICP", "Fe", "Fe_ICP", "Sr", "Sr_ICP", "Zr", "Zr_ICP",
                                     "coh_inc", "dry_mass_pc")

# Import ACE matched ICPMS and ITRAX cps dataset & remove site POB4 data - clr is performed in PCA
ACE_mv0 <- read_csv("Papers_R/2024_DeVleeschouwer/Figure3/Data/clr/ACE_xrf_icp_matched_cps.csv") %>%  
  select(Location:MSE, all_of(xrf_icp_Elements_min_PCA)) %>% 
  filter(!Site =="POB4") #remove POB4 data
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
write.csv(ACE_mv,"Papers_R/2024_DeVleeschouwer/Figure3/Data/clr/ACE_multivariate_cps.csv", row.names = FALSE)

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

# PCA - summary plots - model = ~ princomp(clr(.))

# Fig 3a
MyElements <- xrf_icp_Elements_min_PCA
#MyElements <- xrf_icp_Elements_min_PCA_edited

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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_PCA_key_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Fig 3b
#MyElements <- xrf_icp_Elements_min_PCA
MyElements <- xrf_icp_Elements_min_PCA_edited

# Define titles & labels for plotting
XRF_title <- "/ inc [Ln XRF-CS]"
ICP_title <- " [Ln ICPMS]"
correlation_title <- "Log Correlation"
palette_set <- "jco" # or "npg", "uchicago"

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
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  #scale_color_continuous(type = "viridis", trans = "reverse")
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_PCA_key_edited_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")



