# ITRAX-ICPMS Calibration

# Step 1: ITRAX Quality Control: Spectra & Element Filtering
# Step 2: ITRAX COMPOSITE data conversion in R (validity filtered cps, %cps sum, inc normalised & as Z-scores, cps filtering for first pass clr)
  
# SECTION 1: Set up ------------------------------------------------------------

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
              'ggrepel', 'itraxR','PeriodicTable','errors','patchwork',
              'forecast','directlabels','broom','performance','lmtest','ggpmisc',
              'cowplot','Hmisc')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

# Set up ------------------------------------------------------------------

# Set working directory - Macbook Pro 2013
#setwd("/Users/Steve/Dropbox/BAS/Data/R/")
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()


# SECTION 2A: IMPORT ALL ACE_SHW DATA -----------------------------------------------------

# Import existing composite datasets where validity = 1
# itrax.R qc is carried out on all data except POB4 (measured under different conditions)

# Assign elements to/from PeriodicTable package to 'elementsList' (then unload package)
data(periodicTable)
elementsList <- periodicTable$symb

# list of all possible elements only
data(periodicTable)
allelements <- c(symb(1:117))
rm(periodicTable)

# # Import an existing composite Site and change column names to match itrax.R names/formatting
ACE_COMP <- read_csv("Papers_R/2024_DeVleeschouwer/Data/ACE_SHW_ITRAX_Composite_raw_cps.csv") 
ACE_COMP 

ACE_xrf <- ACE_COMP %>% 
  filter(validity=='1') %>% # filter to remove validity = FALSE rows
  relocate(Mo_inc, .after = U) %>%
  relocate(Mo_coh, .after = Mo_inc) %>% 
  mutate(cps = rowSums(across(Mg:Mo_coh))) %>% #mutate(cps = rowSums(.[df1_rowsums])) %>%
  mutate(Total_scatter = rowSums(across(Mo_inc:Mo_coh))) %>%
  mutate(inc_coh = Mo_inc/Mo_coh) %>%
  mutate(coh_inc = Mo_coh/Mo_inc) %>%
  select(-c(filename, position_mm, `accrate_cmyr-1`,
            strat_depth_top, SH20_min_age_95CI: SH20_median_age, 
            E_gain:F_offset, S1:S3)) %>% 
  rename(`Fe a*2` = D1, SH20_age = SH20_mean_age, label = Section, surface = `sample_surface`) %>% 
  relocate(cps, .before = MSE) %>% 
  relocate(`Fe a*2` , .after = coh_inc)
  
ACE_xrf
write.csv(ACE_xrf,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section2/ACE_xrf_cps.csv", row.names = FALSE)

# SECTION 2B: Create element lists -----------------------------------------------------------

# Remove detector elements and known matrix effect elements from dataframe 
machine_elements <- select(ACE_xrf, c(Ar, Ta, W, Ir, Pt, Au, Hg)) %>% 
  names()
machine_elements

# Create list of all elements (removing machine elements and Fe a*2) 
ACE_elements <- select(ACE_xrf, c(Mg:Mo_coh), -c(all_of(machine_elements), `Fe a*2`)) %>% # , Zr
  names()
ACE_elements

# All elements and scatter parameters
ACE_elements1 <- select(ACE_xrf, c(Mg:Mo_coh, Total_scatter, inc_coh, coh_inc), -c(all_of(machine_elements), `Fe a*2`)) %>% # , Zr
  names()
ACE_elements1

# Summary stats for cps - all elements and scatter parameters
ACE_xrf_stats <- ACE_xrf %>%
  select(any_of(ACE_elements1)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(ACE_xrf_stats)
write.csv(ACE_xrf_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section2/ACE_xrf_cps_stats.csv", row.names = FALSE)

# SECTION 3: Quality control & Element Filtering -------------------------------------------------------------------------

# Import ACE cps dataset and remove POB4 and BI5 data
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
ACE_xrf <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section2/ACE_xrf_cps.csv") %>% 
  filter(!(Site =='POB4')) %>% #remove POB4 data which was measured under different conditions
  filter(!(Site =='BI5')) %>% #site not in the matched dataset
  filter(!(Site =='DRPB')) %>% #site not in the matched dataset
  select(Location:MSE, any_of(ACE_elements), `Fe a*2`, Total_scatter:coh_inc)
ACE_xrf

# SECTION 3A: Quality control

# Apply data filtering criteria and quality control measures - from itrax.R

# MSE and cps filtering --------------------------------------------------

# cps filtering - Fe a*2
Fig3.1A <- 
  ggplot(ACE_xrf, mapping = aes(x = cps, y = `Fe a*2`)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()
print(Fig3.1A)

# cps - 2 std dev is too strict for HER42PB and other peat cores with a combination of low and high count matrices
cps.mean <- mean(ACE_xrf$cps)
cps.sd <- 8*sd(ACE_xrf$cps)
cps.min.thres <- cps.mean - cps.sd 
cps.max.thres <- cps.mean + cps.sd 
cps.min.thres
cps.max.thres

#  OR use this 

# cps tolerance filter 
# cps.min.thres <- 40000
# cps.max.thres <- 90000

Fig3.1B <-ACE_xrf  %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  ggplot(mapping = aes(x = depth, y = cps, col = in_cps_tolerance)) + 
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = c(cps.min.thres, cps.max.thres)) +
  geom_rug(sides = "b", data = . %>% filter(in_cps_tolerance == FALSE)) + 
  theme_bw() + 
  theme(legend.position="top") +
  guides(colour = guide_legend(nrow = 1))
print(Fig3.1B)

# MSE tolerance filter set at 6xSD - need to be careful as MSE can indicate different lithologies 

# MSE tolerance - 2 std dev is too strict when MSE values are very similar but all below 2
MSE.mean <- mean(ACE_xrf$MSE)
MSE.sd <- 6*sd(ACE_xrf$MSE)
MSE.thres <- MSE.mean + MSE.sd 
MSE.thres

#  OR

#MSE.thres <- 2 # use this as an established general threshold

Fig3.1C <- ACE_xrf %>%
  mutate(in_mse_tolerance = ifelse(MSE >=MSE.thres, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth,  y = MSE, col = in_mse_tolerance)) +
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = MSE.thres) +
  geom_rug(sides = "b", data = . %>% filter(in_mse_tolerance == FALSE)) + 
  theme_bw() + 
  theme(legend.position="top") +
  guides(colour = guide_legend(nrow = 1))
print(Fig3.1C)

# Surface slope tolerance filter ------------------------------------------

# Either use this with wide margins for peat cores eg 6 std dev equivalent to +/-0.5 
#slope.min.thres = -0.5
#slope.max.thres = 0.5

#  OR based on mean and SD thresholds 

slope1 <-  ACE_xrf$surface - lag(ACE_xrf$surface)
s1 <- as_tibble(slope1) %>% 
  filter(!if_any(everything(), is.na))
slope.mean <- mean(s1$value)
slope.sd <- 2*sd(s1$value)
slope.min.thres <- slope.mean - slope.sd 
slope.max.thres <- slope.mean + slope.sd 
slope.min.thres
slope.max.thres

Fig3.1D <- ACE_xrf %>%
  mutate(slope = surface - dplyr::lag(surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth, y = slope, col = in_slope_tolerance)) +
  scale_y_continuous(limits = c(-1, 1), oob = scales::squish) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = c(slope.min.thres, slope.max.thres)) +
  geom_rug(data = . %>% filter(validity == FALSE)) +
  scale_x_reverse() +
  theme_bw()  + 
  theme(legend.position="top") +
  guides(colour = guide_legend(nrow = 1))
print(Fig3.1D)

# plotFigure 3.1-3.4:
ggarrange(Fig3.1A, Fig3.1B, Fig3.1C, Fig3.1D, nrow = 2, ncol = 2, labels = c('A', 'B', 'C', 'D'))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.1_QC.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Combining all 'validity' flags  ---------------------------------------------

ACE_xrf_qc_cps <- ACE_xrf %>%
  mutate(slope = surface - dplyr::lag(surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=MSE.thres, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) %>% 
  relocate(qc, .after = validity) %>% 
  filter(qc == TRUE) #to remove from ACE_xrf rows that dont pass QC

# Save QC xrf dataset to use and stats to file
write.csv(ACE_xrf_qc_cps,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/ACE_xrf_qc_cps.csv", row.names = FALSE)

# Example plots ----------------------------------------------------------

# Use Site == to choose example site as an example - here BI10
ACE_xrf_qc_BI10 <- ACE_xrf_qc_cps %>% 
  filter(Site == "BI10")
theme_set(theme_bw(12))
Fig3.2A <- ggplot(data = ACE_xrf_qc_BI10, aes(y = depth, x = Ti, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.2B <- ggplot(data = ACE_xrf_qc_BI10, aes(y = depth, x = Fe, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.2C <- ggplot(data = ACE_xrf_qc_BI10, aes(y = depth, x = Br, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.2D <- ggplot(data = ACE_xrf_qc_BI10, aes(y = depth, x = coh_inc, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

ggarrange(Fig3.2A, Fig3.2B, Fig3.2C, Fig3.2D, ncol = 2, nrow = 2, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.2_ACE_xrf_qc_BI10.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# SECTION 3B: Select dataset to take forward ------------------------------------------

# load qc dataset from above back in
ACE_xrf_qc <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/ACE_xrf_qc_cps.csv")

# Create list of elements - can remove REE and/or Zr (matrix effect in peat/organic seds) at this stage 
ACE_xrf_qc_elements <- select(ACE_xrf_qc, c(Mg:Mo_coh), -c(`Fe a*2`)) %>% #Zr
  names()
ACE_xrf_qc_elements

# SECTION 3C: SNR element filtering ---------------------------------

# Summary stats for xrf_qc dataset
ACE_xrf_qc_stats <- ACE_xrf_qc %>%
  select(kcps:MSE, any_of(ACE_xrf_qc_elements)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(ACE_xrf_qc_stats)

# Signal to noise ratio SNR =𝜇/𝜎 & mean cps >100 and max cps >100
# define variable to remove
remove_var <- c("MSE", "kcps", "cps") # remove unneeded variables in table 
ACE_SNR_stats <- ACE_xrf_qc_stats %>%
  mutate(SNR = mean/sd) %>% 
  #filter(!(element%in%remove_var)) %>% # remove non-element/scatter variables
  filter(SNR>1) #set SNR and cps conditions - can add cps as mean>100 and max>100
ACE_SNR_stats

ACE_SNR_List <-ACE_SNR_stats[,1,2]
ACE_SNR_List 

# Save QC xrf stats to file
write.csv(ACE_xrf_qc_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/ACE_xrf_qc_stats.csv", row.names = FALSE)
write.csv(ACE_SNR_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/ACE_xrf_qc_SNR_stats.csv", row.names = FALSE)

# SECTION 3D: ACF element filtering ------------------------------------------

ACE_xrf_qc <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/ACE_xrf_qc.csv")

# Autocorrelation based filtering of elements -----------------------------

# Use autocorrelation function (acf) and plots to explore noise in a time-series
library(forecast)
library(ggrepel)
library(directlabels)

# Adjust for any element of interest by changing $
# split into two groups below to visualise the most common elements measured by ITRAX
theme_set(theme_bw(8))
Fig3.3 <- ggarrange(
  ggAcf(ACE_xrf_qc$Al) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Si) + ylim(c(NA,1)), 
  ggAcf(ACE_xrf_qc$P) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$S) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Cl) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$K) + ylim(c(NA,1)), 
  ggAcf(ACE_xrf_qc$Ca) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Ti) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$V) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Cr) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Mn) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Fe) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.3)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.3_ACF_pt1.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

Fig3.4 <- ggarrange(
  ggAcf(ACE_xrf_qc$Co) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Ni) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Cu) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Zn) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Se) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Br) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Rb) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Sr) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Zr) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Ba) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Mo_inc) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Mo_coh) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.4)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.4_ACF_pt2.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Filter elements based on acf lag thresholds

# set lag threshold to 20 for whole dataset
# for 1 mm dataset, 20 checks autocorrelation +/-20 mm around each measurement 
# can be performed on different units later
# 0.2 or 0.1 is  minimum threshold
# 0.5 is  maximum threshold - i.e, lag time to half correlation coefficient 

# define filter and lag thresholds
acf_thres_min <- 0.1
acf_thres_max <- 0.5
lag_thres <- 20

# ACF threshold element filtering 
# below uses acf min > 0.1 - can change this to max
Fig3.5a <- apply(ACE_xrf_qc %>% select(any_of(ACE_elements), Mo_inc, Mo_coh), 
                 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", linewidth=0.5, linetype = 3) +
  xlim(0, 50) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig3.5a)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.5a_ACF_all_elements.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# apply acf based filtering to ACE_xrf_qc elmement list - leaving acf elements >0.1 min threshold
apply(ACE_xrf_qc %>% select(any_of(ACE_elements), Mo_inc, Mo_coh), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>% 
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>% 
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% 
                             pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == lag_thres) %>%
  filter(value >= acf_thres_min) %>% #OR acf_thres_max
  pull(elements) %>% 
  ordered() -> acfElements_min # OR acf_thres_max
acfElements_min
#acfElements_max

acfElementsList <- select(ACE_xrf_qc, any_of(acfElements_min)) %>% #OR acfElements_max
  names()
acfElementsList

# Replot with min acf filtered elements only
Fig3.5b <- apply(ACE_xrf_qc %>% select(any_of(acfElements_min)), #OR acfElements_max
                 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", linewidth=0.5, linetype = 3) +
  xlim(0, 50) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig3.5b)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.5b_ACF_elements.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# Summary ACF element plot vs depth
Fig3.6 <- ACE_xrf_qc %>% 
  filter(qc == TRUE) %>% 
  select(any_of(acfElements_min), depth, label) %>% #OR acfElements_max
  pivot_longer(any_of(acfElements_min), names_to = "param", values_to = "element") %>% #OR acfElements_max
  filter(param %in% acfElements_min) %>%
  mutate(param = fct_relevel(param, acfElementsList)) %>%
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
  ggtitle("ACE Composite - all sites: cps, elements ACF >0.1")
print(Fig3.6)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.6_All_Sites_ACFmin.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig3.6, # First row
          ggarrange(Fig3.5a, Fig3.5b, ncol = 2, labels = c("B", "C")), # Second row with two plots
          nrow = 2, 
          labels = "A", common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.7_All_Sites_ACF_min.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# check acf vs SNR element lists
acfElementsList
ACE_SNR_List

# combine lists and remove duplicates
ACE_filterList <- unique(c(acfElementsList, ACE_SNR_List))
ACE_filterList

# SECTION 4: Smoothing -------------------------------------------

# Run these plots with %cps sum data to make sure everything looks OK 
library(tidypaleo)
theme_set(theme_bw(base_size=8))

# BI10
BI10_xrf <- ACE_xrf_qc %>% # can change this to use pc_cps file, inc, or coh/inc normlaised file in here  
  filter(Site == "BI10") %>% 
  filter(qc == TRUE) %>%
  select(all_of(acfElementsList_min), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig3.8 <- full_join(y = BI10_xrf %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acfElementsList_min)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = BI10_xrf %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acfElements_min), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acfElementsList_min)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE BI10: cps, acf elements ")
print(Fig3.8)

# DR171B
DR171B_combined <- c("DRPB_1B.1", "DRPB_1B.2")
DR171B_xrf <- ACE_xrf_qc %>%
  filter(label == DR171B_combined) %>% 
  select(all_of(acfElementsList_min), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig3.9 <- full_join(y = DR171B_xrf %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acfElementsList_min)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = DR171B_xrf  %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acfElementsList_min), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acfElementsList_min)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE DR171B: cps, acf elements ")
print(Fig3.9)

# KER1 - CAT1-S1
CAT1.1_xrf <- ACE_xrf_qc %>%
  filter(Site == "KER1") %>% 
  select(all_of(acfElementsList_min), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig3.10 <- full_join(y = CAT1.1_xrf %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acfElementsList_min)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = CAT1.1_xrf %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acfElementsList_min), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acfElementsList_min)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT1-S1 (KER1): cps, acf elements ")
print(Fig3.10)

# KER3 - CAT2-S3
CAT2.3_xrf <- ACE_xrf_qc %>%
  filter(Site == "KER3") %>% 
  select(all_of(acfElementsList_min), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig3.11 <- full_join(y = CAT2.3_xrf %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acfElementsList_min)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = CAT2.3_xrf %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acfElementsList_min), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acfElementsList_min)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT2-S3 (KER3): cps, acf elements ")
print(Fig3.11)

# HER42PB
H42PB_xrf <- ACE_xrf_qc %>%
  filter(Site == "HER42PB") %>% 
  select(all_of(acfElementsList_min), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig3.12 <- full_join(y = H42PB_xrf %>%
            as_tibble() %>%
            # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
            mutate(across(any_of(c(acfElementsList_min)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = H42PB_xrf %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acfElementsList_min), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acfElementsList_min)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
glimpse() %>%
    ggplot(aes(x = peakarea, y = depth)) +
    tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
    scale_alpha_manual(values = c(0.1, 1)) +
    theme(legend.position="bottom") +
    guides(colour = guide_legend(nrow = 2)) +
    scale_y_reverse() +
    scale_x_continuous(n.breaks = 4) +
    facet_geochem_gridh(vars(elements)) +
    labs(x = "Peak area [cps]", y = "Depth [mm]")+
    ggtitle("ACE HER42PB: cps, acf elements")
print(Fig3.12)

# PB1
PB1_xrf <- ACE_xrf_qc %>%
  filter(Site == "PB1") %>% 
  select(all_of(acfElementsList_min), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig3.13 <- full_join(y = PB1_xrf  %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acfElementsList_min)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = PB1_xrf  %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acfElementsList_min), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acfElementsList_min)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE PB1: cps, acf elements ")
print(Fig3.13)

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig3.8, Fig3.9, nrow = 2, labels = c("A", "B"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section4/Figures/Fig4_A-B.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig3.10, Fig3.11, nrow = 2, labels = c("C", "D"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section4/Figures/Fig4_C-D.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig3.12, Fig3.13, nrow = 2, labels = c("E", "F"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig4_E-F.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# SECTION 5: ACE QC and element filtered ------------------------

# Import orginal cps dataset & element lists
# Set working directory 
ACE_xrf_qc <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section2/ACE_xrf_cps.csv") 

# Create list of elements for acf - also remove Zr which is matrix effect in peat/organic seds at this stage 
ACE_elements <- select(ACE_xrf_qc, c(Mg:Mo_coh), -c(all_of(machine_elements), `Fe a*2`)) %>% # , Zr
  names()
ACE_elements

ACE_elements1 <- select(ACE_xrf_qc, c(Mg:Mo_coh, Total_scatter, inc_coh, coh_inc), -c(all_of(machine_elements), `Fe a*2`)) %>% # , Zr
  names()
ACE_elements1

# Create SNR element list
ACE_xrf_qc_SNR_stats <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/ACE_xrf_qc_SNR_stats.csv")
ACE_xrf_qc_SNR_List <-ACE_xrf_qc_SNR_stats[,1,2]
ACE_xrf_qc_SNR_List 

# selected by itrax.R QC process for HER42PB and PB1, when run separately
acfElements_key <- c("P",  "K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", 
                  "Cu", "Zn", "Rb", "Sr", "Zr", "Pb")
icp_key_Elements <- c("P_31", "K_39", "Ca_43", "Ti_49", "Mn_55", 
                      "Fe_56", "Co_59", "Ni_60", "Cu_65", "Zn_66", 
                      "Rb_85", "Sr_86", "Zr_90", "Pb_206")

# selected by itrax.R QC process for all data without POB4 run together - then matched to icp list above
acfElements_min <- c("Si", "P", "S", "Cl", "K", "Ca", "Ti", "V",
                     "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Se", 
                     "Br", "Rb", "Sr", "Zr", "Ba", "Pb", "Mo_inc", "Mo_coh")
acfElements_max <-  c("S", "Cl", "K", "Ca", "Ti", "V", "Mn", "Fe", 
                      "Zn", "Br", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")

# icp and xrf element matching
acf_icp_Elements_max <- c("K", "Ca", "Ti", "Mn", "Fe", "Zn", 
                      "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh", 
                      "inc_coh", "coh_inc")
icp_Elements_max <- c("K_39", "Ca_43", "Ti_49", "Mn_55", "Fe_56", 
                  "Zn_66", "Rb_85", "Sr_86", "Zr_90")

# ICPMS
# ICPMS elements  defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")
# ICPMS elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# XRF and ICPMS elements combined
xrf_icp_elements <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                      "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                      "Zr", "Zr_ICP", "Mo_coh")
# Same list as above but with Mo_inc inlcuded
xrf_icp_elements1 <- c("K", "K_ICP", "Ca", "Ca_ICP", "Ti", "Ti_ICP", "Mn", "Mn_ICP",
                       "Fe", "Fe_ICP", "Zn", "Zn_ICP", "Rb", "Rb_ICP", "Sr", "Sr_ICP",
                       "Zr", "Zr_ICP", "Mo_inc",  "Mo_coh")

# we need some kind of justification why these were chosen
# eg I'd like to include Br even if it's not considered reliable/mobile etc., and to see if the calibration fails 

# ACF min elements & QC filtering based on multiple QC parameters in Section 4 ----------------------

# Import the qc and element filtered dataset
ACE_xrf_qc <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/ACE_xrf_qc.csv") 
ACE_xrf_qc 

ACE_xrf_qc_acf <- ACE_xrf_qc %>%
  filter(qc == TRUE) %>% #to remove from ACE_xrf rows that dont pass QC
  select(Location:MSE, all_of(acf_icp_Elements_max), Total_scatter:coh_inc, qc)
ACE_xrf_qc_acf

# Summary stats for ACE_xrf_qc
ACE_xrf_qc_acf_stats <- ACE_xrf_qc_acf %>%
  select(kcps:MSE, any_of(ACE_elements), Total_scatter:coh_inc) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(ACE_xrf_qc_acf_stats)

# Save QC xrf dataset to use and stats to file
write.csv(ACE_xrf_qc_acf,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section5/ACE_xrf_qc_acf.csv", row.names = FALSE)
write.csv(ACE_xrf_qc_acf_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section5/ACE_xrf_qc_acf_stats.csv", row.names = FALSE)

# SECTION 6A: Create %cps, /inc, clr datasets for calibration exercise -----------------

# Import the qc and element filtered dataset ------------------------------
ACE_xrf_qc_acf <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section5/ACE_xrf_qc_acf.csv") 
ACE_xrf_qc_acf

# Transform QC cps dataframe to %cps_sum, inc normalised, and CIR  --------
# cps as % of cps_sum
ACE_xrf_qc_acf_pc <- ACE_xrf_qc_acf %>%
  mutate(across(all_of(acf_icp_Elements_max)) / `cps`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  replace(is.na(.), 0) %>%
  mutate(across(all_of(acf_icp_Elements_max)) *100) %>%
  mutate(cps_sum2 = rowSums(across(all_of(acf_icp_Elements_max)))) %>% 
  print()
# check sum <100% & write file
head(ACE_xrf_qc_acf_pc$cps_sum2)
write.csv(ACE_xrf_qc_acf_pc,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_pc.csv", row.names = FALSE)

# %cps_sum summary stats table
ACE_xrf_qc_acf_pc_stats <- ACE_xrf_qc_acf_pc %>%
  select(any_of(acf_icp_Elements_max), Total_scatter:cps_sum2) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_xrf_qc_acf_pc_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_pc_stats.csv", row.names = FALSE)

# Normalise by inc scatter
ACE_xrf_qc_acf_inc <- ACE_xrf_qc_acf %>%
  #select(any_of(ACE_elements)) %>% 
  mutate(across(any_of(ACE_elements)) /`Mo_inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf)))
ACE_xrf_qc_acf_inc
# check inc = 1 & write file
tail(ACE_xrf_qc_acf_inc$Mo_inc)
write.csv(ACE_xrf_qc_acf_inc,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_inc.csv", row.names = FALSE)

# Standardise and centre inc normalised dataframe
ACE_xrf_qc_acf_inc_Z <- ACE_xrf_qc_acf_inc %>% 
  select(any_of(acfElements_key))
ACE_xrf_qc_acf_inc_Z[, acfElements_key] <- scale(ACE_xrf_qc_acf_inc_Z[, acfElements_min], center = TRUE, scale = TRUE)
ACE_xrf_qc_acf_inc_Z

#select columns & remake original table
ACE_xrf_qc_acf_titles <- ACE_xrf_qc_acf_inc %>% 
  select(Location:MSE)
ACE_xrf_qc_acf_scatter <- ACE_xrf_qc_acf_inc %>% 
  select(Total_scatter, inc_coh, coh_inc)
ACE_xrf_qc_acf_inc_Z1 <- bind_cols(ACE_xrf_qc_acf_titles, ACE_xrf_qc_acf_inc_Z, ACE_xrf_qc_acf_scatter)
ACE_xrf_qc_acf_inc_Z1
write.csv(ACE_xrf_qc_acf_inc_Z1,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_inc_Z.csv", row.names = FALSE)

# Normalised by coh/inc ratio (Boyle et al., 2015 method)
ACE_xrf_qc_acf_CIR <- ACE_xrf_qc_acf %>%
  #select(any_of(ACE_elements)) %>% 
  mutate(across(any_of(acfElements_min)) / coh_inc) %>%
  mutate_if(is.numeric, list(~na_if(., Inf)))
ACE_xrf_qc_acf_CIR
write.csv(ACE_xrf_qc_acf_CIR,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_CIR.csv", row.names = FALSE)

#set NA back to 0 for QC cps dataframe for clr and next section 
ACE_xrf_qc_acf <- ACE_xrf_qc_acf %>% 
  replace(is.na(.), 0) 
ACE_xrf_qc_acf

# clr (centred log ratios)
# library(compositions)
# # Cannot run clr with zeroes or NAs - run this with key elements only to match icp later on
# # removing rows with zeros creates too many gaps in the Site for peat cores 
# # so +1 to all cps
# ACE_clr0 <- ACE_xrf_qc_acf %>%
#   replace(is.na(.), 0) %>%
#   select(any_of(acf_icp_Elements_min)) %>% 
#   select(-c("Mo_inc",  "Mo_coh",  "inc_coh", "coh_inc"))
# ACE_clr1 <- ACE_clr0 + 1 #can run without +1 but this reduces the number of elements to 8
# as_tibble(ACE_clr1)
# # Or replace all 0 with NA and run clr
# #ACE_clr1[ACE_clr1 == 0] <- NA
# ACE_clr1 <- as_tibble(ACE_clr1) %>%
#   clr() 
# head(ACE_clr1)
# tail(ACE_clr1)
# # Add scatter columns back into the clr dataframe for plotting in next section
# ACE_xrf_qc_acf_scatter <- ACE_xrf_qc_acf %>% 
#   select(Mo_inc, Mo_coh, Total_scatter, inc_coh, coh_inc)
# ACE_xrf_qc_acf_clr <-  bind_cols(ACE_xrf_qc_acf_titles, ACE_clr1, ACE_xrf_qc_acf_scatter) %>%
# write.csv(ACE_xrf_qc_acf_clr,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_clr.csv", row.names = FALSE)
# ACE_xrf_qc_acf_clr

# clr (centered log ratios) dataset -------------------------------------------------------------------------

# using PLK1_elements1, without Mo_inc or Mo_coh 

library(compositions)

# Cannot run clr with zeroes or NAs - removing rows with zeros creates too many gaps in the record for peat cores 
## Replace zeros in each data column with half minimum value of that column to allow linear modelling to work
# Recommended procedure from Bertrand et al. (2023) - retains dataframe structure
ACE_clr0 <- ACE_xrf_qc_acf %>%
  replace(is.na(.), 0) %>%
  select(any_of(acf_icp_Elements_max)) %>%
  select(-c(Mo_inc, Mo_coh, inc_coh, coh_inc)) %>% 
ACE_clr1 <- ACE_clr0 %>% 
  mutate_at(vars(any_of(acf_icp_Elements_max)), 
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% 
  print()

# Or replace all 0 with NA and run clr
#ACE_clr1[ACE_clr1 == 0] <- NA

ACE_clr1 <- as_tibble(ACE_clr1) %>%
  clr() 
head(ACE_clr1)
tail(ACE_clr1)

# Add scatter columns back into the clr dataframe for plotting in next section
ACE_xrf_qc_acf_scatter <- ACE_xrf_qc_acf %>% 
  select(Mo_inc, Mo_coh, Total_scatter, inc_coh, coh_inc)
ACE_xrf_qc_acf_clr <-  bind_cols(ACE_xrf_qc_acf_titles, ACE_clr1, ACE_xrf_qc_acf_scatter)
write.csv(ACE_xrf_qc_acf_clr,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_clr.csv", row.names = FALSE)
ACE_xrf_qc_acf_clr


# SECTION 6B: ACE - individual Sites plotted with smoothing---------------

# Import the cps qc and element filtered dataset
ACE_xrf_qc_acf_cps <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section5/ACE_xrf_qc_acf.csv") 
ACE_xrf_qc_acf_cps
library(tidypaleo)
theme_set(theme_bw(base_size=8))

# BI10
BI10_xrf2_cps <- ACE_xrf_qc_acf_cps %>% # can change this to use pc_cps file, inc, or coh/inc normlaised file in here  
  filter(Site == "BI10") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.1 <- full_join(y = BI10_xrf2_cps %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = BI10_xrf2_cps %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE BI10: cps, qc & acf elements ")
print(Fig6.1)

# DR171B
DR171B_combined <- c("DRPB_1B.1", "DRPB_1B.2")
DR171B_xrf2_cps <- ACE_xrf_qc_acf_cps %>%
  filter(label == DR171B_combined) %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.2 <- full_join(y = DR171B_xrf2_cps %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = DR171B_xrf2_cps  %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE DR171B: cps, qc & acf elements ")
#print(Fig6.2)

# KER1 - CAT1-S1
CAT1.1_xrf2_cps <- ACE_xrf_qc_acf_cps %>%
  filter(Site == "KER1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.3 <- full_join(y = CAT1.1_xrf2_cps %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = CAT1.1_xrf2_cps %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT1-S1 (KER1): cps, qc & acf elements ")
#print(Fig6.3)

# KER3 - CAT2-S3
CAT2.3_xrf2_cps <- ACE_xrf_qc_acf_cps %>%
  filter(Site == "KER3") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.4 <- full_join(y = CAT2.3_xrf2_cps %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = CAT2.3_xrf2_cps %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT2-S3 (KER3): cps, qc & acf elements ")
#print(Fig6.4)

# HER42PB
H42PB_xrf2_cps <- ACE_xrf_qc_acf_cps %>%
  filter(Site == "HER42PB") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.5 <- full_join(y = H42PB_xrf2_cps %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = H42PB_xrf2_cps %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE HER42PB: cps, qc & acf elements ")
#print(Fig6.5)

# PB1
PB1_xrf2_cps <- ACE_xrf_qc_acf_cps %>%
  filter(Site == "PB1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.6 <- full_join(y = PB1_xrf2_cps  %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = PB1_xrf2_cps  %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]")+
  ggtitle("ACE PB1: cps, qc & acf elements ")
#print(Fig6.6)

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.1, Fig6.2, nrow = 2, labels = c("A", "B"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.1-6.2_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.3, Fig6.4, nrow = 2, labels = c("C", "D"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.3-6.4_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.5, Fig6.6, nrow = 2, labels = c("E", "F"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.5-6.6_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Import the %cps qc and element filtered dataset --------------------------
ACE_xrf_qc_acf_pc <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_pc.csv") 
ACE_xrf_qc_acf_pc
library(tidypaleo)
theme_set(theme_bw(base_size=8))

# BI10
BI10_xrf2_pc <- ACE_xrf_qc_acf_pc %>% # can change this to use pc_cps file, inc, or coh/inc normlaised file in here  
  filter(Site == "BI10") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.1 <- full_join(y = BI10_xrf2_pc %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = BI10_xrf2_pc %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [%cps]", y = "Depth [mm]")+
  ggtitle("ACE BI10: %cps, qc & acf elements ")
#print(Fig6.1)

# DR171B
DR171B_combined <- c("DRPB_1B.1", "DRPB_1B.2")
DR171B_xrf2_pc <- ACE_xrf_qc_acf_pc %>%
  filter(label == DR171B_combined) %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.2 <- full_join(y = DR171B_xrf2_pc %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = DR171B_xrf2_pc  %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [%cps]", y = "Depth [mm]")+
  ggtitle("ACE DR171B: %cps, qc & acf elements ")
#print(Fig6.2)

# KER1 - CAT1-S1
CAT1.1_xrf2_pc <- ACE_xrf_qc_acf_pc %>%
  filter(Site == "KER1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.3 <- full_join(y = CAT1.1_xrf2_pc %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = CAT1.1_xrf2_pc %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [%cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT1-S1 (KER1): %cps, qc & acf elements ")
#print(Fig6.3)

# KER3 - CAT2-S3
CAT2.3_xrf2_pc <- ACE_xrf_qc_acf_pc %>%
  filter(Site == "KER3") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.4 <- full_join(y = CAT2.3_xrf2_pc %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = CAT2.3_xrf2_pc %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [%cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT2-S3 (KER3): %cps, qc & acf elements ")
#print(Fig6.4)

# HER42PB
H42PB_xrf2_pc <- ACE_xrf_qc_acf_pc %>%
  filter(Site == "HER42PB") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.5 <- full_join(y = H42PB_xrf2_pc %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = H42PB_xrf2_pc %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [%cps]", y = "Depth [mm]")+
  ggtitle("ACE HER42PB: %cps, qc & acf elements ")
#print(Fig6.5)

# PB1
PB1_xrf2_pc <- ACE_xrf_qc_acf_pc %>%
  filter(Site == "PB1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.6 <- full_join(y = PB1_xrf2_pc  %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = PB1_xrf2_pc  %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [%cps]", y = "Depth [mm]")+
  ggtitle("ACE PB1: %cps, qc & acf elements ")
#print(Fig6.6)

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.1, Fig6.2, nrow = 2, labels = c("A", "B"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.1-6.2_pc.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.3, Fig6.4, nrow = 2, labels = c("C", "D"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.3-6.4_pc.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.5, Fig6.6, nrow = 2, labels = c("E", "F"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.5-6.6_pc.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Import the CIR qc and element filtered dataset --------------------------
ACE_xrf_qc_acf_CIR <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_CIR.csv") 
ACE_xrf_qc_acf_CIR
library(tidypaleo)
theme_set(theme_bw(base_size=8))

# BI10
BI10_xrf2_CIR <- ACE_xrf_qc_acf_CIR %>% # can change this to use pc_cps file, inc, or coh/inc normlaised file in here  
  filter(Site == "BI10") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.1 <- full_join(y = BI10_xrf2_CIR %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = BI10_xrf2_CIR %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [CIR]", y = "Depth [mm]")+
  ggtitle("ACE BI10: CIR, qc & acf elements ")
#print(Fig6.1)

# DR171B
DR171B_combined <- c("DRPB_1B.1", "DRPB_1B.2")
DR171B_xrf2_CIR <- ACE_xrf_qc_acf_CIR %>%
  filter(label == DR171B_combined) %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.2 <- full_join(y = DR171B_xrf2_CIR %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = DR171B_xrf2_CIR  %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [CIR]", y = "Depth [mm]")+
  ggtitle("ACE DR171B: CIR, qc & acf elements ")
#print(Fig6.2)

# KER1 - CAT1-S1
CAT1.1_xrf2_CIR <- ACE_xrf_qc_acf_CIR %>%
  filter(Site == "KER1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.3 <- full_join(y = CAT1.1_xrf2_CIR %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = CAT1.1_xrf2_CIR %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [CIR]", y = "Depth [mm]")+
  ggtitle("ACE CAT1-S1 (KER1): CIR, qc & acf elements ")
#print(Fig6.3)

# KER3 - CAT2-S3
CAT2.3_xrf2_CIR <- ACE_xrf_qc_acf_CIR %>%
  filter(Site == "KER3") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.4 <- full_join(y = CAT2.3_xrf2_CIR %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = CAT2.3_xrf2_CIR %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [CIR]", y = "Depth [mm]")+
  ggtitle("ACE CAT2-S3 (KER3): CIR, qc & acf elements ")
#print(Fig6.4)

# HER42PB
H42PB_xrf2_CIR <- ACE_xrf_qc_acf_CIR %>%
  filter(Site == "HER42PB") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.5 <- full_join(y = H42PB_xrf2_CIR %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = H42PB_xrf2_CIR %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [CIR]", y = "Depth [mm]")+
  ggtitle("ACE HER42PB: CIR, qc & acf elements ")
#print(Fig6.5)

# PB1
PB1_xrf2_CIR <- ACE_xrf_qc_acf_CIR %>%
  filter(Site == "PB1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)
# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.6 <- full_join(y = PB1_xrf2_CIR  %>%
                      as_tibble() %>%
                      # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                      mutate(across(any_of(c(acf_icp_Elements_max)), 
                                    function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                      )
                      ) %>%
                      mutate(type = "mean"), 
                    x = PB1_xrf2_CIR  %>% 
                      as_tibble() %>% 
                      mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area [CIR]", y = "Depth [mm]")+
  ggtitle("ACE PB1: CIR, qc & acf elements ")
#print(Fig6.6)

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.1, Fig6.2, nrow = 2, labels = c("A", "B"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.1-6.2_CIR.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.3, Fig6.4, nrow = 2, labels = c("C", "D"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.3-6.4_CIR.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.5, Fig6.6, nrow = 2, labels = c("E", "F"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.5-6.6_CIR.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Import the clr qc and element filtered dataset --------------------------
ACE_xrf_qc_acf_clr <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/ACE_xrf_qc_acf_clr.csv") 
ACE_xrf_qc_acf_clr

library(tidypaleo)
theme_set(theme_bw(base_size=8))

# BI10
BI10_xrf2_clr <- ACE_xrf_qc_acf_clr %>% # can change this to use pc_cps file, inc, or coh/inc normlaised file in here  
  filter(Site == "BI10") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.1 <- full_join(y = BI10_xrf2_clr %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acf_icp_Elements_max)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = BI10_xrf2_clr %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area: Elements [clr] Scatter [cps]", y = "Depth [mm]")+
  ggtitle("ACE BI10: clr, qc & acf elements ")
#print(Fig6.1)

# DR171B
DR171B_combined <- c("DRPB_1B.1", "DRPB_1B.2")
DR171B_xrf2_clr <- ACE_xrf_qc_acf_clr %>%
  filter(label == DR171B_combined) %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.2 <- full_join(y = DR171B_xrf2_clr %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acf_icp_Elements_max)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = DR171B_xrf2_clr  %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area: Elements [clr] Scatter [cps]", y = "Depth [mm]")+
  ggtitle("ACE DR171B: clr, qc & acf elements ")
#print(Fig6.2)

# KER1 - CAT1-S1
CAT1.1_xrf2_clr <- ACE_xrf_qc_acf_clr %>%
  filter(Site == "KER1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.3 <- full_join(y = CAT1.1_xrf2_clr %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acf_icp_Elements_max)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = CAT1.1_xrf2_clr %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area: Elements [clr] Scatter [cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT1-S1 (KER1): qc & clr, acf elements ")
#print(Fig6.3)

# KER3 - CAT2-S3
CAT2.3_xrf2_clr <- ACE_xrf_qc_acf_clr %>%
  filter(Site == "KER3") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.4 <- full_join(y = CAT2.3_xrf2_clr %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acf_icp_Elements_max)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = CAT2.3_xrf2_clr %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area: Elements [clr] Scatter [cps]", y = "Depth [mm]")+
  ggtitle("ACE CAT2-S3 (KER3): clr, qc & acf elements ")
#print(Fig6.4)

# HER42PB
H42PB_xrf2_clr <- ACE_xrf_qc_acf_clr %>%
  filter(Site == "HER42PB") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.5 <- full_join(y = H42PB_xrf2_clr %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acf_icp_Elements_max)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = H42PB_xrf2_clr %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area: Elements [clr] Scatter [cps]", y = "Depth [mm]")+
  ggtitle("ACE HER42PB: clr, qc & acf elements ")
#print(Fig6.5)

# PB1
PB1_xrf2_clr <- ACE_xrf_qc_acf_clr %>%
  filter(Site == "PB1") %>% 
  select(all_of(acf_icp_Elements_max), MSE,  depth, label)

# plotting the running mean stratigraphically - and add onto XRFplot for each Site
# with ACF elements only
Fig6.6 <- full_join(y = PB1_xrf2_clr  %>%
                       as_tibble() %>%
                       # uses a 10 point running mean (10 mm for this data); 5 before, 5 after
                       mutate(across(any_of(c(acf_icp_Elements_max)), 
                                     function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
                       )
                       ) %>%
                       mutate(type = "mean"), 
                     x = PB1_xrf2_clr  %>% 
                       as_tibble() %>% 
                       mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>%
  select(all_of(acf_icp_Elements_max), MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", all_of(acf_icp_Elements_max)))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  glimpse() %>%
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "Peak area: Elements [clr] Scatter [cps]", y = "Depth [mm]")+
  ggtitle("ACE PB1: clr, qc & acf elements ")
#print(Fig6.6)

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.1, Fig6.2, nrow = 2, labels = c("A", "B"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.1-6.2_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.3, Fig6.4, nrow = 2, labels = c("C", "D"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.3-6.4_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig6.5, Fig6.6, nrow = 2, labels = c("E", "F"), common.legend = FALSE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section6/Figures/Fig6.5-6.6_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# END








# REMOVED CODE - move to separate file before publication ------------------------------------------------------------

# Old transformation code for cps dataset prior to qc filtering - not really needed - 3/6/23
# ---------------------------------------------------------------------------------------------------------
# SECTION 2: Transform VALID cps dataframe to %cps_sum, inc normalised, and CIR (coh/inc) normalised
# ---------------------------------------------------------------------------------------------------------

# Transformation before qc filtering --------------------------------------

# cps as % of cps_sum

# original cps qc dataframe
ACE_xrf_pc <- ACE_xrf %>%
  mutate(across(Mg:Mo_coh) / `cps`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  replace(is.na(.), 0) %>%
  mutate(across(Mg:Mo_coh) *100) %>%
  select(-c(`Fe a*2`)) %>% 
  mutate(cps_sum = rowSums(across(Mg:Mo_coh))) %>% 
  print()
# check sum = 100% & write file
head(ACE_xrf_pc$cps_sum)
ACE_xrf_pc
write.csv(ACE_xrf_pc,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section2/ACE_xrf_pc.csv", row.names = FALSE)

# %cps_sum summary stats table
ACE_xrf_pc_stats <- ACE_xrf_pc %>%
  select(any_of(ACE_elements), Total_scatter:cps_sum) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(ACE_xrf_pc_stats)
write.csv(ACE_xrf_pc_stats,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section2/ACE_xrf_pc_stats.csv", row.names = FALSE)

# Normalise by inc scatter
ACE_xrf_inc <- ACE_xrf %>%
  #select(any_of(ACE_elements)) %>% 
  mutate(across(any_of(ACE_elements)) /`Mo_inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf)))
ACE_xrf_inc
# check inc = 1 & write file
tail(ACE_xrf_inc$Mo_inc)
write.csv(ACE_xrf_inc,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section2/ACE_xrf_inc.csv", row.names = FALSE)

# Standardise and centre inc normalised dataframe
ACE_xrf_inc_Z <- ACE_xrf_inc %>% 
  select(ACE_elements)
ACE_xrf_inc_Z[, ACE_elements] <- scale(ACE_xrf_inc[, ACE_elements], center = TRUE, scale = TRUE)
ACE_xrf_inc_Z

#select columns & remake original table
ACE_xrf_titles <- ACE_xrf_inc %>% 
  select(Location:MSE)
ACE_xrf_scatter <- ACE_xrf_inc %>% 
  select(Total_scatter:coh_inc)
ACE_xrf_inc_Z1 <- bind_cols(ACE_xrf_titles, ACE_xrf_inc_Z, ACE_xrf_scatter)
ACE_xrf_inc_Z1
write.csv(ACE_xrf_inc_Z1,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section2/ACE_xrf_inc_Z.csv", row.names = FALSE)

# Normalised by inc scatter
ACE_xrf_CIR <- ACE_xrf %>%
  #select(any_of(ACE_elements)) %>% 
  mutate(across(any_of(ACE_elements)) / coh_inc) %>%
  mutate_if(is.numeric, list(~na_if(., Inf)))
ACE_xrf_CIR
write.csv(ACE_xrf_CIR,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section2/ACE_xrf_CIR.csv", row.names = FALSE)

#set NA back to 0 for _xrf dataframe for next section 
ACE_xrf <- ACE_xrf %>% 
  replace(is.na(.), 0) 
ACE_xrf

# Removed code 10/4/2024 - related to inlcusion of POB4

# POB4 only -------------------------------------------------------------------------

# Data filtering criteria and quality control measures follow itrax.R principles
# can skip this section if happy to accept validity = TRUE from ITRAX measurements  
ACE_xrf_POB4 <- filter(ACE_xrf, Site=='POB4') %>% 
  select(Location:MSE, any_of(ACE_elements), `Fe a*2`, Total_scatter:coh_inc)

# Make ACF element xrf dataset to use and stats to file
ACE_xrf_POB4_acf <- ACE_xrf_POB4 %>% 
  select(Location:MSE, any_of(acfElements_min), Total_scatter:coh_inc)

write.csv(ACE_xrf_POB4,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section4/ACE_xrf_POB4.csv", row.names = FALSE)
write.csv(ACE_xrf_POB4_acf,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section4/ACE_xrf_POB4_acf.csv", row.names = FALSE)

# MSE and cps filtering --------------------------------------------------

# cps filtering - Fe a*2
Fig4.3A <- 
  ggplot(ACE_xrf_POB4, mapping = aes(x = cps, y = `Fe a*2`)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()
print(Fig4.3A)

# cps - 2 std dev is too strict for HER42PB and other peat cores with a combination of low and high count matrices
cps.mean <- mean(ACE_xrf_POB4$cps)
cps.sd <- 6*sd(ACE_xrf_POB4$cps)
cps.min.thres <- cps.mean - cps.sd 
cps.max.thres <- cps.mean + cps.sd 

#  OR use this 

# cps tolerance filter 
# cps.min.thres <- 40000
# cps.max.thres <- 90000

Fig4.3B <-ACE_xrf_POB4  %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  ggplot(mapping = aes(x = depth, y = cps, col = in_cps_tolerance)) + 
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = c(cps.min.thres, cps.max.thres)) +
  geom_rug(sides = "b", data = . %>% filter(in_cps_tolerance == FALSE)) + 
  theme_bw() + 
  theme(legend.position="top") +
  guides(colour = guide_legend(nrow = 1))
print(Fig4.3B)

# MSE tolerance filter set at 6xSD - need to be careful as MSE can indicate different lithologies 

# MSE tolerance - 2 std dev is too strict when MSE values are very similar but all below 2
MSE.mean <- mean(ACE_xrf_POB4$MSE)
MSE.sd <- 6*sd(ACE_xrf_POB4$MSE)
MSE.thres <- MSE.mean + MSE.sd 
MSE.thres

#  OR

#MSE.thres <- 2 # use this as an established general threshold

Fig4.3C <- ACE_xrf_POB4 %>%
  mutate(in_mse_tolerance = ifelse(MSE >=MSE.thres, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth,  y = MSE, col = in_mse_tolerance)) +
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = MSE.thres) +
  geom_rug(sides = "b", data = . %>% filter(in_mse_tolerance == FALSE)) + 
  theme_bw() + 
  theme(legend.position="top") +
  guides(colour = guide_legend(nrow = 1))
print(Fig4.3C)

# Surface slope tolerance filter ------------------------------------------

# Either use this with wide margins for peat cores eg 6 std dev equivalent to +/-0.5 
#slope.min.thres = -0.5
#slope.max.thres = 0.5

#  OR based on mean and SD thresholds 

slope1 <-  ACE_xrf_POB4$surface - lag(ACE_xrf_POB4$surface)
s1 <- as_tibble(slope1) %>% 
  filter(!if_any(everything(), is.na))
slope.mean <- mean(s1$value)
slope.sd <- 4*sd(s1$value)
slope.min.thres <- slope.mean - slope.sd 
slope.max.thres <- slope.mean + slope.sd 

Fig4.3D <- ACE_xrf_POB4 %>%
  mutate(slope = surface - dplyr::lag(surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth, y = slope, col = in_slope_tolerance)) +
  scale_y_continuous(limits = c(-1, 1), oob = scales::squish) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = c(slope.min.thres, slope.max.thres)) +
  geom_rug(data = . %>% filter(validity == FALSE)) +
  scale_x_reverse() +
  theme_bw()  + 
  theme(legend.position="top") +
  guides(colour = guide_legend(nrow = 1))
print(Fig4.3D)

# show Figure 3.1 - 3.4 together:
ggarrange(Fig4.3A, Fig4.3B, Fig4.3C, Fig4.3D, nrow = 2, ncol = 2, labels = c('A', 'B', 'C', 'D'))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section4/Figures/Fig4.3_POB4_QC.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Combining all 'validity' flags  ---------------------------------------------
ACE_xrf_POB4_qc <- ACE_xrf_POB4 %>%
  mutate(slope = surface - dplyr::lag(surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=MSE.thres, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) 
#filter(qc == TRUE) #to remove from ACE_xrf rows that dont pass QC

# Save QC xrf POB4 dataset to use and stats to file
write.csv(ACE_xrf_POB4_qc,"Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section4/ACE_xrf_POB4_qc.csv", row.names = FALSE)

# plot summary
# Use Site == POB4 site as an example 
theme_set(theme_bw(12))
Fig4.4A <- ggplot(data = ACE_xrf_POB4_qc, aes(y = depth, x = Ti, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig4.4B <- ggplot(data = ACE_xrf_POB4_qc, aes(y = depth, x = Fe, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig4.4C <- ggplot(data = ACE_xrf_POB4_qc, aes(y = depth, x = Br, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig4.4D <- ggplot(data = ACE_xrf_POB4_qc, aes(y = depth, x = coh_inc, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

ggarrange(Fig4.4A, Fig4.4B, Fig4.4C, Fig4.4D, ncol = 2, nrow = 2, common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section4/Figures/Fig4.4_POB4_QC_example.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# -------------------------------------------------------------------------
# Element filtering - ALL SITES DATA - run once
# -------------------------------------------------------------------------

ACE_xrf_qc <- read_csv("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section4/ACE_xrf_qc.csv")

# Autocorrelation based filtering of elements -----------------------------

# Use autocorrelation function (acf) and plots to explore noise in a time-series
library(forecast)
library(ggrepel)
library(directlabels)

# Adjust for any element of interest by changing $
# split into two groups below to visualise the most common elements measured by ITRAX
theme_set(theme_bw(8))
Fig3.3 <- ggarrange(
  ggAcf(ACE_xrf_qc$Al) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Si) + ylim(c(NA,1)), 
  ggAcf(ACE_xrf_qc$P) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$S) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Cl) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$K) + ylim(c(NA,1)), 
  ggAcf(ACE_xrf_qc$Ca) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Ti) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$V) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Cr) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Mn) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Fe) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.3)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section3/Figures/Fig3.3_ACF_pt1.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

Fig3.4 <- ggarrange(
  ggAcf(ACE_xrf_qc$Co) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Ni) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Cu) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Zn) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Se) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Br) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Rb) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Sr) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Zr) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Ba) + ylim(c(NA,1)),
  ggAcf(ACE_xrf_qc$Mo_inc) + ylim(c(NA,1)), ggAcf(ACE_xrf_qc$Mo_coh) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.4)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section3/Figures/Fig3.4_ACF_pt2.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Filter elements based on acf lag thresholds

# set lag threshold to 20 for whole dataset
# for 1 mm dataset, 20 checks autocorrelation +/-20 mm around each measurement 
# can be performed on different units later
# 0.2 or 0.1 is  minimum threshold
# 0.5 is  maximum threshold - i.e, lag time to half correlation coefficient 

# define filter and lag thresholds
acf_thres_min <- 0.1
acf_thres_max <- 0.5
lag_thres <- 20

# MINIMUM ACF threshold element filtering ACF > 0.1 - ALL SITES DATA ----------------------------------------

Fig3.5a <- apply(ACE_xrf_qc %>% select(any_of(ACE_elements), Mo_inc, Mo_coh), 
                 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", size=0.5, linetype = 3) +
  xlim(0, 50) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig3.5a)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section3/Figures/Fig3.5a_ACF_all_elements.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# apply acf based filtering to ACE_xrf_qc elmement list - leaving acf elements >0.1 min threshold
apply(ACE_xrf_qc %>% select(any_of(ACE_elements), Mo_inc, Mo_coh), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>% 
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>% 
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% 
                             pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == lag_thres) %>%
  filter(value >= acf_thres_min) %>% 
  pull(elements) %>% 
  ordered() -> acfElements_min 
acfElements_min

acfElementsList_min <- select(ACE_xrf_qc, any_of(acfElements_min)) %>% 
  names()
acfElementsList_min

# Replot with min acf filtered elements only
Fig3.5b <- apply(ACE_xrf_qc %>% select(any_of(acfElements_min)), 
                 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", size=0.5, linetype = 3) +
  xlim(0, 50) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig3.5b)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Section3/Figures/Fig3.5b_ACF_elements.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# Summary ACF element plot vs depth
Fig3.6 <- ACE_xrf_qc %>% 
  filter(qc == TRUE) %>% 
  select(any_of(acfElements_min), depth, label) %>%
  pivot_longer(any_of(acfElements_min), names_to = "param", values_to = "element") %>%
  filter(param %in% acfElements_min) %>%
  mutate(param = fct_relevel(param, acfElementsList_min)) %>%
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
  ggtitle("ACE Composite - all sites: cps, elements ACF >0.1")
print(Fig3.6)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.6_All_Sites_ACFmin.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig3.6, # First row
          ggarrange(Fig3.5a, Fig3.5b, ncol = 2, labels = c("B", "C")), # Second row with two plots
          nrow = 2, 
          labels = "A", common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.7_All_Sites_ACF_min.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# MAXIMUM ACF threshold element filtering - all sites POB4 included ----------------------------------------

# ACF > 0.5 - ALL SITES DATA
# apply acf based filtering to ACE_xrf_qc elmement list - leaving acf elements >0.1 min threshold
apply(ACE_xrf_qc %>% select(any_of(ACE_elements), Mo_inc, Mo_coh), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>% 
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>% 
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% 
                             pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == lag_thres) %>%
  filter(value >= acf_thres_max) %>% 
  pull(elements) %>% 
  ordered() -> acfElements_max 
acfElements_max

acfElementsList_max <- select(ACE_xrf_qc, any_of(acfElements_max)) %>% 
  names()
acfElementsList_max

# Replot with max acf filtered elements only
Fig3.8 <- apply(ACE_xrf_qc %>% select(any_of(acfElements_max)), 
                2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", size=0.5, linetype = 3) +
  xlim(0, 50) +
  ylim(0.3, 1) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig3.8)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.8_ACFmax_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Summary acf element plot vs depth
library(tidypaleo)
theme_set(theme_bw(base_size=8))
Fig3.9 <- ACE_xrf_qc %>% 
  #filter(qc == TRUE) %>% 
  select(any_of(acfElements_max), depth, label) %>%
  pivot_longer(any_of(acfElements_max), names_to = "param", values_to = "element") %>%
  filter(param %in% acfElements_max) %>%
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
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.9_ACFmax_key_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig3.9, # First row
          ggarrange(Fig3.5a, Fig3.8, ncol = 2, labels = c("B", "C")), # Second row with two plots
          nrow = 2, 
          labels = "A", common.legend = TRUE)
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/qc_acf/Section3/Figures/Fig3.10_All_Sites_ACF_max.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

