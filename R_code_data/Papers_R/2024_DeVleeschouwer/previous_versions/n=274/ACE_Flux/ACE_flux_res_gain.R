# ITRAX-ICPMS resolution gain & flux data 

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
              'cowplot','Hmisc', 'rbacon', 'rioja', 'psych')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

# Set up ------------------------------------------------------------------

#setwd("/Users/Steve/Dropbox/BAS/Data/R/Projects/RBacon") # Macbook 2013
setwd("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon") # Macbook M2
#check working directory
getwd()


# Workflow for subsample data --------------------------------------------------

# STAGES
# 1) Add 95% CI min & max, mean, median ages for min, max and midpoint depths for each record
# 2) Calculate SAR & DMAR - add to dataset - for each min-max depth interval from mean ages and 95% CI errors from min-max depth age ranges
# 3) Calculate resolution of ACE matched subsample dataset & xrf dataset & resolution gain 
# 4) Interpolate GRD (wet density) to 1 mm interval level - to match ITRAX resolution and allow min-max depth matching process 
# 5) Calibrate GRD (wet density) against Subsample dry density data to calculate ITRAX DMAR & flux rates at 1 mm intervals  
# 6) Compare with coh/inc calibration with wet vs dry density calibration?


# Stage 1 & 2 for each site ----------------------------------------------------

# Import matched ICPMS-XRF subsample dataset to work with ----------------------
ACE_matched_all_depths <- read_csv("/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE/cps/ACE_xrf_icp_matched_cps.csv", 
                                   col_names = TRUE, skip = 0)

# BI10 -------------------------------------------------------------------------

# Extract BI10 data only
BI10_matched_all_depths <- ACE_matched_all_depths %>% 
  filter(Site == "BI10")

# BI10 min depth - ages from BACON model
BI10_min_depth <- BI10_matched_all_depths %>% 
  select(min_depth)
BI10_min_depth
write.table(BI10_min_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/BI10/BI10_depths.txt", 
            row.names = FALSE, col.names = FALSE)

Bacon("BI10",depths.file=TRUE,d.max=520,thick=10, res=c(5),cc=3, 
      postbomb=4,rotate.axes=TRUE,rounded=2,
      mem.mean=0.6, mem.strength=20,acc.mean=10,yr.max=10000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=10, plot.pdf = TRUE, 
      title.location='topright', mgp=c(1.5, 0.7, 0))

BI10_min_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/BI10/BI10_54_ages.txt") %>% 
  rename(c(min_depth = depth, SH20_min_depth_min = min, SH20_min_depth_max = max, 
           SH20_min_depth_median = median, SH20_min_depth_mean = mean)) %>% 
  select(SH20_min_depth_min, SH20_min_depth_max, SH20_min_depth_median, SH20_min_depth_mean)
BI10_min_depth_ages

# BI10 max depth - ages from BACON model
BI10_max_depth <- BI10_matched_all_depths %>% 
  select(max_depth)
BI10_max_depth
write.table(BI10_max_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/BI10/BI10_depths.txt", 
            row.names = FALSE, col.names = FALSE)

Bacon("BI10",depths.file=TRUE,d.max=520,thick=10, res=c(5),cc=3, 
      postbomb=4,rotate.axes=TRUE,rounded=2,
      mem.mean=0.6, mem.strength=20,acc.mean=10,yr.max=10000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=10, plot.pdf = TRUE, 
      title.location='topright', mgp=c(1.5, 0.7, 0))

BI10_max_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/BI10/BI10_54_ages.txt") %>% 
  rename(c(max_depth = depth, SH20_max_depth_min = min, SH20_max_depth_max = max, 
           SH20_max_depth_median = median, SH20_max_depth_mean = mean)) %>% 
  select(SH20_max_depth_min, SH20_max_depth_max, SH20_max_depth_median, SH20_max_depth_mean)
BI10_max_depth_ages

# BI10 midpoint depth - ages from BACON model
BI10_midpoint_depth <- BI10_matched_all_depths %>% 
  select(midpoint)
BI10_midpoint_depth
write.table(BI10_midpoint_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/BI10/BI10_depths.txt", 
            row.names = FALSE, col.names = FALSE)

Bacon("BI10",depths.file=TRUE,d.max=520,thick=10, res=c(5),cc=3, 
      postbomb=4,rotate.axes=TRUE,rounded=2,
      mem.mean=0.6, mem.strength=20,acc.mean=10,yr.max=10000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=10, plot.pdf = TRUE, 
      title.location='topright', mgp=c(1.5, 0.7, 0))

BI10_midpoint_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/BI10/BI10_54_ages.txt") %>% 
  rename(c(midpoint_depth = depth, SH20_midpoint_depth_min = min, SH20_midpoint_depth_max = max, 
           SH20_midpoint_depth_median = median, SH20_midpoint_depth_mean = mean)) %>% 
  select(SH20_midpoint_depth_min, SH20_midpoint_depth_max, SH20_midpoint_depth_median, SH20_midpoint_depth_mean)
BI10_midpoint_depth_ages

# Combine into matched dataset file
BI10_all_depths_ages <- bind_cols(BI10_matched_all_depths, 
                             BI10_min_depth_ages, 
                             BI10_max_depth_ages, 
                             BI10_midpoint_depth_ages) %>% 
  relocate(c(SH20_min_depth_min:SH20_midpoint_depth_mean), .after = SH20_age) %>% 
# Calculate sed acc rate across each subsample min-max depth interval 
  mutate(SAR = (max_depth-min_depth) / (SH20_max_depth_mean - SH20_min_depth_mean)) %>%
  mutate(SAR_err = (max_depth-min_depth) / (SH20_max_depth_max - SH20_min_depth_min)) %>% # 95% CI (max depth max age - min depth min age)
  mutate(SAR_pc_err = (SAR_err/SAR)*100) %>% 
  mutate(DMAR = SAR * density_gcm3) %>%
  mutate(DMAR_err = DMAR*((SAR_err/SAR)+(density_gcm3_err/density_gcm3))) %>%
  mutate(DMAR_pc_err = (DMAR_err/DMAR)*100) %>% 
  print()
write.csv(BI10_all_depths_ages,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/BI10_matched_all_depths_ages_cps.csv", row.names = FALSE)

# HER42PB ------------------------------------------

# Extract HER42PB data only
HER42PB_matched_all_depths <- ACE_matched_all_depths %>% 
  filter(Site == "HER42PB")

# HER42PB min depth - ages from BACON model
HER42PB_min_depth <- HER42PB_matched_all_depths %>% 
  select(min_depth)
HER42PB_min_depth
write.table(HER42PB_min_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/HER42PB/HER42PB_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#HER42PB-M2 - min
Bacon("HER42PB",depths.file=TRUE,d.max=500,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=15000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=20, plot.pdf = TRUE, 
      title.location='topright', mgp=c(1.5, 0.7, 0))

HER42PB_min_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/HER42PB/HER42PB_52_ages.txt") %>% 
  rename(c(min_depth = depth, SH20_min_depth_min = min, SH20_min_depth_max = max, 
           SH20_min_depth_median = median, SH20_min_depth_mean = mean)) %>% 
  select(SH20_min_depth_min, SH20_min_depth_max, SH20_min_depth_median, SH20_min_depth_mean)
HER42PB_min_depth_ages

# HER42PB max depth - ages from BACON model
HER42PB_max_depth <- HER42PB_matched_all_depths %>% 
  select(max_depth)
HER42PB_max_depth
write.table(HER42PB_max_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/HER42PB/HER42PB_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#HER42PB-M2 - max
Bacon("HER42PB",depths.file=TRUE,d.max=500,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=15000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=20, plot.pdf = TRUE, 
      title.location='topright', mgp=c(1.5, 0.7, 0))

HER42PB_max_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/HER42PB/HER42PB_52_ages.txt") %>% 
  rename(c(max_depth = depth, SH20_max_depth_min = min, SH20_max_depth_max = max, 
           SH20_max_depth_median = median, SH20_max_depth_mean = mean)) %>% 
  select(SH20_max_depth_min, SH20_max_depth_max, SH20_max_depth_median, SH20_max_depth_mean)
HER42PB_max_depth_ages

# HER42PB midpoint depth - ages from BACON model
HER42PB_midpoint_depth <- HER42PB_matched_all_depths %>% 
  select(midpoint)
HER42PB_midpoint_depth
write.table(HER42PB_midpoint_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/HER42PB/HER42PB_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#HER42PB-M2 - midpoint 
Bacon("HER42PB",depths.file=TRUE,d.max=500,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=15000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=20, plot.pdf = TRUE, 
      title.location='topright', mgp=c(1.5, 0.7, 0))

HER42PB_midpoint_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/HER42PB/HER42PB_52_ages.txt") %>% 
  rename(c(midpoint_depth = depth, SH20_midpoint_depth_min = min, SH20_midpoint_depth_max = max, 
           SH20_midpoint_depth_median = median, SH20_midpoint_depth_mean = mean)) %>% 
  select(SH20_midpoint_depth_min, SH20_midpoint_depth_max, SH20_midpoint_depth_median, SH20_midpoint_depth_mean)
HER42PB_midpoint_depth_ages

# Combine into matched dataset file
HER42PB_all_depths_ages <- bind_cols(HER42PB_matched_all_depths, 
                                  HER42PB_min_depth_ages, 
                                  HER42PB_max_depth_ages, 
                                  HER42PB_midpoint_depth_ages) %>% 
  relocate(c(SH20_min_depth_min:SH20_midpoint_depth_mean), .after = SH20_age) %>% 
  # Calculate sed acc rate across each subsample min-max depth interval 
  mutate(SAR = (max_depth-min_depth) / (SH20_max_depth_mean - SH20_min_depth_mean)) %>%
  mutate(SAR_err = (max_depth-min_depth) / (SH20_max_depth_max - SH20_min_depth_min)) %>% # 95% CI (max depth max age - min depth min age)
  mutate(SAR_pc_err = (SAR_err/SAR)*100) %>% 
  mutate(DMAR = SAR * density_gcm3) %>%
  mutate(DMAR_err = DMAR*((SAR_err/SAR)+(density_gcm3_err/density_gcm3))) %>%
  mutate(DMAR_pc_err = (DMAR_err/DMAR)*100) %>% 
  print()
write.csv(HER42PB_all_depths_ages,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/HER42PB_matched_all_depths_ages_cps.csv", row.names = FALSE)

# KER1 ------------------------------------------

# Extract KER1 data only
KER1_matched_all_depths <- ACE_matched_all_depths %>% 
  filter(Site == "KER1")

# KER1 min depth - ages from BACON model
KER1_min_depth <- KER1_matched_all_depths %>% 
  select(min_depth)
KER1_min_depth
write.table(KER1_min_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER1-M1/KER1-M1_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#KER1-M1 - midpoint
Bacon("KER1-M1",depths.file=TRUE,thick=5,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=20, mem.mean=0.2,rounded=2,
      mem.strength=20,yr.max=12000,d.min=0,d.max=250,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=10, title.location='topright',mgp=c(1.5, 0.7, 0))

KER1_min_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER1-M1/KER1-M1_52_ages.txt") %>% 
  rename(c(min_depth = depth, SH20_min_depth_min = min, SH20_min_depth_max = max, 
           SH20_min_depth_median = median, SH20_min_depth_mean = mean)) %>% 
  select(SH20_min_depth_min, SH20_min_depth_max, SH20_min_depth_median, SH20_min_depth_mean)
KER1_min_depth_ages

# KER1 max depth - ages from BACON model
KER1_max_depth <- KER1_matched_all_depths %>% 
  select(max_depth)
KER1_max_depth
write.table(KER1_max_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER1-M1/KER1-M1_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#KER1-M1 - midpoint
Bacon("KER1-M1",depths.file=TRUE,thick=5,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=20, mem.mean=0.2,rounded=2,
      mem.strength=20,yr.max=12000,d.min=0,d.max=250,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=10, title.location='topright',mgp=c(1.5, 0.7, 0))

KER1_max_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER1-M1/KER1-M1_52_ages.txt") %>% 
  rename(c(max_depth = depth, SH20_max_depth_min = min, SH20_max_depth_max = max, 
           SH20_max_depth_median = median, SH20_max_depth_mean = mean)) %>% 
  select(SH20_max_depth_min, SH20_max_depth_max, SH20_max_depth_median, SH20_max_depth_mean)
KER1_max_depth_ages

# KER1 midpoint depth - ages from BACON model
KER1_midpoint_depth <- KER1_matched_all_depths %>% 
  select(midpoint)
KER1_midpoint_depth
write.table(KER1_midpoint_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER1-M1/KER1-M1_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#KER1-M1 - midpoint
Bacon("KER1-M1",depths.file=TRUE,thick=5,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=20, mem.mean=0.2,rounded=2,
      mem.strength=20,yr.max=12000,d.min=0,d.max=250,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=10, title.location='topright',mgp=c(1.5, 0.7, 0))

KER1_midpoint_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER1-M1/KER1-M1_52_ages.txt") %>% 
  rename(c(midpoint_depth = depth, SH20_midpoint_depth_min = min, SH20_midpoint_depth_max = max, 
           SH20_midpoint_depth_median = median, SH20_midpoint_depth_mean = mean)) %>% 
  select(SH20_midpoint_depth_min, SH20_midpoint_depth_max, SH20_midpoint_depth_median, SH20_midpoint_depth_mean)
KER1_midpoint_depth_ages

# Combine into matched dataset file
KER1_all_depths_ages <- bind_cols(KER1_matched_all_depths, 
                                     KER1_min_depth_ages, 
                                     KER1_max_depth_ages, 
                                     KER1_midpoint_depth_ages) %>% 
  relocate(c(SH20_min_depth_min:SH20_midpoint_depth_mean), .after = SH20_age) %>% 
  # Calculate sed acc rate across each subsample min-max depth interval 
  mutate(SAR = (max_depth-min_depth) / (SH20_max_depth_mean - SH20_min_depth_mean)) %>%
  mutate(SAR_err = (max_depth-min_depth) / (SH20_max_depth_max - SH20_min_depth_min)) %>% # 95% CI (max depth max age - min depth min age)
  mutate(SAR_pc_err = (SAR_err/SAR)*100) %>% 
  mutate(DMAR = SAR * density_gcm3) %>%
  mutate(DMAR_err = DMAR*((SAR_err/SAR)+(density_gcm3_err/density_gcm3))) %>%
  mutate(DMAR_pc_err = (DMAR_err/DMAR)*100) %>% 
  print()
write.csv(KER1_all_depths_ages,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/KER1_matched_all_depths_ages_cps.csv", row.names = FALSE)

# KER3 ------------------------------------------

# Extract KER3 data only
KER3_matched_all_depths <- ACE_matched_all_depths %>% 
  filter(Site == "KER3")

# KER3 min depth - ages from BACON model
KER3_min_depth <- KER3_matched_all_depths %>% 
  select(min_depth)
KER3_min_depth
write.table(KER3_min_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER3-M2/KER3-M2_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#KER3-M2 - min
Bacon("KER3-M2",depths.file=TRUE,thick=10,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=50, mem.mean=0.2,rounded=2,
      mem.strength=20,yr.max=12000,d.min=0,d.max=220,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=50, title.location='topright',mgp=c(1.5, 0.7, 0))

KER3_min_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER3-M2/KER3-M2_24_ages.txt") %>% 
  rename(c(min_depth = depth, SH20_min_depth_min = min, SH20_min_depth_max = max, 
           SH20_min_depth_median = median, SH20_min_depth_mean = mean)) %>% 
  select(SH20_min_depth_min, SH20_min_depth_max, SH20_min_depth_median, SH20_min_depth_mean)
KER3_min_depth_ages

# KER3 max depth - ages from BACON model
KER3_max_depth <- KER3_matched_all_depths %>% 
  select(max_depth)
KER3_max_depth
write.table(KER3_max_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER3-M2/KER3-M2_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#KER3-M2 - max
Bacon("KER3-M2",depths.file=TRUE,thick=10,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=50, mem.mean=0.2,rounded=2,
      mem.strength=20,yr.max=12000,d.min=0,d.max=220,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=50, title.location='topright',mgp=c(1.5, 0.7, 0))

KER3_max_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER3-M2/KER3-M2_24_ages.txt") %>% 
  rename(c(max_depth = depth, SH20_max_depth_min = min, SH20_max_depth_max = max, 
           SH20_max_depth_median = median, SH20_max_depth_mean = mean)) %>% 
  select(SH20_max_depth_min, SH20_max_depth_max, SH20_max_depth_median, SH20_max_depth_mean)
KER3_max_depth_ages

# KER3 midpoint depth - ages from BACON model
KER3_midpoint_depth <- KER3_matched_all_depths %>% 
  select(midpoint)
KER3_midpoint_depth
write.table(KER3_midpoint_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER3-M2/KER3-M2_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#KER3-M2 - midpoint
Bacon("KER3-M2",depths.file=TRUE,thick=10,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=50, mem.mean=0.2,rounded=2,
      mem.strength=20,yr.max=12000,d.min=0,d.max=220,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=50, title.location='topright',mgp=c(1.5, 0.7, 0))

KER3_midpoint_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/KER3-M2/KER3-M2_24_ages.txt") %>% 
  rename(c(midpoint_depth = depth, SH20_midpoint_depth_min = min, SH20_midpoint_depth_max = max, 
           SH20_midpoint_depth_median = median, SH20_midpoint_depth_mean = mean)) %>% 
  select(SH20_midpoint_depth_min, SH20_midpoint_depth_max, SH20_midpoint_depth_median, SH20_midpoint_depth_mean)
KER3_midpoint_depth_ages

# Combine into matched dataset file
KER3_all_depths_ages <- bind_cols(KER3_matched_all_depths, 
                                  KER3_min_depth_ages, 
                                  KER3_max_depth_ages, 
                                  KER3_midpoint_depth_ages) %>% 
  relocate(c(SH20_min_depth_min:SH20_midpoint_depth_mean), .after = SH20_age) %>% 
  # Calculate sed acc rate across each subsample min-max depth interval 
  mutate(SAR = (max_depth-min_depth) / (SH20_max_depth_mean - SH20_min_depth_mean)) %>%
  mutate(SAR_err = (max_depth-min_depth) / (SH20_max_depth_max - SH20_min_depth_min)) %>% # 95% CI (max depth max age - min depth min age)
  mutate(SAR_pc_err = (SAR_err/SAR)*100) %>% 
  mutate(DMAR = SAR * density_gcm3) %>%
  mutate(DMAR_err = DMAR*((SAR_err/SAR)+(density_gcm3_err/density_gcm3))) %>%
  mutate(DMAR_pc_err = (DMAR_err/DMAR)*100) %>% 
  print()
write.csv(KER3_all_depths_ages,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/KER3_matched_all_depths_ages_cps.csv", row.names = FALSE)


# PB1 ------------------------------------------

# Extract PB1 data only
PB1_matched_all_depths <- ACE_matched_all_depths %>% 
  filter(Site == "PB1")

# PB1 min depth - ages from BACON model
PB1_min_depth <- PB1_matched_all_depths %>% 
  select(min_depth)
PB1_min_depth
write.table(PB1_min_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/PB1/PB1_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#PB1-M1 - min
Bacon("PB1",depths.file=TRUE,thick=5,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=50, mem.mean=0.4,rounded=2,
      mem.strength=20,yr.max=10000,d.min=0,d.max=130,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=10, title.location='topright',mgp=c(1.5, 0.7, 0))

PB1_min_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/PB1/PB1_28_ages.txt") %>% 
  rename(c(min_depth = depth, SH20_min_depth_min = min, SH20_min_depth_max = max, 
           SH20_min_depth_median = median, SH20_min_depth_mean = mean)) %>% 
  select(SH20_min_depth_min, SH20_min_depth_max, SH20_min_depth_median, SH20_min_depth_mean)
PB1_min_depth_ages

# PB1 max depth - ages from BACON model
PB1_max_depth <- PB1_matched_all_depths %>% 
  select(max_depth)
PB1_max_depth
write.table(PB1_max_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/PB1/PB1_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#PB1-M1 - max
Bacon("PB1",depths.file=TRUE,thick=5,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=50, mem.mean=0.4,rounded=2,
      mem.strength=20,yr.max=10000,d.min=0,d.max=130,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=10, title.location='topright',mgp=c(1.5, 0.7, 0))

PB1_max_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/PB1/PB1_28_ages.txt") %>% 
  rename(c(max_depth = depth, SH20_max_depth_min = min, SH20_max_depth_max = max, 
           SH20_max_depth_median = median, SH20_max_depth_mean = mean)) %>% 
  select(SH20_max_depth_min, SH20_max_depth_max, SH20_max_depth_median, SH20_max_depth_mean)
PB1_max_depth_ages

# PB1 midpoint depth - ages from BACON model
PB1_midpoint_depth <- PB1_matched_all_depths %>% 
  select(midpoint)
PB1_midpoint_depth
write.table(PB1_midpoint_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/PB1/PB1_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#PB1-M1 - midpoint
Bacon("PB1",depths.file=TRUE,thick=5,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=50, mem.mean=0.4,rounded=2,
      mem.strength=20,yr.max=10000,d.min=0,d.max=130,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=10, title.location='topright',mgp=c(1.5, 0.7, 0))

PB1_midpoint_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/PB1/PB1_28_ages.txt") %>% 
  rename(c(midpoint_depth = depth, SH20_midpoint_depth_min = min, SH20_midpoint_depth_max = max, 
           SH20_midpoint_depth_median = median, SH20_midpoint_depth_mean = mean)) %>% 
  select(SH20_midpoint_depth_min, SH20_midpoint_depth_max, SH20_midpoint_depth_median, SH20_midpoint_depth_mean)
PB1_midpoint_depth_ages

# Combine into matched dataset file
PB1_all_depths_ages <- bind_cols(PB1_matched_all_depths, 
                                     PB1_min_depth_ages, 
                                     PB1_max_depth_ages, 
                                     PB1_midpoint_depth_ages) %>% 
  relocate(c(SH20_min_depth_min:SH20_midpoint_depth_mean), .after = SH20_age) %>% 
  # Calculate sed acc rate across each subsample min-max depth interval 
  mutate(SAR = (max_depth-min_depth) / (SH20_max_depth_mean - SH20_min_depth_mean)) %>%
  mutate(SAR_err = (max_depth-min_depth) / (SH20_max_depth_max - SH20_min_depth_min)) %>% # 95% CI (max depth max age - min depth min age)
  mutate(SAR_pc_err = (SAR_err/SAR)*100) %>% 
  mutate(DMAR = SAR * density_gcm3) %>%
  mutate(DMAR_err = DMAR*((SAR_err/SAR)+(density_gcm3_err/density_gcm3))) %>%
  mutate(DMAR_pc_err = (DMAR_err/DMAR)*100) %>% 
  print()
write.csv(PB1_all_depths_ages,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/PB1_matched_all_depths_ages_cps.csv", row.names = FALSE)

# POB4 ------------------------------------------

# Extract POB4 data only
POB4_matched_all_depths <- ACE_matched_all_depths %>% 
  filter(Site == "POB4")

# POB4 min depth - ages from BACON model
POB4_min_depth <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_min_depths.txt", col_names = TRUE) %>% 
  select(depth)
POB4_min_depth
write.table(POB4_min_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_depths.txt", 
              row.names = FALSE, col.names = FALSE)

#POB4-M1 - min
Bacon("POB4",depths.file=TRUE,thick=10,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=20, mem.mean=0.4,rounded=2,
      mem.strength=20,yr.max=10000,d.min=0,d.max=350,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=20, title.location='topright',mgp=c(1.5, 0.7, 0))

POB4_min_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_37_ages.txt") %>% 
  rename(c(min_depth = depth, SH20_min_depth_min = min, SH20_min_depth_max = max, 
           SH20_min_depth_median = median, SH20_min_depth_mean = mean)) %>% 
  select(min_depth, SH20_min_depth_min, SH20_min_depth_max, SH20_min_depth_median, SH20_min_depth_mean) %>% 
  filter(!min_depth == "69.3") %>% 
  filter(!min_depth == "89.3") %>% 
  filter(!min_depth == "109.3") %>% 
  filter(!min_depth == "129.3") %>% 
  filter(!min_depth == "149.3") %>% 
  filter(!min_depth == "169.3") %>%
  filter(!min_depth == "189.3") %>% 
  select(!min_depth)
POB4_min_depth_ages

# POB4 max depth - ages from BACON model
POB4_max_depth <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_max_depths.txt", col_names = TRUE) %>% 
  select(depth)
POB4_max_depth
write.table(POB4_max_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_depths.txt", 
            row.names = FALSE, col.names = FALSE)

#POB4-M1 - max
Bacon("POB4",depths.file=TRUE,thick=10,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=20, mem.mean=0.4,rounded=2,
      mem.strength=20,yr.max=10000,d.min=0,d.max=350,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=20, title.location='topright',mgp=c(1.5, 0.7, 0))

POB4_max_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_37_ages.txt") %>% 
  rename(c(max_depth = depth, SH20_max_depth_min = min, SH20_max_depth_max = max, 
           SH20_max_depth_median = median, SH20_max_depth_mean = mean)) %>% 
  select(max_depth, SH20_max_depth_min, SH20_max_depth_max, SH20_max_depth_median, SH20_max_depth_mean) %>% 
  filter(!max_depth == "70.7") %>% 
  filter(!max_depth == "90.7") %>% 
  filter(!max_depth == "110.7") %>% 
  filter(!max_depth == "130.7") %>% 
  filter(!max_depth == "150.7") %>% 
  filter(!max_depth == "170.7") %>%
  filter(!max_depth == "190.7") %>% 
  select(!max_depth)
POB4_max_depth_ages

# POB4 midpoint depth - ages from BACON model
POB4_midpoint_depth <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_midpoint_depths.txt", col_names = TRUE) %>% 
  select(depth)
POB4_midpoint_depth
write.table(POB4_midpoint_depth, "/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_depths.txt", 
            row.names = FALSE, col.names = FALSE)


#POB4-M1 - midpoint
Bacon("POB4",depths.file=TRUE,thick=10,rotate.axes=TRUE,cc=3, 
      postbomb=4,acc.mean=20, mem.mean=0.4,rounded=2,
      mem.strength=20,yr.max=10000,d.min=0,d.max=350,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=20, title.location='topright',mgp=c(1.5, 0.7, 0))

POB4_midpoint_depth_ages <- read_table("/Users/sjro/Dropbox/BAS/Data/R/Projects/RBacon/Bacon_runs/POB4/POB4_37_ages.txt") %>% 
  rename(c(midpoint_depth = depth, SH20_midpoint_depth_min = min, SH20_midpoint_depth_max = max, 
           SH20_midpoint_depth_median = median, SH20_midpoint_depth_mean = mean)) %>% 
  select(midpoint_depth, SH20_midpoint_depth_min, SH20_midpoint_depth_max, SH20_midpoint_depth_median, SH20_midpoint_depth_mean) %>% 
  filter(!midpoint_depth == "70") %>% 
  filter(!midpoint_depth == "90") %>% 
  filter(!midpoint_depth == "110") %>% 
  filter(!midpoint_depth == "130") %>% 
  filter(!midpoint_depth == "150") %>% 
  filter(!midpoint_depth == "170") %>%
  filter(!midpoint_depth == "190") %>% 
  select(!midpoint_depth)
POB4_midpoint_depth_ages 

# Combine into matched dataset file
POB4_all_depths_ages <- bind_cols(POB4_matched_all_depths, 
                                 POB4_min_depth_ages, 
                                 POB4_max_depth_ages, 
                                 POB4_midpoint_depth_ages) %>% 
  relocate(c(SH20_min_depth_min:SH20_midpoint_depth_mean), .after = SH20_age) %>% 
  # Calculate sed acc rate across each subsample min-max depth interval 
  mutate(SAR = (max_depth-min_depth) / (SH20_max_depth_mean - SH20_min_depth_mean)) %>%
  mutate(SAR_err = (max_depth-min_depth) / (SH20_max_depth_max - SH20_min_depth_min)) %>% # 95% CI (max depth max age - min depth min age)
  mutate(SAR_pc_err = (SAR_err/SAR)*100) %>% 
  mutate(DMAR = SAR * density_gcm3) %>%
  mutate(DMAR_err = DMAR*((SAR_err/SAR)+(density_gcm3_err/density_gcm3))) %>%
  mutate(DMAR_pc_err = (DMAR_err/DMAR)*100) %>% 
  print()
write.csv(POB4_all_depths_ages,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/POB4_matched_all_depths_ages_cps.csv", row.names = FALSE)




# Combine All Sites into Matched Dataset ---------------------------------------
ACE_xrf_icp_matched_cps_flux <- bind_rows(BI10_all_depths_ages,
                                  HER42PB_all_depths_ages,
                                  KER1_all_depths_ages,
                                  KER3_all_depths_ages,
                                  PB1_all_depths_ages,
                                  POB4_all_depths_ages) %>% 
  print()
write.csv(ACE_xrf_icp_matched_cps_flux,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/ACE_xrf_icp_matched_cps_flux.csv", row.names = FALSE)


# Stage 3: Resolution gain -------------------------------------------------------

# Import ACE flux dataset
ACE_flux <- read_csv("/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/ACE_xrf_icp_matched_cps_flux.csv", 
col_names = TRUE, skip = 0) %>% 
  filter(!(Site == 'POB4')) %>% 
  mutate(subsample_period = SH20_max_depth_mean - SH20_min_depth_min) %>% 
  mutate(subsample_period_err = ((SH20_max_depth_max - SH20_max_depth_mean) + 
                                   (SH20_max_depth_mean - SH20_max_depth_min) + 
                                   (SH20_min_depth_max - SH20_min_depth_mean) + 
                                   (SH20_min_depth_mean - SH20_min_depth_min))/4)  %>%
  mutate(subsample_period_err_pc = (subsample_period_err/subsample_period)*100)  %>%
  print()

ACE_flux$subsample_period_err_pc
  
ACE_flux_stats <- ACE_flux %>%
  select (SAR:subsample_period_err_pc) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_flux,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/ACE_xrf_icp_matched_flux_period.csv", row.names = FALSE)
write.csv(ACE_flux_stats,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Flux/ACE_xrf_icp_matched_flux_period_stats.csv", row.names = FALSE)


# Import ITRAX QC & ACF cps dataframe
ACE_xrf_period0 <- read_csv("/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_Matching/Output/itrax_Composite/Matching_mean/ACE_xrf_qc.csv", 
                      col_names = TRUE, skip = 0)

BI10_xrf_period <- ACE_xrf_period0 %>% 
  filter(Site == "BI10") %>% 
  mutate(spect_period = SH20_age - lag(SH20_age)) %>%
  mutate(spect_period_err = ((SH20_min_age_95CI - lag(SH20_min_age_95CI)) + 
                               (SH20_max_age_95CI - lag(SH20_max_age_95CI)))/2) %>%
  mutate(spect_period_err_pc = (spect_period_err/spect_period)*100)

HER42PB_xrf_period <- ACE_xrf_period0 %>% 
  filter(Site == "HER42PB") %>% 
  mutate(spect_period = SH20_age - lag(SH20_age)) %>%
  mutate(spect_period_err = ((SH20_min_age_95CI - lag(SH20_min_age_95CI)) + 
                               (SH20_max_age_95CI - lag(SH20_max_age_95CI)))/2) %>%
  mutate(spect_period_err_pc = (spect_period_err/spect_period)*100)

KER1_xrf_period <- ACE_xrf_period0 %>% 
  filter(Site == "KER1") %>% 
  mutate(spect_period = SH20_age - lag(SH20_age)) %>%
  mutate(spect_period_err = ((SH20_min_age_95CI - lag(SH20_min_age_95CI)) + 
                               (SH20_max_age_95CI - lag(SH20_max_age_95CI)))/2) %>%
  mutate(spect_period_err_pc = (spect_period_err/spect_period)*100)

KER3_xrf_period <- ACE_xrf_period0 %>% 
  filter(Site == "KER3") %>% 
  mutate(spect_period = SH20_age - lag(SH20_age)) %>%
  mutate(spect_period_err = ((SH20_min_age_95CI - lag(SH20_min_age_95CI)) + 
                               (SH20_max_age_95CI - lag(SH20_max_age_95CI)))/2) %>%
  mutate(spect_period_err_pc = (spect_period_err/spect_period)*100)

PB1_xrf_period <- ACE_xrf_period0 %>% 
  filter(Site == "PB1") %>% 
  mutate(spect_period = SH20_age - lag(SH20_age)) %>%
  mutate(spect_period_err = ((SH20_min_age_95CI - lag(SH20_min_age_95CI)) + 
                               (SH20_max_age_95CI - lag(SH20_max_age_95CI)))/2) %>%
  mutate(spect_period_err_pc = (spect_period_err/spect_period)*100)

# Combine all sites
ACE_xrf_period <- bind_rows(BI10_xrf_period,HER42PB_xrf_period,KER1_xrf_period,
                          KER3_xrf_period, PB1_xrf_period) %>% 
  print()

# Summary stats
ACE_xrf_period_stats <- ACE_xrf_period %>%
  select (spect_period:spect_period_err_pc) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_xrf_period,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_flux/ACE_xrf_period.csv", row.names = FALSE)
write.csv(ACE_xrf_period,"/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/ACE_flux/ACE_xrf_period_stats.csv", row.names = FALSE)



# Stage 4 - Wet MSCL GRD vs Dry subsample density calibration for each site ---------------- TO FINISH FOR DUST PAPER


# Import ACE GEOTEK_MSCL dataset to work with ----------------------
ACE_MSCL <- read_csv("/Users/sjro/Dropbox/BAS/Data/R/Papers_R/2024_DeVleeschouwer/Data/ACE_SHW_MSCL_Composite.csv", 
                                   col_names = TRUE, skip = 0)

# BI10 0.1 cm res wet density data ---------------------------------------------

# Extract BI10 data only
BI10_MSCL <- ACE_MSCL %>% 
  filter(Site == "BI10") %>% 
  select(depth, Den1_SAT)

# create dataframe for bas R
BI10_MSCL_df <- data.frame(BI10_MSCL)

# define variables to interpolate at 0.1 cm resolution
den <- BI10_MSCL_df$Den1_SAT
depth_data <- BI10_MSCL_df$depth
start_depth <- 1
end_depth <- 500

# Target depth resolution (0.1 cm) for depth of record 
target_depth <- seq(start_depth, end_depth, by = 0.1)  # Depth data at 0.1 cm resolution

# Perform linear interpolation
interpolated_den <- approx(depth_data, values, xout = target_depth)$y

# Output the interpolated values
BI10_Den_0.1cm <- data.frame(Depth = target_depth, den_interpolated = interpolated_values)






