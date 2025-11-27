# ACE subsample average interval - depth and time

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

# Load libraries ---------------------------------------------------------------
packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'purrr', 'readr', 'ggpubr', 'patchwork',
              'gridExtra', 'cowplot', 'vegan', 'rioja', 'ellipse', 'factoextra',
              'reshape2', 'GGally', 'ggsci', 'ggdendro', 'dendextend', 'dynamicTreeCut',
              'colorspace', 'cluster', 'magrittr', 'mgcv', 'gtable', 'repr',
              'bestNormalize','sjmisc', 'chemometrics', 'compositions', 
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis' )
lapply(packages, library, character.only=TRUE)


# Depth Intervals --------------------------------------------------------------

# ACE - Depth & interval stats -------------------------------------------------
ACE_xrf_icpms_matched_depth <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, sample, min_depth, max_depth, depth, depth_err, SH20_mean_age, SH20_mean_95CI) %>% 
  mutate(thickness_cm = max_depth - min_depth) %>% 
# calculate depth difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    depth_next = lead(depth), # puts next row value in new column
    depth_next_err = lead(depth_err), # puts next row value in new column
    interval_cm = depth_next - depth,# calculates next row - row
    interval_cm_err = sqrt(depth_err^2 + depth_next_err^2) #propagates error between rows in calculation
  ) %>%
  ungroup()
ACE_xrf_icpms_matched_depth
write.csv(ACE_xrf_icpms_matched_depth,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/ACE_interval_depth.csv", 
          row.names = FALSE)

# Average and stdev
ACE_xrf_icpms_depth_stats <- ACE_xrf_icpms_matched_depth %>%
  summarise(
    Site = "ACE",
    n_rows = n(),
    mean_thickness = mean(thickness_cm, na.rm = TRUE),
    sd_thickness = sd(thickness_cm, na.rm = TRUE),
    se_thickness = sd_thickness/sqrt(n_rows),
    mean_interval_cm = mean(interval_cm, na.rm = TRUE),
    mean_interval_cm_err = mean(interval_cm_err, na.rm = TRUE),
    sd_interval_cm = sd(interval_cm, na.rm = TRUE),
    se_interval_cm = sd(interval_cm, na.rm = TRUE)/sqrt(n_rows),
    min_interval_cm = min(interval_cm, na.rm = TRUE),
    max_interval_cm = max(interval_cm, na.rm = TRUE)
  ) 
ACE_xrf_icpms_depth_stats
write.csv(ACE_xrf_icpms_depth_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/ACE_interval_stats_depth.csv", 
          row.names = FALSE)

# BI10 - Depth & interval stats ------------------------------------------------
BI10_xrf_icpms_matched_depth <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "BI10") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, 
         SH20_mean_age, SH20_mean_95CI) %>% 
  mutate(thickness_cm = max_depth - min_depth) %>% 
  # calculate depth difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    depth_next = lead(depth), # puts next row value in new column
    depth_next_err = lead(depth_err), # puts next row value in new column
    interval_cm = depth_next - depth,# calculates next row - row
    interval_cm_err = sqrt(depth_err^2 + depth_next_err^2) #propagates error between rows in calculation
  ) %>%
  ungroup()
BI10_xrf_icpms_matched_depth
write.csv(BI10_xrf_icpms_matched_depth,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/BI10_interval_depth.csv", 
          row.names = FALSE)

# Average and stdev
BI10_xrf_icpms_depth_stats <- BI10_xrf_icpms_matched_depth %>%
  summarise(
    Site = "BI10",
    n_rows = n(),
    mean_thickness = mean(thickness_cm, na.rm = TRUE),
    sd_thickness = sd(thickness_cm, na.rm = TRUE),
    se_thickness = sd_thickness/sqrt(n_rows),
    mean_interval_cm = mean(interval_cm, na.rm = TRUE),
    mean_interval_cm_err = mean(interval_cm_err, na.rm = TRUE),
    sd_interval_cm = sd(interval_cm, na.rm = TRUE),
    se_interval_cm = sd(interval_cm, na.rm = TRUE)/sqrt(n_rows),
    min_interval_cm = min(interval_cm, na.rm = TRUE),
    max_interval_cm = max(interval_cm, na.rm = TRUE)
  ) 
BI10_xrf_icpms_depth_stats
write.csv(BI10_xrf_icpms_depth_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/BI10_interval_stats_depth.csv", 
          row.names = FALSE)

# HER42PB - Depth & interval stats ------------------------------------------------
HER42PB_xrf_icpms_matched_depth <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "HER42PB") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, 
         SH20_mean_age, SH20_mean_95CI) %>% 
  mutate(thickness_cm = max_depth - min_depth) %>% 
  # calculate depth difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    depth_next = lead(depth), # puts next row value in new column
    depth_next_err = lead(depth_err), # puts next row value in new column
    interval_cm = depth_next - depth,# calculates next row - row
    interval_cm_err = sqrt(depth_err^2 + depth_next_err^2) #propagates error between rows in calculation
  ) %>%
  ungroup()
HER42PB_xrf_icpms_matched_depth
write.csv(HER42PB_xrf_icpms_matched_depth,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/HER42PB_interval_depth.csv", 
          row.names = FALSE)

# Average and stdev
HER42PB_xrf_icpms_depth_stats <- HER42PB_xrf_icpms_matched_depth %>%
  summarise(
    Site = "HER42PB",
    n_rows = n(),
    mean_thickness = mean(thickness_cm, na.rm = TRUE),
    sd_thickness = sd(thickness_cm, na.rm = TRUE),
    se_thickness = sd_thickness/sqrt(n_rows),
    mean_interval_cm = mean(interval_cm, na.rm = TRUE),
    mean_interval_cm_err = mean(interval_cm_err, na.rm = TRUE),
    sd_interval_cm = sd(interval_cm, na.rm = TRUE),
    se_interval_cm = sd(interval_cm, na.rm = TRUE)/sqrt(n_rows),
    min_interval_cm = min(interval_cm, na.rm = TRUE),
    max_interval_cm = max(interval_cm, na.rm = TRUE)
  ) 
HER42PB_xrf_icpms_depth_stats
write.csv(HER42PB_xrf_icpms_depth_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/HER42PB_interval_stats_depth.csv", 
          row.names = FALSE)
# KER1 - Depth & interval stats ------------------------------------------------
KER1_xrf_icpms_matched_depth <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "KER1") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, 
         SH20_mean_age, SH20_mean_95CI) %>% 
  mutate(thickness_cm = max_depth - min_depth) %>% 
  # calculate depth difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    depth_next = lead(depth), # puts next row value in new column
    depth_next_err = lead(depth_err), # puts next row value in new column
    interval_cm = depth_next - depth,# calculates next row - row
    interval_cm_err = sqrt(depth_err^2 + depth_next_err^2) #propagates error between rows in calculation
  ) %>%
  ungroup()
KER1_xrf_icpms_matched_depth
write.csv(KER1_xrf_icpms_matched_depth,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER1_interval_depth.csv", 
          row.names = FALSE)

# Average and stdev
KER1_xrf_icpms_depth_stats <- KER1_xrf_icpms_matched_depth %>%
  summarise(
    Site = "KER1",
    n_rows = n(),
    mean_thickness = mean(thickness_cm, na.rm = TRUE),
    sd_thickness = sd(thickness_cm, na.rm = TRUE),
    se_thickness = sd_thickness/sqrt(n_rows),
    mean_interval_cm = mean(interval_cm, na.rm = TRUE),
    mean_interval_cm_err = mean(interval_cm_err, na.rm = TRUE),
    sd_interval_cm = sd(interval_cm, na.rm = TRUE),
    se_interval_cm = sd(interval_cm, na.rm = TRUE)/sqrt(n_rows),
    min_interval_cm = min(interval_cm, na.rm = TRUE),
    max_interval_cm = max(interval_cm, na.rm = TRUE)
  ) 
KER1_xrf_icpms_depth_stats
write.csv(KER1_xrf_icpms_depth_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER1_interval_stats_depth.csv", 
          row.names = FALSE)

# KER3 - Depth & interval stats ------------------------------------------------
KER3_xrf_icpms_matched_depth <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "KER3") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, 
         SH20_mean_age, SH20_mean_95CI) %>% 
  mutate(thickness_cm = max_depth - min_depth) %>% 
  # calculate depth difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    depth_next = lead(depth), # puts next row value in new column
    depth_next_err = lead(depth_err), # puts next row value in new column
    interval_cm = depth_next - depth,# calculates next row - row
    interval_cm_err = sqrt(depth_err^2 + depth_next_err^2) #propagates error between rows in calculation
  ) %>%
  ungroup()
KER3_xrf_icpms_matched_depth
write.csv(KER3_xrf_icpms_matched_depth,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER3_interval_depth.csv", 
          row.names = FALSE)

# Average and stdev
KER3_xrf_icpms_depth_stats <- KER3_xrf_icpms_matched_depth %>%
  summarise(
    Site = "KER3",
    n_rows = n(),
    mean_thickness = mean(thickness_cm, na.rm = TRUE),
    sd_thickness = sd(thickness_cm, na.rm = TRUE),
    se_thickness = sd_thickness/sqrt(n_rows),
    mean_interval_cm = mean(interval_cm, na.rm = TRUE),
    mean_interval_cm_err = mean(interval_cm_err, na.rm = TRUE),
    sd_interval_cm = sd(interval_cm, na.rm = TRUE),
    se_interval_cm = sd(interval_cm, na.rm = TRUE)/sqrt(n_rows),
    min_interval_cm = min(interval_cm, na.rm = TRUE),
    max_interval_cm = max(interval_cm, na.rm = TRUE)
  ) 
KER3_xrf_icpms_depth_stats
write.csv(KER3_xrf_icpms_depth_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER3_interval_stats_depth.csv", 
          row.names = FALSE)

# PB1 - Depth & interval stats ------------------------------------------------
PB1_xrf_icpms_matched_depth <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "PB1") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, 
         SH20_mean_age, SH20_mean_95CI) %>% 
  mutate(thickness_cm = max_depth - min_depth) %>% 
  # calculate depth difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    depth_next = lead(depth), # puts next row value in new column
    depth_next_err = lead(depth_err), # puts next row value in new column
    interval_cm = depth_next - depth,# calculates next row - row
    interval_cm_err = sqrt(depth_err^2 + depth_next_err^2) #propagates error between rows in calculation
  ) %>%
  ungroup()
PB1_xrf_icpms_matched_depth
write.csv(PB1_xrf_icpms_matched_depth,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/PB1_interval_depth.csv", 
          row.names = FALSE)

# Average and stdev
PB1_xrf_icpms_depth_stats <- PB1_xrf_icpms_matched_depth %>%
  summarise(
    Site = "PB1",
    n_rows = n(),
    mean_thickness = mean(thickness_cm, na.rm = TRUE),
    sd_thickness = sd(thickness_cm, na.rm = TRUE),
    se_thickness = sd_thickness/sqrt(n_rows),
    mean_interval_cm = mean(interval_cm, na.rm = TRUE),
    mean_interval_cm_err = mean(interval_cm_err, na.rm = TRUE),
    sd_interval_cm = sd(interval_cm, na.rm = TRUE),
    se_interval_cm = sd(interval_cm, na.rm = TRUE)/sqrt(n_rows),
    min_interval_cm = min(interval_cm, na.rm = TRUE),
    max_interval_cm = max(interval_cm, na.rm = TRUE)
  ) 
PB1_xrf_icpms_depth_stats
write.csv(PB1_xrf_icpms_depth_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/PB1_interval_stats_depth.csv", 
          row.names = FALSE)


# Create ACE Depth & intervals stats file --------------------------------------

ACE_depth_interval_stats_all <- bind_rows(ACE_xrf_icpms_depth_stats, 
                                      BI10_xrf_icpms_depth_stats,
                                        HER42PB_xrf_icpms_depth_stats,
                                        KER1_xrf_icpms_depth_stats,
                                        KER3_xrf_icpms_depth_stats,
                                        PB1_xrf_icpms_depth_stats) %>%
  rename(n = n_rows) %>% 
  select(-se_thickness) %>%  # remove column - values are less than 1 dp
  mutate(across(
    .cols = -c(Site, n),                       # exclude Site and n
    .fns = ~ sprintf("%.1f", round(.x, 1))    # round and format with 1 decimal so that eg 6 becomes 6.0
  ))

ACE_depth_interval_stats_all
tail(ACE_depth_interval_stats_all)
write.csv(ACE_depth_interval_stats_all,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/ACE_depth_interval_stats_all.csv", 
          row.names = FALSE)



# Age Intervals ----------------------------------------------------------------

# ACE - Age & interval stats  --------------------------------------------------
ACE_xrf_icpms_matched_age <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  select(Site, sample, min_depth, max_depth, depth, depth_err, SH20_min_age_95CI, SH20_max_age_95CI, SH20_mean_age, SH20_mean_95CI, accrate, accrate_err) %>% 
  mutate(max_95CI_range_yrs = SH20_max_age_95CI - SH20_min_age_95CI) %>% 
  # calculate age difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    age_next = lead(SH20_mean_age), # puts next row value in new column
    age_next_err = lead(SH20_mean_95CI), # puts next row value in new column
    interval_yrs = age_next - SH20_mean_age,# calculates next row - row
    interval_yrs_err = abs(age_next_err - SH20_mean_95CI),
  ) %>%
  ungroup()
ACE_xrf_icpms_matched_age
write.csv(ACE_xrf_icpms_matched_age,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/ACE_interval_age.csv", 
          row.names = FALSE)


# Average and stdev
ACE_xrf_icpms_age_stats <- ACE_xrf_icpms_matched_age %>%
  summarise(
    Site = "ACE",
    n_rows = n(),
    mean_accrate = mean(accrate, na.rm = TRUE),
    sd_accrate = sd(accrate, na.rm = TRUE),
    se_accrate = sd_accrate/sqrt(n_rows),
    mean_max_95CI_range = mean(max_95CI_range_yrs, na.rm = TRUE),
    sd_max_95CI_range = sd(max_95CI_range_yrs, na.rm = TRUE),
    se_max_95CI_range = sd_max_95CI_range/sqrt(n_rows),
    mean_interval_yrs = mean(interval_yrs, na.rm = TRUE),
    mean_interval_yrs_err = mean(interval_yrs_err, na.rm = TRUE),
    sd_interval_yrs = sd(interval_yrs, na.rm = TRUE),
    se_interval_yrs = sd(interval_yrs, na.rm = TRUE)/sqrt(n_rows),
    min_interval_yrs = min(interval_yrs, na.rm = TRUE),
    max_interval_yrs = max(interval_yrs, na.rm = TRUE)
  ) 
ACE_xrf_icpms_age_stats
write.csv(ACE_xrf_icpms_age_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/ACE_interval_stats_age.csv", 
          row.names = FALSE)

# BI10 - Age & interval stats  -------------------------------------------------
BI10_xrf_icpms_matched_age <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "BI10") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, SH20_min_age_95CI, 
         SH20_max_age_95CI, SH20_mean_age, SH20_mean_95CI, accrate, accrate_err) %>% 
  mutate(max_95CI_range_yrs = SH20_max_age_95CI - SH20_min_age_95CI) %>% 
  # calculate age difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    age_next = lead(SH20_mean_age), # puts next row value in new column
    age_next_err = lead(SH20_mean_95CI), # puts next row value in new column
    interval_yrs = age_next - SH20_mean_age,# calculates next row - row
    interval_yrs_err = abs(age_next_err - SH20_mean_95CI),
  ) %>%
  ungroup()
BI10_xrf_icpms_matched_age
write.csv(BI10_xrf_icpms_matched_age,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/BI10_interval_age.csv", 
          row.names = FALSE)

# Average and stdev
BI10_xrf_icpms_age_stats <- BI10_xrf_icpms_matched_age %>%
  summarise(
    Site = "BI10",
    n_rows = n(),
    mean_accrate = mean(accrate, na.rm = TRUE),
    sd_accrate = sd(accrate, na.rm = TRUE),
    se_accrate = sd_accrate/sqrt(n_rows),
    mean_max_95CI_range = mean(max_95CI_range_yrs, na.rm = TRUE),
    sd_max_95CI_range = sd(max_95CI_range_yrs, na.rm = TRUE),
    se_max_95CI_range = sd_max_95CI_range/sqrt(n_rows),
    mean_interval_yrs = mean(interval_yrs, na.rm = TRUE),
    mean_interval_yrs_err = mean(interval_yrs_err, na.rm = TRUE),
    sd_interval_yrs = sd(interval_yrs, na.rm = TRUE),
    se_interval_yrs = sd(interval_yrs, na.rm = TRUE)/sqrt(n_rows),
    min_interval_yrs = min(interval_yrs, na.rm = TRUE),
    max_interval_yrs = max(interval_yrs, na.rm = TRUE),
  ) 
BI10_xrf_icpms_age_stats
write.csv(BI10_xrf_icpms_age_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/BI10_interval_stats_age.csv", 
          row.names = FALSE)

# HER42PB - Age & interval stats  -------------------------------------------------
HER42PB_xrf_icpms_matched_age <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "HER42PB") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, SH20_min_age_95CI, 
         SH20_max_age_95CI, SH20_mean_age, SH20_mean_95CI, accrate, accrate_err) %>% 
  mutate(max_95CI_range_yrs = SH20_max_age_95CI - SH20_min_age_95CI) %>% 
  # calculate age difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    age_next = lead(SH20_mean_age), # puts next row value in new column
    age_next_err = lead(SH20_mean_95CI), # puts next row value in new column
    interval_yrs = age_next - SH20_mean_age,# calculates next row - row
    interval_yrs_err = abs(age_next_err - SH20_mean_95CI),
  ) %>%
  ungroup()
HER42PB_xrf_icpms_matched_age
write.csv(HER42PB_xrf_icpms_matched_age,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/HER42PB_interval_age.csv", 
          row.names = FALSE)

# Average and stdev
HER42PB_xrf_icpms_age_stats <- HER42PB_xrf_icpms_matched_age %>%
  summarise(
    Site = "HER42PB",
    n_rows = n(),
    mean_accrate = mean(accrate, na.rm = TRUE),
    sd_accrate = sd(accrate, na.rm = TRUE),
    se_accrate = sd_accrate/sqrt(n_rows),
    mean_max_95CI_range = mean(max_95CI_range_yrs, na.rm = TRUE),
    sd_max_95CI_range = sd(max_95CI_range_yrs, na.rm = TRUE),
    se_max_95CI_range = sd_max_95CI_range/sqrt(n_rows),
    mean_interval_yrs = mean(interval_yrs, na.rm = TRUE),
    mean_interval_yrs_err = mean(interval_yrs_err, na.rm = TRUE),
    sd_interval_yrs = sd(interval_yrs, na.rm = TRUE),
    se_interval_yrs = sd(interval_yrs, na.rm = TRUE)/sqrt(n_rows),
    min_interval_yrs = min(interval_yrs, na.rm = TRUE),
    max_interval_yrs = max(interval_yrs, na.rm = TRUE),
  ) 
HER42PB_xrf_icpms_age_stats
write.csv(HER42PB_xrf_icpms_age_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/HER42PB_interval_stats_age.csv", 
          row.names = FALSE)


# KER1 - Age & interval stats  -------------------------------------------------
KER1_xrf_icpms_matched_age <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "KER1") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, SH20_min_age_95CI, 
         SH20_max_age_95CI, SH20_mean_age, SH20_mean_95CI, accrate, accrate_err) %>% 
  mutate(max_95CI_range_yrs = SH20_max_age_95CI - SH20_min_age_95CI) %>% 
  # calculate age difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    age_next = lead(SH20_mean_age), # puts next row value in new column
    age_next_err = lead(SH20_mean_95CI), # puts next row value in new column
    interval_yrs = age_next - SH20_mean_age,# calculates next row - row
    interval_yrs_err = abs(age_next_err - SH20_mean_95CI),
  ) %>%
  ungroup()
KER1_xrf_icpms_matched_age
write.csv(KER1_xrf_icpms_matched_age,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER1_interval_age.csv", 
          row.names = FALSE)

# Average and stdev
KER1_xrf_icpms_age_stats <- KER1_xrf_icpms_matched_age %>%
  summarise(
    Site = "KER1",
    n_rows = n(),
    mean_accrate = mean(accrate, na.rm = TRUE),
    sd_accrate = sd(accrate, na.rm = TRUE),
    se_accrate = sd_accrate/sqrt(n_rows),
    mean_max_95CI_range = mean(max_95CI_range_yrs, na.rm = TRUE),
    sd_max_95CI_range = sd(max_95CI_range_yrs, na.rm = TRUE),
    se_max_95CI_range = sd_max_95CI_range/sqrt(n_rows),
    mean_interval_yrs = mean(interval_yrs, na.rm = TRUE),
    mean_interval_yrs_err = mean(interval_yrs_err, na.rm = TRUE),
    sd_interval_yrs = sd(interval_yrs, na.rm = TRUE),
    se_interval_yrs = sd(interval_yrs, na.rm = TRUE)/sqrt(n_rows),
    min_interval_yrs = min(interval_yrs, na.rm = TRUE),
    max_interval_yrs = max(interval_yrs, na.rm = TRUE),
  ) 
KER1_xrf_icpms_age_stats
write.csv(KER1_xrf_icpms_age_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER1_interval_stats_age.csv", 
          row.names = FALSE)


# KER3 - Age & interval stats  -------------------------------------------------
KER3_xrf_icpms_matched_age <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "KER3") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, SH20_min_age_95CI, 
         SH20_max_age_95CI, SH20_mean_age, SH20_mean_95CI, accrate, accrate_err) %>% 
  mutate(max_95CI_range_yrs = SH20_max_age_95CI - SH20_min_age_95CI) %>% 
  # calculate age difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    age_next = lead(SH20_mean_age), # puts next row value in new column
    age_next_err = lead(SH20_mean_95CI), # puts next row value in new column
    interval_yrs = age_next - SH20_mean_age,# calculates next row - row
    interval_yrs_err = abs(age_next_err - SH20_mean_95CI),
  ) %>%
  ungroup()
KER3_xrf_icpms_matched_age
write.csv(KER3_xrf_icpms_matched_age,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER3_interval_age.csv", 
          row.names = FALSE)

# Average and stdev
KER3_xrf_icpms_age_stats <- KER3_xrf_icpms_matched_age %>%
  summarise(
    Site = "KER3",
    n_rows = n(),
    mean_accrate = mean(accrate, na.rm = TRUE),
    sd_accrate = sd(accrate, na.rm = TRUE),
    se_accrate = sd_accrate/sqrt(n_rows),
    mean_max_95CI_range = mean(max_95CI_range_yrs, na.rm = TRUE),
    sd_max_95CI_range = sd(max_95CI_range_yrs, na.rm = TRUE),
    se_max_95CI_range = sd_max_95CI_range/sqrt(n_rows),
    mean_interval_yrs = mean(interval_yrs, na.rm = TRUE),
    mean_interval_yrs_err = mean(interval_yrs_err, na.rm = TRUE),
    sd_interval_yrs = sd(interval_yrs, na.rm = TRUE),
    se_interval_yrs = sd(interval_yrs, na.rm = TRUE)/sqrt(n_rows),
    min_interval_yrs = min(interval_yrs, na.rm = TRUE),
    max_interval_yrs = max(interval_yrs, na.rm = TRUE),
  ) 
KER3_xrf_icpms_age_stats
write.csv(KER3_xrf_icpms_age_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/KER3_interval_stats_age.csv", 
          row.names = FALSE)


# PB1 - Age & interval stats  -------------------------------------------------
PB1_xrf_icpms_matched_age <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_subsample_icp_xrf_matched_cps.csv") %>%
  filter(Site == "PB1") %>% 
  select(Site, sample, min_depth, max_depth, depth, depth_err, SH20_min_age_95CI, 
         SH20_max_age_95CI, SH20_mean_age, SH20_mean_95CI, accrate, accrate_err) %>% 
  mutate(max_95CI_range_yrs = SH20_max_age_95CI - SH20_min_age_95CI) %>% 
  # calculate age difference to next subsample in row and its propagated error
  arrange(Site, sample) %>% # ensure proper ordering
  group_by(Site) %>%
  mutate(
    age_next = lead(SH20_mean_age), # puts next row value in new column
    age_next_err = lead(SH20_mean_95CI), # puts next row value in new column
    interval_yrs = age_next - SH20_mean_age,# calculates next row - row
    interval_yrs_err = abs(age_next_err - SH20_mean_95CI),
  ) %>%
  ungroup()
PB1_xrf_icpms_matched_age
write.csv(PB1_xrf_icpms_matched_age,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/PB1_interval_age.csv", 
          row.names = FALSE)

# Average and stdev
PB1_xrf_icpms_age_stats <- PB1_xrf_icpms_matched_age %>%
  summarise(
    Site = "PB1",
    n_rows = n(),
    mean_accrate = mean(accrate, na.rm = TRUE),
    sd_accrate = sd(accrate, na.rm = TRUE),
    se_accrate = sd_accrate/sqrt(n_rows),
    mean_max_95CI_range = mean(max_95CI_range_yrs, na.rm = TRUE),
    sd_max_95CI_range = sd(max_95CI_range_yrs, na.rm = TRUE),
    se_max_95CI_range = sd_max_95CI_range/sqrt(n_rows),
    mean_interval_yrs = mean(interval_yrs, na.rm = TRUE),
    mean_interval_yrs_err = mean(interval_yrs_err, na.rm = TRUE),
    sd_interval_yrs = sd(interval_yrs, na.rm = TRUE),
    se_interval_yrs = sd(interval_yrs, na.rm = TRUE)/sqrt(n_rows),
    min_interval_yrs = min(interval_yrs, na.rm = TRUE),
    max_interval_yrs = max(interval_yrs, na.rm = TRUE),
  ) 
PB1_xrf_icpms_age_stats
write.csv(PB1_xrf_icpms_age_stats,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/PB1_interval_stats_age.csv", 
          row.names = FALSE)


# Create ACE Age & intervals stats file --------------------------------------

ACE_age_interval_stats_all <- bind_rows(ACE_xrf_icpms_age_stats, 
                                          BI10_xrf_icpms_age_stats,
                                          HER42PB_xrf_icpms_age_stats,
                                          KER1_xrf_icpms_age_stats,
                                          KER3_xrf_icpms_age_stats,
                                          PB1_xrf_icpms_age_stats) %>%
  rename(n = n_rows) %>% 
  select(-c(se_accrate, mean_max_95CI_range, sd_max_95CI_range, se_max_95CI_range, min_interval_yrs, max_interval_yrs)) %>% 
  mutate(across(
    .cols = -c(Site, n, mean_accrate, sd_accrate),   # exclude these columns
    .fns = ~ round(.x, 0)                            # round to nearest integer
  ))

ACE_age_interval_stats_all
tail(ACE_age_interval_stats_all)
write.csv(ACE_age_interval_stats_all,"Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Interval_stats/ACE_age_interval_stats_all.csv", 
          row.names = FALSE)



# END ---------------------------------------------------------------------
