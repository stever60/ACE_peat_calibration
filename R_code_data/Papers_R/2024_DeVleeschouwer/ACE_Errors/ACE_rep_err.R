# ITRAX-ICPMS Reproducibility analysis - cps

# Clear previous console
remove (list = ls())
# Set working directory - Macbook Pro M2
setwd("/Users/sjro/Dropbox/BAS/Data/R/")
getwd()
# clear plot window
dev.off()

# Load libraries and define elements list ----------------------------------

packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'readr', 'ggpubr', 'patchwork',
              'gridExtra', 'cowplot', 'vegan', 'rioja', 'ellipse', 'factoextra',
              'reshape2', 'GGally', 'ggsci', 'ggdendro', 'dendextend', 'dynamicTreeCut',
              'colorspace', 'cluster', 'magrittr', 'mgcv', 'gtable', 'repr',
              'bestNormalize','sjmisc', 'chemometrics', 'compositions', 
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis',
              'ggrepel', 'itraxR','PeriodicTable','errors','patchwork',
              'forecast','directlabels','broom','performance','lmtest','ggpmisc',
              'cowplot','Hmisc', 'errors')
lapply(packages, library, character.only=TRUE)
options(scipen = 999)

# key elements
elementsList <- c("Ca", "Ti", "Mn", "Fe","Sr", "Zr")

# Marion Island - AW/DH core
# MI1A --------------------------------------------------------------------
core_title = "MI1A"
# Import data and remove validity = FALSE rows
MI1A_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1A_rep1/Results_ACE.txt")
  max_position = max(MI1A_REP2$position, na.rm = TRUE)
MI1A_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1A_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1A_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1A_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1A_rep_err <- list(MI1A_REP1, MI1A_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1A_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1A_rep_err.txt")
print(MI1A_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1A <- MI1A_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1A
# MI1B --------------------------------------------------------------------
core_title = "MI1B"
# Import data and remove validity = FALSE rows
MI1B_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1B_rep1/Results_ACE.txt")
max_position = max(MI1B_REP2$position, na.rm = TRUE)
MI1B_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1B_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1B_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1B_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1B_rep_err <- list(MI1B_REP1, MI1B_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1B_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1B_rep_err.txt")
print(MI1B_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1B <- MI1B_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1B
# MI1C --------------------------------------------------------------------
core_title = "MI1C"
# Import data and remove validity = FALSE rows
MI1C_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1C_rep1/Results_ACE.txt")
max_position = max(MI1C_REP2$position, na.rm = TRUE)
MI1C_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1C_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1C_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1C_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1C_rep_err <- list(MI1C_REP1, MI1C_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1C_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1C_rep_err.txt")
print(MI1C_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1C <- MI1C_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1C
# MI1D --------------------------------------------------------------------
core_title = "MI1D"
# Import data and remove validity = FALSE rows
MI1D_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1D_rep1/Results_ACE.txt")
max_position = max(MI1D_REP2$position, na.rm = TRUE)
MI1D_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1D_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1D_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1D_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1D_rep_err <- list(MI1D_REP1, MI1D_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1D_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1D_rep_err.txt")
print(MI1D_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1D <- MI1D_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1D
# MI1E --------------------------------------------------------------------
core_title = "MI1E"
# Import data and remove validity = FALSE rows
MI1E_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1E_rep1/Results_ACE.txt")
max_position = max(MI1E_REP2$position, na.rm = TRUE)
MI1E_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1E_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1E_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1E_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1E_rep_err <- list(MI1E_REP1, MI1E_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1E_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1E_rep_err.txt")
print(MI1E_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1E <- MI1E_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1E
# MI1F --------------------------------------------------------------------
core_title = "MI1F"
# Import data and remove validity = FALSE rows
MI1F_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1F_rep1/Results_ACE.txt")
max_position = max(MI1F_REP2$position, na.rm = TRUE)
MI1F_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1F_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1F_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1F_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1F_rep_err <- list(MI1F_REP1, MI1F_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1F_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1F_rep_err.txt")
print(MI1F_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1F <- MI1F_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1F
# MI1G --------------------------------------------------------------------
core_title = "MI1G"
# Import data and remove validity = FALSE rows
MI1G_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1G_rep1/Results_ACE.txt")
max_position = max(MI1G_REP2$position, na.rm = TRUE)
MI1G_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1G_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1G_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1G_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1G_rep_err <- list(MI1G_REP1, MI1G_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1G_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1G_rep_err.txt")
print(MI1G_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1G <- MI1G_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1G
 # xrf and rep1 missing from currently missing from dataset
# MI1H --------------------------------------------------------------------
core_title = "MI1H"
# Import data and remove validity = FALSE rows
MI1H_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1H_rep1/Results_ACE.txt")
max_position = max(MI1H_REP2$position, na.rm = TRUE)
MI1H_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1H_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1H_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1H_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1H_rep_err <- list(MI1H_REP1, MI1H_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1H_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1H_rep_err.txt")
print(MI1H_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1H <- MI1H_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1H
# MI1I --------------------------------------------------------------------
core_title = "MI1I"
# Import data and remove validity = FALSE rows
MI1I_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1I_rep1/Results_ACE.txt")
max_position = max(MI1I_REP2$position, na.rm = TRUE)
MI1I_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1I_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1I_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1I_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1I_rep_err <- list(MI1I_REP1, MI1I_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1I_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1I_rep_err.txt")
print(MI1I_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1I <- MI1I_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1I
# MI1J --------------------------------------------------------------------
core_title = "MI1J"
# Import data and remove validity = FALSE rows
MI1J_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1J_rep1/Results_ACE.txt")
max_position = max(MI1J_REP2$position, na.rm = TRUE)
MI1J_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1J_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1J_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1J_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1J_rep_err <- list(MI1J_REP1, MI1J_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1J_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1J_rep_err.txt")
print(MI1J_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1J <- MI1J_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1J
# MI1K --------------------------------------------------------------------
core_title = "MI1K"
# Import data and remove validity = FALSE rows
MI1K_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1K_rep1/Results_ACE.txt")
max_position = max(MI1K_REP2$position, na.rm = TRUE)
MI1K_REP1 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1K_xrf/Results_ACE.txt") %>% 
  filter(position <= max_position) %>% 
  filter(!validity == FALSE) %>% 
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan1")
MI1K_REP2 <- itrax_import("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Input/Repeat_err/Marion/MI1K_rep1/Results_ACE.txt") %>%
  filter(!validity == FALSE) %>%
  select(any_of(c(elementsList, "position", "Mo inc", "Mo coh"))) %>%
  mutate(scan = "scan2")

MI1K_rep_err <- list(MI1K_REP1, MI1K_REP2) %>%
  reduce(full_join) %>%
  select(scan, position, everything()) %>%
  group_by(position) %>%
  summarise(across(any_of(c(elementsList, "Mo inc", "Mo coh")), 
                   function(x){set_errors(x = mean(x, na.rm = TRUE), 
                                          value = sd(x, na.rm = TRUE))}))
print(MI1K_rep_err, n = 50)
sink(file = "Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/MI1K_rep_err.txt")
print(MI1K_rep_err, n = 50)
sink(file = NULL)

# print summary graph
MI1K <- MI1K_rep_err %>%
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  select(elementsList, `Mo inc`, `Mo coh`, `coh/inc`, position) %>%
  pivot_longer(!c("position"), names_to = "elements", values_to = "peakarea") %>% 
  drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "Mo inc", "Mo coh", "coh/inc"))) %>%
  ggplot(aes(x = peakarea, y = position)) +
  scale_y_reverse() + 
  geom_ribbonh(aes(xmin = errors_min(peakarea), xmax = errors_max(peakarea)), fill = "grey80") +
  geom_lineh() + 
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(vars(elements), scales = "free_x", nrow = 1) + 
  ggtitle(paste(core_title)) +
  theme_paleo()
MI1K
 # rep1 missing from currently missing from dataset

# Summary Figure ----------------------------------------------------------
ggarrange(MI1C, MI1F, MI1J, ncol = 1, 
          nrow = 3, labels = c("a", "b", "c"))
ggsave("Papers_R/2024_DeVleeschouwer/ACE_Errors/Data/Output/Repeat_err/Marion/Marion_rep_err.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")



