# Figure 3 - initial data and regression analysis

# a) LM tests - log inc

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

# Define elements to use ----------------------------------------------------------

# ICP elements as defined by Francois
icp_Elements_fdv <- c("P_ICP", "K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "As_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP", "Pb_ICP", "dry_mass_pc")

# Below were defined by autocorrelation (acf) analysis of ACE09 composite ITRAX dataframe analysis file - see other R files for mor info

# XRF elements defined by ITRAX acf and matched to Francois ICPMS element list above
acf_icp_Elements_min <- c("K", "Ca", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", 
                          "Zn", "Rb", "Sr", "Zr", "Mo_coh")
acf_icp_Elements_min_sd <- c("K_sd", "Ca_sd", "Ti_sd", "Mn_sd", "Fe_sd", "Co_sd", "Ni_sd", "Cu_sd", 
                             "Zn_sd", "Rb_sd", "Sr_sd", "Zr_sd", "Mo_coh_sd")

# ICP elements defined by Francois & ITRAX acf
icp_Elements_min <- c("K_ICP", "Ca_ICP", "Ti_ICP", "Mn_ICP", "Fe_ICP", "Co_ICP", "Ni_ICP", "Cu_ICP", 
                      "Zn_ICP", "Rb_ICP", "Sr_ICP", "Zr_ICP")
icp_Elements_min_sd <- c("K_ICP_sd", "Ca_ICP_sd", "Ti_ICP_sd", "Mn_ICP_sd", "Fe_ICP_sd", "Co_ICP_sd", "Ni_ICP_sd", "Cu_ICP_sd", 
                         "Zn_ICP_sd", "Rb_ICP_sd", "Sr_ICP_sd", "Zr_ICP_sd")

# Import original ACE inc normalised matched file and log transform ----------------------------------------------

ACE_xrf_icp_matched_input <-read_csv("Papers_R/2024_DeVleeschouwer/Figure3/Data/Input/log_inc/ACE_xrf_icp_matched_inc.csv")
is.na(ACE_xrf_icp_matched)<-sapply(ACE_xrf_icp_matched, is.infinite) # replace any infinite values with NA
ACE_xrf_icp_matched_input

# ICPMS data - change errors supplied by FDV originally to CRM % values in Table 2 
ACE_xrf_icp_matched <- ACE_xrf_icp_matched_input %>%
  rename(K_ICP_sd_FDV2023 = K_ICP_sd) %>% 
  mutate(K_ICP_sd = K_ICP*0.18) %>%
  relocate(K_ICP_sd, .after = K_ICP) %>%
  rename(Ca_ICP_sd_FDV2023 = Ca_ICP_sd) %>% 
  mutate(Ca_ICP_sd = Ca_ICP*0.13) %>% 
  relocate(Ca_ICP_sd, .after = Ca_ICP) %>%
  rename(Ti_ICP_sd_FDV2023 = Ti_ICP_sd) %>% 
  mutate(Ti_ICP_sd = Ti_ICP*0.14) %>%
  relocate(Ti_ICP_sd, .after = Ti_ICP) %>%
  rename(Mn_ICP_sd_FDV2023 = Mn_ICP_sd) %>% 
  mutate(Mn_ICP_sd = Mn_ICP*0.12) %>% 
  relocate(Mn_ICP_sd, .after = Mn_ICP) %>%
  rename(Fe_ICP_sd_FDV2023 = Fe_ICP_sd) %>% 
  mutate(Fe_ICP_sd = Fe_ICP*0.18) %>%
  relocate(Fe_ICP_sd, .after = Fe_ICP) %>%
  rename(Co_ICP_sd_FDV2023 = Co_ICP_sd) %>% 
  mutate(Co_ICP_sd = Co_ICP*0.07) %>%
  relocate(Co_ICP_sd, .after = Co_ICP) %>%
  rename(Ni_ICP_sd_FDV2023 = Ni_ICP_sd) %>% 
  mutate(Ni_ICP_sd = Ni_ICP*0.54) %>%
  relocate(Ni_ICP_sd, .after = Ni_ICP) %>%
  rename(Cu_ICP_sd_FDV2023 = Cu_ICP_sd) %>% 
  mutate(Cu_ICP_sd = Cu_ICP*0.14) %>%
  relocate(Cu_ICP_sd, .after = Cu_ICP) %>%
  rename(Zn_ICP_sd_FDV2023 = Zn_ICP_sd) %>% 
  mutate(Zn_ICP_sd = Zn_ICP*0.09) %>% 
  relocate(Zn_ICP_sd, .after = Zn_ICP) %>%
  rename(Rb_ICP_sd_FDV2023 = Rb_ICP_sd) %>% 
  mutate(Rb_ICP_sd = Rb_ICP*0.16) %>%
  relocate(Rb_ICP_sd, .after = Rb_ICP) %>%
  rename(Sr_ICP_sd_FDV2023 = Sr_ICP_sd) %>% 
  mutate(Sr_ICP_sd = Sr_ICP*0.02) %>% 
  relocate(Sr_ICP_sd, .after = Sr_ICP) %>% 
  rename(Zr_ICP_sd_FDV2023 = Zr_ICP_sd) %>% 
  mutate(Zr_ICP_sd = Zr_ICP*0.05) %>% 
  relocate(Zr_ICP_sd, .after = Zr_ICP) %>% 
  #calculate error in log inc ITRAX data
  mutate(K_sd = K*0.1) %>%
  relocate(K_sd, .after = K) %>%
  mutate(Ca_sd = Ca*0.1) %>% 
  relocate(Ca_sd, .after = Ca) %>%
  mutate(Ti_sd = Ti*0.1) %>%
  relocate(Ti_sd, .after = Ti) %>%
  mutate(Mn_sd = Mn*0.1) %>% 
  relocate(Mn_sd, .after = Mn) %>%
  mutate(Fe_sd = Fe*0.1) %>%
  relocate(Fe_sd, .after = Fe) %>%
  mutate(Co_sd = Co*0.1) %>%
  relocate(Co_sd, .after = Co) %>%
  mutate(Ni_sd = Ni*0.1) %>%
  relocate(Ni_sd, .after = Ni) %>%
  mutate(Cu_sd = Cu*0.1) %>%
  relocate(Cu_sd, .after = Cu) %>%
  mutate(Zn_sd = Zn*0.1) %>% 
  relocate(Zn_sd, .after = Zn) %>%
  mutate(Rb_sd = Rb*0.1) %>%
  relocate(Rb_sd, .after = Rb) %>%
  mutate(Sr_sd = Sr*0.1) %>% 
  relocate(Sr_sd, .after = Sr) %>% 
  mutate(Zr_sd = Zr*0.1) %>% 
  relocate(Zr_sd, .after = Zr) 
write.csv(ACE_xrf_icp_matched,"Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/log_inc/ACE_matched_xrf_icp_inc.csv", row.names = FALSE)

# Replace zeros with half min value then take natural log of all elements and element ratios
ACE_all <- ACE_xrf_icp_matched %>% 
  mutate(across(all_of(c(icp_Elements_min, icp_Elements_min_sd)), ~ ifelse(.x < 0, 0, .x))) %>%  # replace any ICPMS values <0 with zero
  mutate_at(vars(all_of(c(icp_Elements_min, acf_icp_Elements_min))), ## Replace zeros with half minimum value to allow linear modelling to work
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% # Recommended procedure from Bertrand et al. (submitted) - retains dataframe structure
  mutate_at(vars(all_of(c(icp_Elements_min_sd, acf_icp_Elements_min_sd))), 
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% 
  select(Location:MSE, all_of(acf_icp_Elements_min), coh_inc, all_of(acf_icp_Elements_min_sd), coh_inc_sd,
         all_of(icp_Elements_min), dry_mass_pc, all_of(icp_Elements_min_sd), dry_mass_err) %>% 
  mutate(across(all_of(c(acf_icp_Elements_min, acf_icp_Elements_min_sd, icp_Elements_min, icp_Elements_min_sd)), log)) %>% #log all xrf and icp data
  na.omit() %>% #remove rows with NAs - in this case there is only one at the end of CAT1-S1-1F
  print()
write.csv(ACE_all,"Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/log_inc/ACE_xrf_icp_matched_log_inc.csv", row.names = FALSE)

# Define dataset to use for Linear modelling  -----------------------------

ACE_LM1 <- ACE_all %>%
  filter(!Site =="POB4") %>% #remove POB4 data from ACE dataset - has different ITRAX run parameters
  print()
write.csv(ACE_LM1,"Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/log_inc/ACE_xrf_icp_matched_log_inc_noPOB4.csv", row.names = FALSE)

ACE_LM1_stats <- ACE_LM1 %>%
  select(kcps: dry_mass_err) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(ACE_LM1_stats,"Papers_R/2024_DeVleeschouwer/Figure3/Data/Output/log_inc/ACE_LM1_log_inc_stats.csv", row.names = FALSE)

# Convert Site to use as a grouping variable
ACE_LM1$Site <- as.factor(ACE_LM1$Site)

# Linear modelling - OLS regression tests ----------------------------------------------- 

# Define titles & labels for plotting
XRF_title <- "/ inc [Ln XRF-CS]"
ICP_title <- " [Ln ICPMS]"
correlation_title <- "Log Correlation"
palette_set <- "jco" # or "jco", or "npg", "uchicago"

# Primary Matched Elements shown in Figure 3 ----------------------------------------------------
# Ca ----------------------------------------------------------------------
element_title <- "Ca"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ca/Ca_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
Ca_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Ca, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Ca_XRF_violin_boxplot <- Ca_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
ggdraw(Ca_XRF_violin)

# b) ICPMS data
theme_set(theme_classic(10))
Ca_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Ca_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Ca_ICP_violin_boxplot <- Ca_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
ggdraw(Ca_ICP_violin)

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Ca_pmain <- ggplot(ACE_LM1, aes(x = Ca, y = Ca_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Ca, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Ca_ydens <- axis_canvas(Ca_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Ca_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Ca_p1 <- insert_xaxis_grob(Ca_pmain, Ca_xdens, grid::unit(.2, "null"), position = "top")
Ca_corr2 <- insert_yaxis_grob(Ca_p1, Ca_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Ca_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Ca_corr3 <- ggscatter(ACE_LM1, x = "Ca", y = "Ca_ICP", size = 1,
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
ggdraw(Ca_corr3)

# Fig3.2: 4-part plot  
ggarrange(Ca_XRF_violin_boxplot, Ca_ICP_violin_boxplot, Ca_corr2, Ca_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ca/Ca_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Fig3.3 - OLS Linear Model

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Ca
y1.reg <- ACE_LM1$Ca_ICP
x1 <- "Ca"
y1 <- "Ca_ICP"

# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)

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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ca/Ca_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggdraw(Ca_predict)
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ca/Ca_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3 - OLS tests - 4 part plot
ggarrange(Ca_predict,  Ca_qq.x1, Ca_res, Ca_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ca/Ca_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Ca, ACE_LM1$Ca_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Ca, Ca_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Ca, ACE_LM1$Ca_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ca/Ca_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Ti ----------------------------------------------------------------------
element_title <- "Ti"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ti/Ti_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
Ti_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Ti, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Ti_XRF_violin_boxplot <- Ti_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
ggdraw(Ti_XRF_violin)

# b) ICPMS data
theme_set(theme_classic(10))
Ti_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Ti_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Ti_ICP_violin_boxplot <- Ti_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
ggdraw(Ti_ICP_violin)

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Ti_pmain <- ggplot(ACE_LM1, aes(x = Ti, y = Ti_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Ti, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Ti_ydens <- axis_canvas(Ti_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Ti_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Ti_p1 <- insert_xaxis_grob(Ti_pmain, Ti_xdens, grid::unit(.2, "null"), position = "top")
Ti_corr2 <- insert_yaxis_grob(Ti_p1, Ti_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Ti_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Ti_corr3 <- ggscatter(ACE_LM1, x = "Ti", y = "Ti_ICP", size = 1,
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
ggdraw(Ti_corr3)

# Fig3.2: 4-part plot  
ggarrange(Ti_XRF_violin_boxplot, Ti_ICP_violin_boxplot, Ti_corr2, Ti_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ti/Ti_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests - Ti

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Ti
y1.reg <- ACE_LM1$Ti_ICP
x1 <- "Ti"

# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)

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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ti/Ti_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ti/Ti_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3 - 4 part figure
ggarrange(Ti_predict,  Ti_qq.x1, Ti_res, Ti_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ti/Ti_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Ti, ACE_LM1$Ti_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Ti, Ti_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Ti, ACE_LM1$Ti_ICP)
round(cov_LM1,2)

# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)

# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ti/Ti_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Mn ----------------------------------------------------------------------
element_title <- "Mn"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Mn/Mn_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS
# a) XRF-CS data distribution
theme_set(theme_classic(10))
Mn_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Mn, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Mn_XRF_violin_boxplot <- Mn_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Mn_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Mn_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Mn_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Mn_ICP_violin_boxplot <- Mn_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Mn_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Mn_pmain <- ggplot(ACE_LM1, aes(x = Mn, y = Mn_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Mn, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Mn_ydens <- axis_canvas(Mn_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Mn_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Mn_p1 <- insert_xaxis_grob(Mn_pmain, Mn_xdens, grid::unit(.2, "null"), position = "top")
Mn_corr2 <- insert_yaxis_grob(Mn_p1, Mn_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Mn_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Mn_corr3 <- ggscatter(ACE_LM1, x = "Mn", y = "Mn_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Mn_XRF_violin_boxplot, Mn_ICP_violin_boxplot, Mn_corr2, Mn_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Mn/Mn_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests - Mn

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Mn
y1.reg <- ACE_LM1$Mn_ICP
x1 <- "Mn"
y1 <- "Mn_ICP"
# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Mn/Mn_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Mn/Mn_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Mn_predict,  Mn_qq.x1, Mn_res, Mn_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Mn/Mn_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Mn, ACE_LM1$Mn_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Mn, Mn_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Mn, ACE_LM1$Mn_ICP)
round(cov_LM1,2)

# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Mn/Mn_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Fe ----------------------------------------------------------------------
element_title <- "Fe"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Fe/Fe_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS
# a) XRF-CS data distribution
theme_set(theme_classic(10))
Fe_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Fe, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Fe_XRF_violin_boxplot <- Fe_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Fe_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Fe_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Fe_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = , face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Fe_ICP_violin_boxplot <- Fe_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Fe_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Fe_pmain <- ggplot(ACE_LM1, aes(x = Fe, y = Fe_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Fe, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Fe_ydens <- axis_canvas(Fe_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Fe_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Fe_p1 <- insert_xaxis_grob(Fe_pmain, Fe_xdens, grid::unit(.2, "null"), position = "top")
Fe_corr2 <- insert_yaxis_grob(Fe_p1, Fe_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Fe_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Fe_corr3 <- ggscatter(ACE_LM1, x = "Fe", y = "Fe_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Fe_XRF_violin_boxplot, Fe_ICP_violin_boxplot, Fe_corr2, Fe_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Fe/Fe_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests - Fe

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Fe
y1.reg <- ACE_LM1$Fe_ICP
x1 <- "Fe"
y1 <- "Fe_ICP"

# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)

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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Fe/Fe_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Fe/Fe_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Fe_predict,  Fe_qq.x1, Fe_res, Fe_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Fe/Fe_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Fe, ACE_LM1$Fe_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Fe, Fe_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Fe, ACE_LM1$Fe_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Fe/Fe_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Sr ----------------------------------------------------------------------
element_title <- "Sr"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Sr/Sr_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
Sr_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Sr, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Sr_XRF_violin_boxplot <- Sr_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Sr_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Sr_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Sr_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Sr_ICP_violin_boxplot <- Sr_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Sr_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Sr_pmain <- ggplot(ACE_LM1, aes(x = Sr, y = Sr_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Sr, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Sr_ydens <- axis_canvas(Sr_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Sr_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Sr_p1 <- insert_xaxis_grob(Sr_pmain, Sr_xdens, grid::unit(.2, "null"), position = "top")
Sr_corr2 <- insert_yaxis_grob(Sr_p1, Sr_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Sr_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Sr_corr3 <- ggscatter(ACE_LM1, x = "Sr", y = "Sr_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Sr_XRF_violin_boxplot, Sr_ICP_violin_boxplot, Sr_corr2, Sr_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Sr/Sr_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Sr
y1.reg <- ACE_LM1$Sr_ICP
x1 <- "Sr"
y1 <- "Sr_ICP"
# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Sr/Sr_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Sr/Sr_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Sr_predict,  Sr_qq.x1, Sr_res, Sr_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Sr/Sr_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Sr, ACE_LM1$Sr_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Sr, Sr_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Sr, ACE_LM1$Sr_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Sr/Sr_Summary_stats.txt")
print(correlate_LM1)

print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Zr ----------------------------------------------------------------------
element_title <- "Zr"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zr/Zr_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
Zr_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Zr, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Zr_XRF_violin_boxplot <- Zr_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Zr_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Zr_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Zr_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Zr_ICP_violin_boxplot <- Zr_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Zr_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Zr_pmain <- ggplot(ACE_LM1, aes(x = Zr, y = Zr_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Zr, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Zr_ydens <- axis_canvas(Zr_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Zr_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Zr_p1 <- insert_xaxis_grob(Zr_pmain, Zr_xdens, grid::unit(.2, "null"), position = "top")
Zr_corr2 <- insert_yaxis_grob(Zr_p1, Zr_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Zr_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Zr_corr3 <- ggscatter(ACE_LM1, x = "Zr", y = "Zr_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Zr_XRF_violin_boxplot, Zr_ICP_violin_boxplot, Zr_corr2, Zr_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zr/Zr_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Zr
y1.reg <- ACE_LM1$Zr_ICP
x1 <- "Zr"
y1 <- "Zr_ICP"
# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zr/Zr_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
Zr_predict
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zr/Zr_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Zr_predict,  Zr_qq.x1, Zr_res, Zr_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zr/Zr_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Zr, ACE_LM1$Zr_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Zr, Zr_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Zr, ACE_LM1$Zr_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zr/Zr_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)




# Secondary Matched Elements shown in Supplementary Material ----------------------------------------------
# K -----------------------------------------------------------------------
element_title <- "K"




# Fig3.1 - Correlation plot with 68% CI ellipses 
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
  ggtitle(paste("ACE (OLS): Ln ", element_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=10, face="bold"))
K_corr1
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/K/K_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
K_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=K, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
K_XRF_violin_boxplot <- K_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
ggdraw(K_XRF_violin)

# b) ICPMS data
theme_set(theme_classic(10))
K_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=K_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
K_ICP_violin_boxplot <- K_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
ggdraw(K_ICP_violin)

# c) Scatter & data density plot
theme_set(theme_classic(10))
K_pmain <- ggplot(ACE_LM1, aes(x = K, y = K_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = K, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
K_ydens <- axis_canvas(K_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = K_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
K_p1 <- insert_xaxis_grob(K_pmain, K_xdens, grid::unit(.2, "null"), position = "top")
K_corr2 <- insert_yaxis_grob(K_p1, K_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(K_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
K_corr3 <- ggscatter(ACE_LM1, x = "K", y = "K_ICP", size = 1,
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
ggdraw(K_corr3)

# Fig3.2: 4-part plot  
ggarrange(K_XRF_violin_boxplot, K_ICP_violin_boxplot, K_corr2, K_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/K/K_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Build Linear OLS Set up x and y variables
x1.reg <- ACE_LM1$K
y1.reg <- ACE_LM1$K_ICP
x1 <- "K"
y1 <- "K_ICP"

# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)

# Look at output
model_1
summary(model_1)
confint(model_1)

# Create predicted values and upper/lower CI to check model is working OK
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")
# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- bind_cols(ACE_LM1, pred.int1)
data_1_out <- bind_cols(ACE_LM1, pred.int1) %>%
  select(c(Location:midpoint, K, K_sd, K_ICP, K_ICP_sd, fit, lwr, upr)) %>% 
  print()
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/K/K_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
theme_set(theme_classic(base_size=10))
options(scipen = 999)

# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/K/K_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3 - OLS tests - 4 part plot
ggarrange(K_predict,  K_qq.x1, K_res, K_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/K/K_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$K, ACE_LM1$K_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(K, K_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$K, ACE_LM1$K_ICP)
round(cov_LM1,2)

# Shapiro-Wilk test all models
# Test normality of one variable (univariate) at a time 
# If the p-value > 0.05 it implies that the distribution of the data 
# are not significantly different from normal distribution - i.e., can assume normality
# < 0.05 means data are not normally distributed
# Normal distribution is hard to achieve in  geochemical datasets!
shapiro.test(x1.reg)
shapiro.test(y1.reg)

# Durbin-Watson Test
# Use  to test for correlation among residuals
# H0 (null hypothesis) = no correlation among the residuals
# HA (alternative hypothesis) = residuals are autocorrelated
# Value of 2.0 indicates there is no autocorrelation detected in the sample. 
# 0-2 = positive autocorrelation; 2-4 = negative autocorrelation.
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/K/K_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Co ----------------------------------------------------------------------
element_title <- "Co"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Co/Co_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS
# a) XRF-CS data distribution
theme_set(theme_classic(10))
Co_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Co, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Co_XRF_violin_boxplot <- Co_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Co_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Co_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Co_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Co_ICP_violin_boxplot <- Co_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Co_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Co_pmain <- ggplot(ACE_LM1, aes(x = Co, y = Co_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Co, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Co_ydens <- axis_canvas(Co_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Co_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Co_p1 <- insert_xaxis_grob(Co_pmain, Co_xdens, grid::unit(.2, "null"), position = "top")
Co_corr2 <- insert_yaxis_grob(Co_p1, Co_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Co_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Co_corr3 <- ggscatter(ACE_LM1, x = "Co", y = "Co_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Co_XRF_violin_boxplot, Co_ICP_violin_boxplot, Co_corr2, Co_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Co/Co_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Co
y1.reg <- ACE_LM1$Co_ICP
x1 <- "Co"
y1 <- "Co_ICP"
# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Co/Co_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Co/Co_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Co_predict,  Co_qq.x1, Co_res, Co_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Co/Co_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Co, ACE_LM1$Co_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Co, Co_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Co, ACE_LM1$Co_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Co/Co_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Ni ----------------------------------------------------------------------
element_title <- "Ni"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ni/Ni_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS
# a) XRF-CS data distribution
theme_set(theme_classic(10))
Ni_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Ni, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Ni_XRF_violin_boxplot <- Ni_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Ni_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Ni_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Ni_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Ni_ICP_violin_boxplot <- Ni_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Ni_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Ni_pmain <- ggplot(ACE_LM1, aes(x = Ni, y = Ni_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Ni, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Ni_ydens <- axis_canvas(Ni_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Ni_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Ni_p1 <- insert_xaxis_grob(Ni_pmain, Ni_xdens, grid::unit(.2, "null"), position = "top")
Ni_corr2 <- insert_yaxis_grob(Ni_p1, Ni_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Ni_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Ni_corr3 <- ggscatter(ACE_LM1, x = "Ni", y = "Ni_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Ni_XRF_violin_boxplot, Ni_ICP_violin_boxplot, Ni_corr2, Ni_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ni/Ni_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Ni
y1.reg <- ACE_LM1$Ni_ICP
x1 <- "Ni"
y1 <- "Ni_ICP"

# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ni/Ni_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Ni_predict,  Ni_qq.x1, Ni_res, Ni_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ni/Ni_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Ni, ACE_LM1$Ni_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Ni, Ni_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Ni, ACE_LM1$Ni_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Ni/Ni_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Cu ----------------------------------------------------------------------
element_title <- "Cu"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Cu/Cu_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS
# a) XRF-CS data distribution
theme_set(theme_classic(10))
Cu_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Cu, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Cu_XRF_violin_boxplot <- Cu_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Cu_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Cu_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Cu_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Cu_ICP_violin_boxplot <- Cu_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Cu_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Cu_pmain <- ggplot(ACE_LM1, aes(x = Cu, y = Cu_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Cu, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Cu_ydens <- axis_canvas(Cu_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Cu_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Cu_p1 <- insert_xaxis_grob(Cu_pmain, Cu_xdens, grid::unit(.2, "null"), position = "top")
Cu_corr2 <- insert_yaxis_grob(Cu_p1, Cu_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Cu_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Cu_corr3 <- ggscatter(ACE_LM1, x = "Cu", y = "Cu_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Cu_XRF_violin_boxplot, Cu_ICP_violin_boxplot, Cu_corr2, Cu_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Cu/Cu_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Cu
y1.reg <- ACE_LM1$Cu_ICP
x1 <- "Cu"
y1 <- "Cu_ICP"

# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Cu/Cu_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Cu/Cu_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Cu_predict,  Cu_qq.x1, Cu_res, Cu_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Cu/Cu_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Cu, ACE_LM1$Cu_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Cu, Cu_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Cu, ACE_LM1$Cu_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test 
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Cu/Cu_Summary_stats.txt")
print(correlate_LM1)

print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Zn ----------------------------------------------------------------------
element_title <- "Zn"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zn/Zn_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
Zn_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Zn, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Zn_XRF_violin_boxplot <- Zn_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Zn_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Zn_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Zn_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Zn_ICP_violin_boxplot <- Zn_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Zn_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Zn_pmain <- ggplot(ACE_LM1, aes(x = Zn, y = Zn_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Zn, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Zn_ydens <- axis_canvas(Zn_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Zn_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Zn_p1 <- insert_xaxis_grob(Zn_pmain, Zn_xdens, grid::unit(.2, "null"), position = "top")
Zn_corr2 <- insert_yaxis_grob(Zn_p1, Zn_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Zn_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Zn_corr3 <- ggscatter(ACE_LM1, x = "Zn", y = "Zn_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Zn_XRF_violin_boxplot, Zn_ICP_violin_boxplot, Zn_corr2, Zn_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zn/Zn_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Zn
y1.reg <- ACE_LM1$Zn_ICP
x1 <- "Zn"
y1 <- "Zn_ICP"
# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zn/Zn_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zn/Zn_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Zn_predict,  Zn_qq.x1, Zn_res, Zn_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zn/Zn_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Zn, ACE_LM1$Zn_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Zn, Zn_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Zn, ACE_LM1$Zn_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Zn/Zn_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Rb ----------------------------------------------------------------------
element_title <- "Rb"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Rb/Rb_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
Rb_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=Rb, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
Rb_XRF_violin_boxplot <- Rb_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
Rb_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
Rb_ICP_violin <- ggplot(ACE_LM1, aes(x=Site, y=Rb_ICP, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
Rb_ICP_violin_boxplot <- Rb_ICP_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
Rb_ICP_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
Rb_pmain <- ggplot(ACE_LM1, aes(x = Rb, y = Rb_ICP, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = Rb, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
Rb_ydens <- axis_canvas(Rb_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = Rb_ICP, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
Rb_p1 <- insert_xaxis_grob(Rb_pmain, Rb_xdens, grid::unit(.2, "null"), position = "top")
Rb_corr2 <- insert_yaxis_grob(Rb_p1, Rb_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(Rb_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
Rb_corr3 <- ggscatter(ACE_LM1, x = "Rb", y = "Rb_ICP", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(Rb_XRF_violin_boxplot, Rb_ICP_violin_boxplot, Rb_corr2, Rb_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Rb/Rb_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$Rb
y1.reg <- ACE_LM1$Rb_ICP
x1 <- "Rb"
y1 <- "Rb_ICP"
# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Rb/Rb_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Rb/Rb_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(Rb_predict,  Rb_qq.x1, Rb_res, Rb_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Rb/Rb_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$Rb, ACE_LM1$Rb_ICP)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(Rb, Rb_ICP))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$Rb, ACE_LM1$Rb_ICP)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/Rb/Rb_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# Dry mass (DM) ----------------------------------------------------------------------
element_title <- "DM"

# Fig3.1 - Correlation plot
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/DM/DM_Fig3.1_Correlation_OLS_all_sites.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.2 - Violin + box plots, scatterplot & density, Individual site LM-OLS

# a) XRF-CS data distribution
theme_set(theme_classic(10))
DM_XRF_violin <- ggplot(ACE_LM1, aes(x=Site, y=coh_inc, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, XRF_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
DM_XRF_violin_boxplot <- DM_XRF_violin + geom_boxplot(width=0.1, fill="white") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  ggtitle(paste("ACE:", element_title, XRF_title))
DM_XRF_violin

# b) ICPMS data
theme_set(theme_classic(10))
DM_pc_violin <- ggplot(ACE_LM1, aes(x=Site, y=dry_mass_pc, fill=Site)) + 
  geom_violin(trim=FALSE) + 
  ggpubr::fill_palette(palette_set) +
  labs(y = paste(element_title, ICP_title)) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10), plot.title = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8, face="italic"), 
        legend.justification = c(0, 1), legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE))
DM_pc_violin_boxplot <- DM_pc_violin + geom_boxplot(width=0.1, fill="white") +
  ggtitle(paste("ACE:", element_title, ICP_title))
DM_pc_violin

# c) Scatterplot and density plot
theme_set(theme_classic(10))
DM_pmain <- ggplot(ACE_LM1, aes(x = coh_inc, y = dry_mass_pc, color = Site)) +
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
  geom_density(data = ACE_LM1, aes(x = coh_inc, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  ggpubr::fill_palette(palette_set)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
DM_ydens <- axis_canvas(DM_pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ACE_LM1, aes(x = dry_mass_pc, fill = Site),
               alpha = 0.7, linewidth = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette_set)
DM_p1 <- insert_xaxis_grob(DM_pmain, DM_xdens, grid::unit(.2, "null"), position = "top")
DM_corr2 <- insert_yaxis_grob(DM_p1, DM_ydens, grid::unit(.2, "null"), position = "right")
ggdraw(DM_corr2)

# d) GLM-OLS per site
theme_set(theme_classic(10))
DM_corr3 <- ggscatter(ACE_LM1, x = "coh_inc", y = "dry_mass_pc", size = 1,
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

# Fig3.2: 4-part plot  
ggarrange(DM_XRF_violin_boxplot, DM_pc_violin_boxplot, DM_corr2, DM_corr3, ncol = 2, 
          nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/DM/DM_Fig3.2_Data_OLS_Sites.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Figure 3 - Linear Model (OLS) initial tests

# Linear OLS Model - Set up x and y variables
x1.reg <- ACE_LM1$coh_inc
y1.reg <- ACE_LM1$dry_mass_pc
x1 <- "coh_inc"
y1 <- "dry_mass_pc"
# Build linear regression model
model_1 <- lm(y1.reg~x1.reg, data = ACE_LM1)
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
write.csv(data_1_out,"Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/DM/DM_Model_predict.csv", row.names = FALSE)

# LM-OLS & Q-Q plots and residuals plot to assess normality
library(ggpubr)
theme_set(theme_classic(base_size=10))
# Fig3.3a - OLS Linear model
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
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/DM/DM_Fig3.3a_LM-OLS_predict.pdf", 
       height = c(16), width = c(16), dpi = 600, units = "cm")

# Fig3.3b - Residuals vs fitted values plot
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

#Fig3.3c, d
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

# Fig3.3
ggarrange(DM_predict,  DM_qq.x1, DM_res, DM_qq.y1,
          ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
ggsave("Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/DM/DM_Fig3.3_OLS_tests.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Statistics tests

# Hypothesis test 
ctest_LM1 <- cor.test(ACE_LM1$coh_inc, ACE_LM1$dry_mass_pc)
# Correlation stats by site
correlate_LM1 <- ACE_LM1 %>%
  group_by(Site) %>% 
  summarise(r = cor(coh_inc, dry_mass_pc))
# Examine co-variance in all models
cov_LM1 <- cov(ACE_LM1$coh_inc, ACE_LM1$dry_mass_pc)
round(cov_LM1,2)
# Shapiro-Wilk test all models
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# Durbin-Watson Test
dwtest(model_1)

# Write summary stats/normality & residual checks output to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/Figure3/Plots/Fig3_LM_Tests_log_inc/DM/DM_Summary_stats.txt")
print(correlate_LM1)
print(ctest_LM1)
print(summary(model_1))
print(confint(model_1))
print(dwtest(model_1))
print(shapiro.test(x1.reg))
print(shapiro.test(y1.reg))
sink(file = NULL)

# END ---------------------------------------------------------------------










