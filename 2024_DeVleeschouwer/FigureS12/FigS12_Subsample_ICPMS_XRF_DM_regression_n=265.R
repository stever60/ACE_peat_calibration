# Figure S12: Dry Mass vs coh/inc regression n = 265 

# Uses ACE_Subsample_ICPMS_XRF_matched dataset made in Figure 2
# n = 265 defined by size of ICPMS dataset
# reduced from n = 272 due to new subsample composite dataset

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
packages <- c('tidyverse', 'tidypaleo', 'dplyr', 'readr', 'ggpubr', 'patchwork',
              'gridExtra', 'cowplot', 'vegan', 'rioja', 'ellipse', 'factoextra',
              'reshape2', 'GGally', 'ggsci', 'ggdendro', 'dendextend', 'dynamicTreeCut',
              'colorspace', 'cluster', 'magrittr', 'mgcv', 'gtable', 'repr',
              'bestNormalize','sjmisc', 'chemometrics', 'compositions', 
              'RColorBrewer', 'ggsci', 'wesanderson', 'viridis', 
              'ggrepel', 'itraxR', 'PeriodicTable', 'errors', 'forecast', 'broom',
              'directlabels', 'performance', 'lmtest', 'ggpmisc', 'cowplot', 'Hmisc','car')
lapply(packages, library, character.only=TRUE)

# Dry Mass vs coh/inc regression n = 265 ---------------------------------------

ACE_Subsample_ICPMS_XRF_matched <- read_csv("Papers_R/2024_DeVleeschouwer/Figure2/Data/Output/ACE_Subsample_ICPMS_XRF_matched.csv") %>% 
  print()

# ACE OLS & WLS ----------------------------------------------------------------

# 1) Unweighted OLS (Ordinary Least Squares) - linear model & checks
ACE_DM_lm <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched)
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_DM_lm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_summary.txt")
summary(ACE_DM_lm)
glance(ACE_DM_lm)
model_performance(ACE_DM_lm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_DM_lm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 2) Weighted OLS linear model & checks
ACE_wlm <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weight = 1/(Dry_mass_err)^2)
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wlm) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_wt_summary.txt")
summary(ACE_wlm)
glance(ACE_wlm)
model_performance(ACE_wlm) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wlm) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_DM_lm) # Performance package summary check for heteroscedasticity
icc(ACE_DM_lm) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 3) Unweighted Linear Regression (WLS) model
ACE_model <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched) # define model
ACE_wt <- 1 / lm(abs(ACE_model$residuals) ~ ACE_model$fitted.values)$fitted.values^2 #define weights to use
ACE_wls <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weights=ACE_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls) # summary stats
glance(ACE_wls) # summary stats including AIC
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_summary.txt")
summary(ACE_wls)
glance(ACE_wls)
model_performance(ACE_wls) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls) # Performance package summary check for heteroscedasticity
icc(ACE_wls) # check for random efDMcts - returns NULL if none present
sink(file = NULL)

# 4) Weighted Linear Regression (WLS) - ICPMS error weighted  - model
ACE_model_wt <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weight = 1/Dry_mass_err^2) # define model
ACE_wt_wt <- 1 / lm(abs(ACE_model_wt$residuals) ~ ACE_model_wt$fitted.values)$fitted.values^2 #define weights to use
ACE_wls_wt <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weights=ACE_wt_wt) #perform weighted least squares regression
# Checks
summary(ACE_wls_wt) # summary stats
glance(ACE_wls_wt) # summary stats including AIC
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
check_model(ACE_wls_wt) # summary check plots
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_wt_performance.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Write summary stats/checks to txt file
sink(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_wt_summary.txt")
summary(ACE_wls_wt)
glance(ACE_wls_wt)
model_performance(ACE_wls_wt) # AIC - Akaike information criterion; BIC Bayesian IC: lower = better fit for both
bptest(ACE_wls_wt) # lmtest package check for heteroscedasticity - p <0.05 = reject null hypothesis - heteroscedasticity present
check_heteroscedasticity(ACE_wls_wt) # Performance package summary check for heteroscedasticity
icc(ACE_wls_wt) # check for random effects - returns NULL if none present
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
ACE_DM_lm_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_DM_lm_influential_names,] # outliers only using of index values
ACE_DM_lm_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_DM_lm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_DM_lm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_DM_lm_lev_bar.pdf")
barplot(hatvalues(ACE_DM_lm), 
        col = "aquamarine3")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_DM_lm_lev.pdf")
leveragePlots(ACE_DM_lm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_DM_lm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_DM_lm)
dev.off()

# 2) Weighted OLS linear model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wlm_hats <- as.data.frame(hatvalues(ACE_wlm))
ACE_wlm_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wlm_cooksD <- cooks.distance(ACE_wlm)
ACE_wlm_influential <- ACE_wlm_cooksD[(ACE_wlm_cooksD > (3 * mean(ACE_wlm_cooksD, na.rm = TRUE)))]
ACE_wlm_influential
ACE_wlm_influential_names <- names(ACE_wlm_influential)
ACE_wlm_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_wlm_influential_names,] # outliers only using of index values
ACE_wlm_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_wlm_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wlm_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_OLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wlm_lev_bar.pdf")
barplot(hatvalues(ACE_wlm), 
        col = "red")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wlm_lev.pdf")
leveragePlots(ACE_wlm,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wlm_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wlm)
dev.off()

# 3) Unweighted Linear Regression (WLS) model: Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_hats <- as.data.frame(hatvalues(ACE_wls))
ACE_wls_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_cooksD <- cooks.distance(ACE_wls)
ACE_wls_influential <- ACE_wls_cooksD[(ACE_wls_cooksD > (3 * mean(ACE_wls_cooksD, na.rm = TRUE)))]
ACE_wls_influential
ACE_wls_influential_names <- names(ACE_wls_influential)
ACE_wls_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_wls_influential_names,] # outliers only using of index values
ACE_wls_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_wls_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_lev_bar.pdf")
barplot(hatvalues(ACE_wls), 
        col = "blue")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_lev.pdf")
leveragePlots(ACE_wls,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls)
dev.off()

# 4) Weighted Linear Regression (WLS): Leverage & Cooks distance - influence tests
# Leverage - 2x difference from mean consider removal due to high leverage
ACE_wls_wt_hats <- as.data.frame(hatvalues(ACE_wls_wt))
ACE_wls_wt_hats
# Cooks distance - 2-3 x difference from mean 
ACE_wls_wt_cooksD <- cooks.distance(ACE_wls_wt)
ACE_wls_wt_influential <- ACE_wls_wt_cooksD[(ACE_wls_wt_cooksD > (3 * mean(ACE_wls_wt_cooksD, na.rm = TRUE)))]
ACE_wls_wt_influential
ACE_wls_wt_influential_names <- names(ACE_wls_wt_influential)
ACE_wls_wt_outliers <- ACE_Subsample_ICPMS_XRF_matched[ACE_wls_wt_influential_names,] # outliers only using of index values
ACE_wls_wt_no_outliers <- ACE_Subsample_ICPMS_XRF_matched %>% anti_join(ACE_wls_wt_outliers) # generates a new dataset with outliers removed
write.csv(ACE_wls_wt_no_outliers,"Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_WLS_wt_no_outliers.csv", row.names = FALSE)

#Summary Leverage & Cooks distnce plots - write to file 
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_wt_lev_bar.pdf")
barplot(hatvalues(ACE_wls_wt), 
        col = "green")
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_wt_lev.pdf")
leveragePlots(ACE_wls_wt,layout=c(2,2)) 
dev.off()
pdf(file = "Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/DM_wls_wt_lev_cooks.pdf")
par(mfrow = c(2, 2))
plot(ACE_wls_wt)
dev.off()

# Linear model summary elemental plots
site_title = "ACE"
element_title <- "coh_inc"
theme_set(theme_classic(10))
ACE <- ggplot(ACE_Subsample_ICPMS_XRF_matched, aes(x = coh_inc, y = Dry_mass)) +
  geom_errorbar(aes(ymin=Dry_mass-Dry_mass_err, ymax=Dry_mass+Dry_mass_err), width=0, colour = "grey", alpha = 0.7) +
  geom_errorbar(aes(xmin=coh_inc-coh_inc_sd, xmax=coh_inc+coh_inc_sd), width=0, colour = "grey", alpha = 0.7) +
  geom_point(aes(fill = Site, color = Site), size = 2) +
  scale_color_jco() +
  #geom_point(fill = "Site", color = "Site", size = 2) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", colour = "blue") +
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2,
  #            aes(weight = 1/Dry_mass_err^2), colour = "blue") +
  stat_poly_eq(formula=y~x, use_label(c("n")), colour = "black", label.y = "top", label.x = "right") + 
  #stat_poly_line(formula=y~x, colour = "red", linetype = "dashed") + # unweighted line
  #stat_poly_eq(formula=y~x, use_label(c("eq", "R2", "p")), colour = "red", label.y = 0.85, label.x = -Inf, hjust = -0.18) + # unweighted stats; also "R2.confint", "adj.R2",
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
              linetype = "solid", aes(weight = ACE_wt), colour="black") + # WLS regression
  #geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=FALSE, alpha = 0.2, 
  #            aes(weight = ACE_wt_wt), colour="black") + # weighted WLS regression
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
  #xlim(0.15, 0.4) +
  #ylim(0, 80) +
  ggtitle(paste("Site: ", site_title, ": Dry mass vs ", element_title))
ACE
# Define p value, OLS equation & R2 as a string to add to plot

ACE_DM_lm_p <- function(ACE_DM_lm) {
  f <- summary(ACE_DM_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_DM_lm_p(ACE_DM_lm)

ACE_DM_lm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_DM_lm_p(ACE_DM_lm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted OLS equation & R2 as a string to add to plot

ACE_wlm_p <- function(ACE_wlm) {
  f <- summary(ACE_wlm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wlm_p(ACE_wlm)

ACE_wlm_eqn <- function(df){
  m <- lm(Dry_mass ~ coh_inc, data = ACE_Subsample_ICPMS_XRF_matched, weight = 1/(Dry_mass_err)^2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wlm_p(ACE_wlm), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, WLS equation & R2 as a string to add to plot

ACE_wls_p <- function(ACE_wls) {
  f <- summary(ACE_wls)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_p(ACE_wls)

ACE_wls_eqn <- function(df){
  m <- ACE_wls;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_p(ACE_wls), digits = 2)))
  as.character(as.expression(eq));
}

# Define p value, weighted WLS equation & R2 as a string to add to plot

ACE_wls_wt_p <- function(ACE_wls_wt) {
  f <- summary(ACE_wls_wt)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
ACE_wls_wt_p(ACE_wls_wt)

ACE_wls_wt_eqn <- function(df){
  m <- ACE_wls_wt;
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(ACE_wls_wt_p(ACE_wls_wt), digits = 2)))
  as.character(as.expression(eq));
}

# Add weighted OLS & WLS R2, equation p-value to top right of plot
ACE_final <- ACE + 
  geom_text(data=data.frame(), aes(label=paste("WLS_wt: ", ACE_wls_wt_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.1, vjust = 2, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS_wt: ", ACE_wlm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.1, vjust = 3, size = 3) +
  geom_text(data=data.frame(), aes(label = paste("WLS: ", ACE_wls_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "black", hjust = -0.15, vjust = 4, size = 3) +
  geom_text(data=data.frame(), aes(label=paste("OLS: ", ACE_DM_lm_eqn(df)), x = -Inf, y = Inf),
            parse = TRUE, colour = "blue", hjust = -0.15, vjust = 5, size = 3)
ACE_final
ggsave("Papers_R/2024_DeVleeschouwer/FigureS12/Data/Output/OLS_WLS/ACE/FigS12_ACE_DM_OLS_WLS_summary.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")




