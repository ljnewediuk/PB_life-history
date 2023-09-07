
library(tidyverse)
library(tidybayes)
library(brms)

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  mutate(YearFirstRepro = Born + FirstRepro) %>%
  filter(Born %in% 1980:2000)

sel_mod <- function(dat, year, sex) {
  
  sel_dat <- dat %>%
    filter(YearFirstRepro == year & Sex == sex)
  
  sel_lm <- lm(LRS ~ FirstRepro, data = sel_dat)
  
  sel_coeff <- as.numeric(sel_lm$coefficients[2])
  sel_sigma <- as.numeric(sqrt(diag(vcov(sel_lm)))[2])
  
  return(c(sel_coeff, sel_sigma))
  
}

male_selection <- data.frame()
for(yr in unique(lh_pop_dat$YearFirstRepro)) {
  
  beta_afr <- sel_mod(lh_pop_dat, yr, 'M')
  beta_afr_yr <- data.frame(Year = yr, Beta = beta_afr[1], Sigma = beta_afr[2])
  male_selection <- rbind(male_selection, beta_afr_yr)
  
}

female_selection <- data.frame()
for(yr in unique(lh_pop_dat$YearFirstRepro)) {
  
  beta_afr <- sel_mod(lh_pop_dat, yr, 'F')
  beta_afr_yr <- data.frame(Year = yr, Beta = beta_afr[1], Sigma = beta_afr[2])
  female_selection <- rbind(female_selection, beta_afr_yr)
  
}

male_selection_no <- male_selection %>%
  filter(Beta > -3 & Beta < 3)

ggplot() +
  geom_point(data = male_selection_no, aes(x = Year, y = Beta)) +
  geom_errorbar(data = male_selection_no, aes(x = Year, ymin = Beta - Sigma, ymax = Beta + Sigma))

ggplot() +
  geom_point(data = female_selection, aes(x = Year, y = Beta)) +
  geom_errorbar(data = female_selection, aes(x = Year, ymin = Beta - Sigma, ymax = Beta + Sigma))





