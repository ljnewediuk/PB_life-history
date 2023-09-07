
library(tidyverse)
library(tidybayes)
library(brms)

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T)))))
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  mutate(across(c(Born:LRS, AgeAccel), list(sc = function(x) as.vector(scale(x, center = T)))))
# Limit life history data to only reliable sampling dates
lh_epi_dat_trunc <- lh_epi_dat %>%
  filter(Born %in% 1980:2000)
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)
# Load aging data
epi_dat <- readRDS('input/PB_clock_ages2.rds') %>%
  rename('Year' = yr) %>%
  mutate(across(c(Year, AgeAccel), list(sc = function(x) as.vector(scale(x, center = T)))))

# Load growth rates and bind to age acceleration
# Load estimated growth rates
growth_rates <- readRDS('models/size_age_id_draws_F.rds') %>%
  rbind(readRDS('models/size_age_id_draws_M.rds')) %>%
  right_join(epi_dat) %>%
  mutate(across(c(mean_BearSlope, AgeAccel), list(sc = function(x) as.vector(scale(x, center = T)))))

# Load ice breakup dates and surface temperature data
# Ice data
ice_dat <- read.csv('input/WH_ice_breakup_dates.csv') %>%
  mutate(IceFree = jday_freezeup - jday_breakup) %>%
  rename('Year' = yr) %>%
  mutate(across(c(Year, IceFree), list(sc = function(x) as.vector(scale(x, center = T)))))
# Temperature data
temp_dat <- read.table('input/global_temp_data.txt', header = T) %>%
  rename('Temp' = No_Smoothing) %>%
  mutate(across(c(Year, Temp), list(sc = function(x) as.vector(scale(x, center = T)))))
  filter(Year %in% c(1979:2018))

# Model ice-free days ~ time
ice_free_mod <- brm(IceFree_sc ~ Year_sc,
                      data = ice_dat, family = gaussian, 
                      iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.99, max_treedepth = 18),
                      backend = 'cmdstanr')

# Model temperature ~ time
temp_mod <- brm(Temp_sc ~ Year_sc,
                data = temp_dat, family = gaussian, 
                iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                prior = prior(normal(0,1), class = b),
                control = list(adapt_delta = 0.99, max_treedepth = 18),
                backend = 'cmdstanr')

# Model age at first repro ~ birth year
repro_year_mod <- brm(FirstRepro_sc ~ Born_sc*Sex,
           data = lh_pop_dat, family = gaussian, 
           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
           prior = prior(normal(0,1), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           backend = 'cmdstanr')

# Model epigenetic acceleration ~ first repro
accel_fr_mod <- brm(AgeAccel_sc ~ FirstRepro_sc*Sex + (1 | BearID),
           data = lh_epi_dat, family = gaussian, 
           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
           prior = prior(normal(0,1), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           backend = 'cmdstanr')

# Model epigenetic acceleration ~ year
# *** Note: small number of divergent transitions. Otherwise seems ok.
accel_year_mod <-  brm(AgeAccel_sc ~ Year_sc + (Year_sc | BearID),
                     data = epi_dat, family = gaussian, 
                     iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                     prior = prior(normal(0,1), class = b),
                     control = list(adapt_delta = 0.99, max_treedepth = 20),
                     backend = 'cmdstanr')

# Model growth rates ~ epigenetic acceleration
accel_gr_mod <- brm(AgeAccel_sc ~ mean_BearSlope_sc*Sex + (1 | BearID),
                    data = growth_rates, family = gaussian, 
                    iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                    prior = prior(normal(0,1), class = b),
                    control = list(adapt_delta = 0.99, max_treedepth = 18),
                    backend = 'cmdstanr')

# Model lifetime reproductive success ~ epigenetic acceleration
lrs_accel_mod <-  brm(LRS ~ AgeAccel,
                       data = lh_epi_dat_trunc, family = negbinomial(link = 'log', link_shape = 'log'), 
                       iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend = 'cmdstanr')

# Model lifetime reproductive success ~ age at first reproduction
lrs_fr_mod <-  brm(LRS ~ FirstRepro_sc*Born_sc,
                      data = lh_pop_dat_trunc, family = negbinomial(link = 'log', link_shape = 'log'), 
                      iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.99, max_treedepth = 20),
                      backend = 'cmdstanr')

# Model lifetime reproductive success ~ year of birth
lrs_yr_mod <-  brm(LRS ~ Born,
                      data = lh_pop_dat_trunc, family = negbinomial(link = 'log', link_shape = 'log'), 
                      iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.99, max_treedepth = 20),
                      backend = 'cmdstanr')

# Save models
saveRDS(temp_mod, 'models/temp_mod.rds')
saveRDS(ice_free_mod, 'models/ice_free_mod.rds')
saveRDS(accel_fr_mod, 'models/accel_fr_mod.rds')
saveRDS(repro_year_mod, 'models/repro_yr_mod.rds')
saveRDS(accel_year_mod, 'models/accel_yr_mod.rds')
saveRDS(accel_gr_mod, 'models/accel_gr_mod.rds')
saveRDS(lrs_accel_mod, 'models/lrs_accel_mod.rds')
saveRDS(lrs_fr_mod, 'models/lrs_fr_mod.rds')
saveRDS(lrs_yr_mod, 'models/lrs_yr_mod.rds')
