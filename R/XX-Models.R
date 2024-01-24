
library(tidyverse)
library(tidybayes)
library(brms)

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T)))))
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  mutate(across(c(Born:LRS, AgeAccel), list(sc = function(x) as.vector(scale(x, center = T)))))
# Limit life history data to only reliable sampling dates
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)
# Load aging data
epi_dat <- readRDS('input/PB_clock_ages2.rds') %>%
  rename('Year' = yr) %>%
  mutate(Born = Year - round(Age),
         across(c(Year, AgeAccel, Age, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Model epigenetic acceleration ~ birth year
accel_born_mod <-  brm(AgeAccel_sc ~ Born_sc + Sex + (Born_sc + Sex | BearID),
                       data = epi_dat, family = gaussian, 
                       iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend = 'cmdstanr')

# Model epigenetic acceleration ~ first repro
accel_fr_mod <- brm(AgeAccel_sc ~ FirstRepro_sc + Sex + (FirstRepro_sc | BearID),
                    data = lh_epi_dat, family = gaussian, 
                    iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                    prior = prior(normal(0,1), class = b),
                    control = list(adapt_delta = 0.99, max_treedepth = 18),
                    backend = 'cmdstanr')

# Model age at first repro ~ birth year
repro_year_mod <- brm(FirstRepro_sc ~ Born_sc*Sex,
                      data = lh_pop_dat_trunc, family = gaussian, 
                      iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.99, max_treedepth = 18),
                      backend = 'cmdstanr')

# Model lifetime reproductive success ~ age at first reproduction
lrs_fr_mod <-  brm(LRS ~ FirstRepro_sc*Born_sc + Sex,
                      data = lh_pop_dat_trunc, family = negbinomial(link = 'log', link_shape = 'log'), 
                      iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.99, max_treedepth = 20),
                      backend = 'cmdstanr')

# Save models
saveRDS(accel_born_mod, 'models/accel_born_mod.rds')
saveRDS(accel_fr_mod, 'models/accel_fr_mod.rds')
saveRDS(repro_year_mod, 'models/repro_yr_mod.rds')
saveRDS(lrs_fr_mod, 'models/lrs_fr_mod.rds')
