
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds')
lh_epi_dat <- readRDS('output/lh_info_epi.rds')
epi_dat <- readRDS('input/PB_clock_ages2.rds')

# Model age at first repro ~ birth year
repro_year_mod <- brm(FirstRepro ~ Born*Sex,
           data = lh_pop_dat, family = gaussian, 
           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
           prior = prior(normal(5,1), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           backend = 'cmdstanr')

# Model epigenetic acceleration ~ first repro
accel_fr_mod <- brm(AgeAccel ~ FirstRepro*Sex + (FirstRepro*Sex | BearID),
           data = lh_epi_dat, family = gaussian, 
           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
           prior = prior(normal(5,1), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           backend = 'cmdstanr')

# Model epigenetic acceleration ~ year
# *** Note: small number of divergent transitions. Otherwise seems ok.
accel_year_mod <-  brm(AgeAccel ~ yr + (yr | BearID),
                     data = epi_dat, family = gaussian, 
                     iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                     prior = prior(normal(5,1), class = b),
                     control = list(adapt_delta = 0.99, max_treedepth = 20),
                     backend = 'cmdstanr')

# Get predicted draws from posteriors for plotting
accel_fr_draws <- lh_epi_dat |>
  data_grid(FirstRepro = seq_range(FirstRepro, n = 1000), Sex = c('M', 'F'), BearID = NA)  |>
  add_epred_draws(accel_fr_mod, ndraws = 1000)
repro_year_draws <- lh_pop_dat |>
  data_grid(Born = seq_range(Born, n = 1000), Sex = c('M', 'F'))  |>
  add_epred_draws(repro_year_mod, ndraws = 1000)
accel_year_draws <- epi_dat |>
  data_grid(yr = seq_range(yr, n = 1000), BearID = NA)  |>
  add_epred_draws(accel_year_mod, ndraws = 1000)

# Save models
saveRDS(accel_fr_mod, 'models/accel_fr_mod.rds')
saveRDS(repro_year_mod, 'models/repro_yr_mod.rds')
saveRDS(accel_year_mod, 'models/accel_yr_mod.rds')

# Save posterior draws
saveRDS(accel_fr_draws, 'models/accel_fr_draws.rds')
saveRDS(repro_year_draws, 'models/repro_yr_draws.rds')
saveRDS(accel_year_draws, 'models/accel_yr_draws.rds')

