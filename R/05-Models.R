
# 05 - Models ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Fit Bayesian GLMs to test hypotheses in MS
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(tidybayes)
library(brms)

# Fit models for MS:
# 
#   2) Age acceleration ~ birth year
#   3) Age accel ~ birth year separated by sex
#       a) Males
#       b) Females
#   4) Age acceleration ~ age at first reproduction
#   5) Lifetime reproduction ~ age at first reproduction

# 1 Load data ====

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T)))))

# With epigenetic data
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  mutate(across(c(Born:LRS, AgeAccel), list(sc = function(x) as.vector(scale(x, center = T)))))

# Limit life history data to only reliable sampling dates
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)

# Load aging data
epi_dat <- readRDS('output/PB_clock_ages.rds') %>%
  rename('Year' = yr) %>%
  mutate(Born = Year - round(Age),
         across(c(Year, AgeAccel, Age, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# For males and females separately
epi_dat_M <- epi_dat %>%
  filter(Sex == 'M')
epi_dat_F <- epi_dat %>%
  filter(Sex == 'F')

# 2 Model epigenetic acceleration ~ birth year ====
accel_born_mod <-  brm(AgeAccel_sc ~ Born_sc + Sex + (Born_sc + Sex | BearID),
                       data = epi_dat, family = gaussian, 
                       iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend = 'cmdstanr')

# 3a Model epigenetic acceleration ~ birth year for males ====
accel_born_mod_M <-  brm(AgeAccel_sc ~ Born_sc + (Born_sc | BearID),
                         data = epi_dat_M, family = gaussian, 
                         iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                         prior = prior(normal(0,1), class = b),
                         control = list(adapt_delta = 0.99, max_treedepth = 20),
                         backend = 'cmdstanr')

# 4b Model epigenetic acceleration ~ birth year for females ====
accel_born_mod_F <-  brm(AgeAccel_sc ~ Born_sc + (Born_sc | BearID),
                         data = epi_dat_F, family = gaussian, 
                         iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend = 'cmdstanr')

# 4 Model epigenetic acceleration ~ first repro ====
accel_fr_mod <- brm(AgeAccel_sc ~ FirstRepro_sc + Sex + (FirstRepro_sc | BearID),
                    data = lh_epi_dat, family = gaussian, 
                    iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                    prior = prior(normal(0,1), class = b),
                    control = list(adapt_delta = 0.99, max_treedepth = 18),
                    backend = 'cmdstanr')

# 5 Model lifetime reproductive success ~ age at first reproduction ====
lrs_fr_mod <-  brm(LRS ~ FirstRepro_sc*Born_sc + Sex,
                      data = lh_pop_dat_trunc, family = negbinomial(link = 'log', link_shape = 'log'), 
                      iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.99, max_treedepth = 20),
                      backend = 'cmdstanr')

# Save models
saveRDS(accel_born_mod, 'models/accel_born_mod.rds')
saveRDS(accel_fr_mod, 'models/accel_fr_mod.rds')
saveRDS(accel_born_mod_M, 'models/accel_born_mod_M.rds')
saveRDS(accel_born_mod_F, 'models/accel_born_mod_F.rds')
saveRDS(lrs_fr_mod, 'models/lrs_fr_mod.rds')
