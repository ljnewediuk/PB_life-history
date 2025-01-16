
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
#   1) Age acceleration ~ birth year
#   2) Age acceleration ~ age at first reproduction
#   3) Lifetime reproduction ~ age at first reproduction

# 1 Load data ====

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  filter(! is.infinite(FirstRepro) | ! is.infinite(LastRepro)) %>%
  # Scale and centre variables
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T))))) %>%
  # Filter individuals born after 2000 (we might not have captured their full
  # reproductive lifespan) and before 1980 (not enough comprehensive sampling)
  filter(Born <= 2000 & Born >=1980)

# With epigenetic data
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  # Summarize by individual
  group_by(BearID) %>%
  summarize(AgeAccel = mean(AgeAccel), FirstRepro = mean(FirstRepro)) %>%
  # Scale and centre variables
  filter(! is.infinite(FirstRepro)) %>%
  mutate(across(c(AgeAccel, FirstRepro), list(sc = function(x) as.vector(scale(x, center = T)))))

# Load aging data
epi_dat <- readRDS('output/WH_combined_ages.rds') %>%
  # Summarize by individual
  group_by(BearID) %>%
  summarize(AgeAccel = mean(AgeAccel), Born = mean(Born)) %>%
  # Scale and centre variables
  mutate(across(c(AgeAccel, Born),
                list(sc = function(x) as.vector(scale(x, center = T)))))

# 1 Model epigenetic acceleration ~ birth year ====
accel_born_mod <-  brm(AgeAccel_sc ~ Born_sc,
                       data = epi_dat, family = gaussian, 
                       iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend = 'cmdstanr')

# 2 Model epigenetic acceleration ~ first repro ====
accel_fr_mod <- brm(AgeAccel_sc ~ FirstRepro_sc,
                    data = lh_epi_dat, family = gaussian, 
                    iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                    prior = prior(normal(0,1), class = b),
                    control = list(adapt_delta = 0.99, max_treedepth = 18),
                    backend = 'cmdstanr')

# 5 Model lifetime reproductive success ~ age at first reproduction ====
lrs_fr_mod <-  brm(LRS ~ FirstRepro_sc*Born_sc,
                   data = lh_pop_dat, 
                   family = negbinomial(link = 'log', link_shape = 'log'), 
                   iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                   prior = prior(normal(0,1), class = b),
                   control = list(adapt_delta = 0.99, max_treedepth = 20),
                   backend = 'cmdstanr')

# Save models
saveRDS(accel_born_mod, 'models/accel_born_mod.rds')
saveRDS(accel_fr_mod, 'models/accel_fr_mod.rds')
saveRDS(lrs_fr_mod, 'models/lrs_fr_mod.rds')
