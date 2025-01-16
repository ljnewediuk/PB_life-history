
# 04 - Fitted effects ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Extract fitted values from Bayesian GLMs
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(ggdist)

# 1 Load data ====

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  filter(! is.infinite(FirstRepro) | ! is.infinite(LastRepro)) %>%
  # Scale and centre variables
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T))))) %>%
  # Filter individuals born after 2000 (we might not have captured their full
  # reproductive lifespan)
  filter(Born <= 2000 & Born >= 1980)

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

# 2 Load models ====

# List of models to load
mods <- list.files('models/', pattern = 'mod')

# Load each model and assign name
for(n in 1:length(mods)) {
  mod <- readRDS(paste0('models/', mods[n]))
  assign(str_extract(mods[n], '[^.]+'), mod)
}

# 2 Fitted effects for age acceleration ~ year of birth ====

# Make new data for predicting
nd_accel_born <- expand_grid(
  Born_sc = seq(from = min(epi_dat$Born_sc), 
                to = max(epi_dat$Born_sc),
                by = 0.05))

# Extract fitted values
f_accel_born <- fitted(accel_born_mod,
            newdata = nd_accel_born,
            probs = c(0.025, 0.975),
            summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_accel_born)) %>%
  # Rename and unscale
  mutate(Born = Born_sc * sd(epi_dat$Born) + mean(epi_dat$Born),
         AgeAccel = value * sd(epi_dat$AgeAccel) + mean(epi_dat$AgeAccel)) %>%
  select(Born, AgeAccel)

# 3 Fitted effects for LRS ~ age at first reproduction model by decade ====

# New data for model (-0.45, 0.45, and 1.3 correspond to 1980, 1990, 2000)
nd_lrs_fr <- expand_grid(Born_sc = c(-1.3, 0.3, 1.9),
                         BearID = NA,
                         Sex = NA,
                         FirstRepro_sc = seq(
                           from = min(lh_pop_dat$FirstRepro_sc), 
                           to = max(lh_pop_dat$FirstRepro_sc),
                           by = 0.05))
# Fitted effects
f_lrs_fr <- fitted(lrs_fr_mod,
                   newdata = nd_lrs_fr,
                   probs = c(0.025, 0.975),
                   summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_lrs_fr)) %>%
  # Rename and unscale
  mutate(Born = floor(Born_sc * sd(lh_pop_dat$Born) + mean(lh_pop_dat$Born)),
         FirstRepro = FirstRepro_sc * sd(lh_pop_dat$FirstRepro) + mean(lh_pop_dat$FirstRepro))

# 4 Fitted effects for age acceleration ~ age at first reproduction ====

# New data for model
nd_accel_fr <- expand_grid(
  FirstRepro_sc = seq(
    from = min(lh_epi_dat$FirstRepro_sc), 
    to = max(lh_epi_dat$FirstRepro_sc),
    by = 0.05))

# Fitted effects
f_accel_fr <- fitted(accel_fr_mod,
                   newdata = nd_accel_fr,
                   probs = c(0.025, 0.975),
                   summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_accel_fr)) %>%
  # Rename and unscale
  mutate(AgeAccel = value * sd(lh_epi_dat$AgeAccel) + mean(lh_epi_dat$AgeAccel),
         FirstRepro = FirstRepro_sc * sd(lh_epi_dat$FirstRepro) + mean(lh_epi_dat$FirstRepro))

# 5 Save fitted effects ====

saveRDS(f_accel_born, 'output/f_effects_accel_born.rds')
saveRDS(f_accel_fr, 'output/f_effects_accel_fr.rds')
saveRDS(f_lrs_fr, 'output/f_effects_lrs_fr.rds')


