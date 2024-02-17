
library(tidyverse)
library(ggdist)

# 07 - Extract fitted values from models

# 1 Load data ====

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T)))))
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  mutate(across(c(Born:LRS, AgeAccel), list(sc = function(x) as.vector(scale(x, center = T)))))

# Limit life history data to only reliable sampling dates
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)

# Load aging data
epi_dat <- readRDS('output/PB_clock_ages.rds') %>%
  # Add necessary columns
  rename('Year' = yr) %>%
  mutate(Born = Year - round(Age),
         across(c(Year, AgeAccel, Age, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Make aging data for males and females
epi_dat_M <- epi_dat %>%
  filter(Sex == 'M')
epi_dat_F <- epi_dat %>%
  filter(Sex == 'F')

# List of models to load
mods <- list.files('models/', pattern = 'mod')

# Load each model and assign name
for(n in 1:length(mods)) {
  mod <- readRDS(paste0('models/', mods[n]))
  assign(str_extract(mods[n], '[^.]+'), mod)
}

# 2 Fitted effects for age acceleration ~ year of birth (M & F separate) ====

# Make new data for predicting
# Female model
nd_accel_born_F <- expand_grid(BearID = NA,
                               Born_sc = seq(from = min(epi_dat_F$Born_sc), 
                                             to = max(epi_dat_F$Born_sc),
                                             by = 0.05))
# Male model
nd_accel_born_M <- expand_grid(BearID = NA,
                               Born_sc = seq(from = min(epi_dat_M$Born_sc), 
                                             to = max(epi_dat_M$Born_sc),
                                             by = 0.05))

# Extract fitted values
# Female
f_accel_born_F <- fitted(accel_born_mod_F,
            newdata = nd_accel_born_F,
            probs = c(0.025, 0.975),
            summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_accel_born_F)) %>%
  # Rename and unscale
  mutate(Born = Born_sc * sd(epi_dat$Born) + mean(epi_dat$Born),
         AgeAccel = value * sd(epi_dat$AgeAccel) + mean(epi_dat$AgeAccel)) %>%
  select(Born, AgeAccel)
# Male
f_accel_born_M <- fitted(accel_born_mod_M,
                         newdata = nd_accel_born_M,
                         probs = c(0.025, 0.975),
                         summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_accel_born_M)) %>%
  # Rename and unscale
  mutate(Born = Born_sc * sd(epi_dat$Born) + mean(epi_dat$Born),
         AgeAccel = value * sd(epi_dat$AgeAccel) + mean(epi_dat$AgeAccel)) %>%
  select(Born, AgeAccel)


# 3 Fitted effects for LRS ~ age at first reproduction model by decade ====

# New data for model (-0.45, 0.45, and 1.35 correspond to 1980, 1990, 2000)
nd_lrs_fr <- expand_grid(Born_sc = c(-0.45, 0.45, 1.35),
                         BearID = NA,
                         Sex = NA,
                         FirstRepro_sc = seq(from = min(lh_pop_dat_trunc$FirstRepro_sc), 
                                             to = max(lh_pop_dat_trunc$FirstRepro_sc),
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
nd_accel_fr_F <- expand_grid(Sex = 'F',
                   BearID = NA,
                   FirstRepro_sc = seq(from = min(lh_epi_dat$FirstRepro_sc), 
                                       to = max(lh_epi_dat$FirstRepro_sc),
                                       by = 0.05))
nd_accel_fr_M <- nd_accel_fr_F %>%
  mutate(Sex = 'M')

# Fitted effects
# Female
f_accel_fr_F <- fitted(accel_fr_mod,
                   newdata = nd_accel_fr_F,
                   probs = c(0.025, 0.975),
                   summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_accel_fr_F)) %>%
  # Rename and unscale
  mutate(AgeAccel = value * sd(lh_epi_dat$AgeAccel) + mean(lh_epi_dat$AgeAccel),
         FirstRepro = FirstRepro_sc * sd(lh_epi_dat$FirstRepro) + mean(lh_epi_dat$FirstRepro))
# Male
f_accel_fr_M <- fitted(accel_fr_mod,
                       newdata = nd_accel_fr_M,
                       probs = c(0.025, 0.975),
                       summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_accel_fr_M)) %>%
  # Rename and unscale
  mutate(AgeAccel = value * sd(lh_epi_dat$AgeAccel) + mean(lh_epi_dat$AgeAccel),
         FirstRepro = FirstRepro_sc * sd(lh_epi_dat$FirstRepro) + mean(lh_epi_dat$FirstRepro))

# 5 Save fitted effects ====

saveRDS(f_accel_born_F, 'output/f_effects_accel_born_F.rds')
saveRDS(f_accel_born_M, 'output/f_effects_accel_born_M.rds')
saveRDS(f_accel_fr_M, 'output/f_effects_accel_fr_M.rds')
saveRDS(f_accel_fr_F, 'output/f_effects_accel_fr_F.rds')
saveRDS(f_lrs_fr, 'output/f_effects_lrs_fr.rds')
saveRDS(f_repro_year, 'output/f_effects_repro_year.rds')


