
library(tidyverse)
library(ggdist)

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

# List of models to load
mods <- list.files('models/', pattern = 'mod')

# Load each model and assign name
for(n in 1:length(mods)) {
  mod <- readRDS(paste0('models/', mods[n]))
  assign(str_extract(mods[n], '[^.]+'), mod)
}


# Fitted effects for age accel ~ year born mod

# New data for model
# Female
nd_accel_born_F <- expand_grid(BearID = NA,
                               Sex = 'F',
                               Born_sc = seq(from = min(epi_dat$Born_sc), 
                                             to = max(epi_dat$Born_sc),
                                             by = 0.05))
# Male
nd_accel_born_M <- nd_accel_born_F %>%
  mutate(Sex = 'M')

# Fitted effects
# Female
f_accel_born_F <- fitted(accel_born_mod,
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
f_accel_born_M <- fitted(accel_born_mod,
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

# Fitted effects for age at first repro ~ year model

# New data for model
nd_repro_year <- expand_grid(Sex = c('M', 'F'),
                             BearID = NA,
                             Born_sc = seq(from = min(lh_pop_dat_trunc$Born_sc), 
                                           to = max(lh_pop_dat_trunc$Born_sc),
                                           by = 0.05))

# Fitted effects
f_repro_year <- fitted(repro_yr_mod,
                       newdata = nd_repro_year,
                       probs = c(0.025, 0.975),
                       summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_repro_year)) %>%
  # Rename and unscale
  mutate(Born = Born_sc * sd(lh_pop_dat$Born) + mean(lh_pop_dat$Born),
         FirstRepro = value * sd(lh_pop_dat$FirstRepro) + mean(lh_pop_dat$FirstRepro)) %>%
  select(Born, FirstRepro)

# Plot
f_repro_year %>%
  ggplot(aes(x = Born, y = FirstRepro)) +
  stat_lineribbon(.width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = 'grey') +
  geom_jitter(data = lh_pop_dat_trunc, aes(x = Born, y = FirstRepro))

# Fitted effects for lifetime repro success ~ age at first repro

# New data for model
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

# Plot
ggplot(f_lrs_fr, aes(x = FirstRepro, y = value)) +
  stat_lineribbon(data = f_lrs_fr[f_lrs_fr$Born == 2000,],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#193A82') +
  stat_lineribbon(data = f_lrs_fr[f_lrs_fr$Born == 1990,],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#6B84C0') +
  stat_lineribbon(data = f_lrs_fr[f_lrs_fr$Born == 1980,],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#BDD4FF') +
  scale_colour_gradient(low = '#BDD4FF', high = '#193A82',
                        breaks = c(-0.45, 0.45, 1.34), 
                        labels = c(1980, 1990, 2000)) +
  geom_jitter(data = lh_pop_dat_trunc, aes(x = FirstRepro, y = LRS, colour = Born_sc))

# Fitted effects for age acceleration ~ age at first reproduction

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

# Plot
ggplot(f_accel_fr_M, aes(x = FirstRepro, y = AgeAccel)) +
  stat_lineribbon(data = f_accel_fr_M,
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#6B84C0') +
  stat_lineribbon(data = f_accel_fr_F,
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = 'green') +
  geom_jitter(data = lh_epi_dat, aes(x = FirstRepro, y = AgeAccel, colour = Sex))

# Save fitted effects
saveRDS(f_accel_born_F, 'output/f_effects_accel_born_F.rds')
saveRDS(f_accel_born_M, 'output/f_effects_accel_born_M.rds')
saveRDS(f_accel_fr_M, 'output/f_effects_accel_fr_M.rds')
saveRDS(f_accel_fr_F, 'output/f_effects_accel_fr_F.rds')
saveRDS(f_lrs_fr, 'output/f_effects_lrs_fr.rds')
saveRDS(f_repro_year, 'output/f_effects_repro_year.rds')


