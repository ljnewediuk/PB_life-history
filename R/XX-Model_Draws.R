
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)
library(bayestestR)
library(performance)

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds')
lh_epi_dat <- readRDS('output/lh_info_epi.rds')
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

# Ice-free days ~ time
ice_free_draws <- ice_free_mod %>%
  spread_draws(b_Year_sc) %>%
  median_qi(b_Year_sc, .width = c(0.9, 0.95)) %>%
  mutate(cl = c(0.9, 0.95)) %>%
  rename('slope' = b_Year_sc) %>%
  select(slope, .lower, .upper, cl) %>%
  pivot_wider(names_from = cl, values_from = c(.lower, .upper))

# Temperature ~ time
temp_draws <- temp_mod %>%
  spread_draws(b_Year_sc) %>%
  median_qi(b_Year_sc, .width = c(0.9, 0.95)) %>%
  mutate(cl = c(0.9, 0.95)) %>%
  rename('slope' = b_Year_sc) %>%
  select(slope, .lower, .upper, cl) %>%
  pivot_wider(names_from = cl, values_from = c(.lower, .upper))

# Age acceleration ~ year
accel_yr_draws <- accel_yr_mod %>%
  spread_draws(b_Year_sc)  %>%
  median_qi(b_Year_sc, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_Year_sc) %>%
  select(slope, .lower, .upper) %>%
  mutate(cl = c(0.9, 0.95)) %>%
  pivot_wider(names_from = cl, values_from = c(.lower, .upper))

# Age acceleration ~ YOB
accel_born_fixed <- accel_born_mod %>%
  spread_draws(b_Born_sc, `b_Born_sc:SexM`) %>%
  mutate(b_Born_sc_M = b_Born_sc + `b_Born_sc:SexM`)
accel_born_M <- accel_born_fixed %>%
  median_qi(b_Born_sc_M, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_Born_sc_M) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'M', cl = c(0.9, 0.95))
accel_born_draws <- accel_born_fixed %>%
  median_qi(b_Born_sc, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_Born_sc) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'F', cl = c(0.9, 0.95)) %>%
  rbind(accel_born_M) %>%
  pivot_wider(names_from = cl, values_from = c(.lower, .upper))

# Age acceleration ~ age at first reproduction
accel_fr_fixed <- accel_fr_mod %>%
  spread_draws(b_FirstRepro_sc, `b_FirstRepro_sc:SexM`) %>%
  mutate(b_FirstRepro_sc_M = b_FirstRepro_sc + `b_FirstRepro_sc:SexM`)
accel_fr_M <- accel_fr_fixed %>%
  median_qi(b_FirstRepro_sc_M, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_FirstRepro_sc_M) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'M', cl = c(0.9, 0.95))
accel_fr_draws <- accel_fr_fixed %>%
  median_qi(b_FirstRepro_sc, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_FirstRepro_sc) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'F', cl = c(0.9, 0.95)) %>%
  rbind(accel_fr_M) %>%
  pivot_wider(names_from = cl, values_from = c(.lower, .upper))

# Age at first reproduction ~ time
repro_yr_fixed <- repro_yr_mod %>%
  spread_draws(b_Born_sc, `b_Born_sc:SexM`) %>%
  mutate(b_Born_sc_M = b_Born_sc + `b_Born_sc:SexM`)
repro_yr_M <- repro_yr_fixed %>%
  median_qi(b_Born_sc_M, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_Born_sc_M) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'M', cl = c(0.9, 0.95))
repro_yr_draws <- repro_yr_fixed %>%
  median_qi(b_Born_sc, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_Born_sc) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'F', cl = c(0.9, 0.95)) %>%
  rbind(repro_yr_M) %>%
  pivot_wider(names_from = cl, values_from = c(.lower, .upper))

# Age acceleration ~ growth rate
accel_gr_fixed <- accel_gr_mod %>%
  spread_draws(b_mean_BearSlope_sc, `b_mean_BearSlope_sc:SexM`) %>%
  mutate(b_mean_BearSlope_sc_M = b_mean_BearSlope_sc + `b_mean_BearSlope_sc:SexM`)
accel_gr_M <- accel_gr_fixed %>%
  median_qi(b_mean_BearSlope_sc_M, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_mean_BearSlope_sc_M) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'M', cl = c(0.9, 0.95))
accel_gr_draws <- accel_gr_fixed %>%
  median_qi(b_mean_BearSlope_sc, .width = c(0.9, 0.95)) %>%
  rename('slope' = b_mean_BearSlope_sc) %>%
  select(slope, .lower, .upper) %>%
  mutate(sex = 'F', cl = c(0.9, 0.95)) %>%
  rbind(accel_gr_M) %>%
  pivot_wider(names_from = cl, values_from = c(.lower, .upper))

# Lifetime reproductive success ~ age acceleration
lrs_accel_draws <- lh_epi_dat |>
  data_grid(AgeAccel = seq_range(AgeAccel, n = 1000))  |>
  add_epred_draws(object = lrs_accel_mod, ndraws = 1000)

# Lifetime reproductive success ~ time
lrs_yr_draws <- lh_pop_dat_trunc |>
  data_grid(Born = seq_range(Born, n = 1000))  |>
  add_epred_draws(object = lrs_yr_mod, ndraws = 1000)

# Calculate the contrast in age acceleration between a bear born in 2010 and 
# one born in 1960 for males and females
# Get posterior predictive draws from model (scaled vals) and unscale
for(sex in c('F', 'M')) {
  
  ppd_1965 <- posterior_predict(accel_born_mod, 
                                newdata = data.frame(BearID = NA, 
                                                     Sex = sex, 
                                                     Born_sc = -2.43), 
                                re.form = NA) * sd(epi_dat$AgeAccel) + mean(epi_dat$AgeAccel) 
  ppd_2013 <- posterior_predict(accel_born_mod, 
                                newdata = data.frame(BearID = NA, 
                                                     Sex = sex, 
                                                     Born_sc = 2.58), 
                                re.form = NA) * sd(epi_dat$AgeAccel) + mean(epi_dat$AgeAccel) 
  
  # Make data frame with difference
  contrast_accel <- data.frame(years = '2013-1965', 
                               accel = (ppd_2013 - ppd_1965)[,1])
  # Get 95% percentiles for contrast and print
  cat(sex, quantile(contrast_accel$accel, probs = c(0.025, 0.5, 0.975)), '\n')
  
}

# Calculate change in temp and ice-free days 88-2016
# 1988 ice
ice_free_88 <- posterior_predict(ice_free_mod, 
                              newdata = data.frame(Year_sc = -0.9)) * sd(ice_dat$IceFree) + mean(ice_dat$IceFree)
# 2016 ice
ice_free_16 <- posterior_predict(ice_free_mod, 
                                 newdata = data.frame(Year_sc = 1.5)) * sd(ice_dat$IceFree) + mean(ice_dat$IceFree)
# Get 95% percentiles for contrast and print
quantile((ice_free_16 - ice_free_88)[,1], probs = c(0.025, 0.5, 0.975))

# 1988 temp
temp_88 <- posterior_predict(temp_mod, 
                                 newdata = data.frame(Year_sc = -0.9)) * sd(temp_dat$Temp) + mean(temp_dat$Temp)
# 2016 ice
temp_16 <- posterior_predict(temp_mod, 
                                 newdata = data.frame(Year_sc = 1.59)) * sd(temp_dat$Temp) + mean(temp_dat$Temp)
# Get 95% percentiles for contrast and print
quantile((temp_16 - temp_88)[,1], probs = c(0.025, 0.5, 0.975))

# Calculate change in LRS from age acceleration -5 to 5 years
# -5 age acceleration
young_accel <- posterior_predict(lrs_accel_mod, 
                                 newdata = data.frame(AgeAccel = 0))
# +5 age acceleration
old_accel <- posterior_predict(lrs_accel_mod, 
                                 newdata = data.frame(AgeAccel = 10))
# Get 95% percentiles for contrast and print
quantile((old_accel - young_accel)[,1], probs = c(0.025, 0.5, 0.975))

# Get posterior predictions for interaction effects in LRS ~ afr*year of birth model and LRS ~ age accel*year of birth model

# Make levels for conditional effects
condition_lvs_afr <- list(
  # Bears born in 1980, 1990, and 2000
  Born_sc = c(-0.45, 0.45, 1.34),
  # AFR ranging from 4-15
  FirstRepro_sc = seq(from = -1.75, to = 1.15, by = 0.005)
)

condition_lvs_accel <- list(
  # Bears born in 1980, 1990, and 2000
  Born_sc = c(-0.78, 0.33, 1.5),
  # AFR ranging from 4-15
  AgeAccel_sc = seq(from = -2.2, to = 2.8, by = 0.005)
)

# Conditional effects from LRS ~ afr*year of birth
lrs_afr_coneffs <- conditional_effects(lrs_fr_mod, prob = 0.95, int_conditions = condition_lvs_afr)$`FirstRepro_sc:Born_sc` %>%
  # Unscale birth dates and ages at first reproduction
  mutate(Born = round(Born_sc * sd(lh_pop_dat$Born) + mean(lh_pop_dat$Born)),
         FirstRepro = round(FirstRepro_sc * sd(lh_pop_dat$FirstRepro) + mean(lh_pop_dat$FirstRepro)))

# Conditional effects from LRS ~ age accel*year of birth
lrs_accel_coneffs <- conditional_effects(lrs_accel_yr_mod, prob = 0.95, int_conditions = condition_lvs_accel)$`AgeAccel_sc:Born_sc` %>%
  # Unscale birth dates and ages at first reproduction
  mutate(Born = round(Born_sc * sd(lh_epi_dat$Born) + mean(lh_epi_dat$Born)),
         AgeAccel = round(AgeAccel_sc * sd(lh_epi_dat$AgeAccel) + mean(lh_epi_dat$AgeAccel)))
  
ggplot(lrs_afr_coneffs, aes(x = FirstRepro_sc, y = estimate__, fill = factor(Born), colour = factor(Born))) + 
  scale_x_continuous(breaks = c(-1.38, -0.0725, 1.13), labels = c(5, 10, 15)) +
  geom_line() + 
  geom_ribbon(aes(x = FirstRepro_sc, ymin = lower__, ymax = upper__), alpha = 0.5, colour = NA)

ggplot(lrs_accel_coneffs, aes(x = AgeAccel_sc, y = estimate__, fill = factor(Born), colour = factor(Born))) + 
  # scale_x_continuous(breaks = c(-1.38, -0.0725, 1.13), labels = c(5, 10, 15)) +
  geom_line() + 
  geom_ribbon(aes(x = AgeAccel_sc, ymin = lower__, ymax = upper__), alpha = 0.5, colour = NA)

# Save posterior draws
saveRDS(temp_draws, 'models/temp_draws.rds')
saveRDS(ice_free_draws, 'models/ice_free_draws.rds')
saveRDS(accel_fr_draws, 'models/accel_fr_draws.rds')
saveRDS(repro_yr_draws, 'models/repro_yr_draws.rds')
saveRDS(accel_yr_draws, 'models/accel_yr_draws.rds')
saveRDS(accel_born_draws, 'models/accel_born_draws.rds')
saveRDS(accel_gr_draws, 'models/accel_gr_draws.rds')
saveRDS(lrs_accel_draws, 'models/lrs_accel_draws.rds')
saveRDS(lrs_accel_coneffs, 'models/lrs_accel_coneffs.rds')
saveRDS(lrs_yr_draws, 'models/lrs_yr_draws.rds')
saveRDS(lrs_afr_coneffs, 'models/lrs_afr_coneffs.rds')

