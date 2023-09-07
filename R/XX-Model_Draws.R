
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
lh_epi_dat_trunc <- lh_epi_dat %>%
  filter(Born %in% 1980:2000)
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)

# List of models to load
mods <- list.files('models/', pattern = 'mod')

# Load each model and assign name
for(n in 1:length(mods)) {
  mod <- readRDS(paste0('models/', mods[n]))
  assign(str_extract(mods[n], '[^.]+'), mod)
}

# Get predicted draws from posteriors for plotting
# Acceleration ~ age at first repro
# accel_fr_draws <- lh_epi_dat |>
#   data_grid(FirstRepro = seq_range(FirstRepro, n = 1000), Sex = c('F', 'M'), BearID = NA)  |>
#   add_epred_draws(accel_fr_mod, ndraws = 1000)
# # Acceleration ~ year
# accel_year_draws <- epi_dat |>
#   data_grid(yr = seq_range(yr, n = 1000), Sex = c('M', 'F'), BearID = NA)  |>
#   add_epred_draws(accel_year_mod, ndraws = 1000)
# # Acceleration ~ growth
# accel_growthr_draws <- growth_rates |>
#   data_grid(mean_BearSlope = seq_range(mean_BearSlope, n = 1000), Sex = c('F', 'M'), BearID = NA)  |>
#   add_epred_draws(accel_growthr_mod, ndraws = 1000)

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
lrs_accel_draws <- lh_epi_dat_trunc |>
  data_grid(AgeAccel = seq_range(AgeAccel, n = 1000))  |>
  add_epred_draws(object = lrs_accel_mod, ndraws = 1000)

# Lifetime reproductive success ~ time
lrs_yr_draws <- lh_pop_dat_trunc |>
  data_grid(Born = seq_range(Born, n = 1000))  |>
  add_epred_draws(object = lrs_yr_mod, ndraws = 1000)

# Calculate the contrast in age acceleration between a bear born in 2023 and 
# one born in 1988
# Get posterior predictive draws from model for 1988 and 2023
ppd_88 <- posterior_predict(accel_yr_mod, newdata = data.frame(BearID = NA, Sex = NA, yr = 1988))
ppd_23 <- posterior_predict(accel_yr_mod, newdata = data.frame(BearID = NA, Sex = NA, yr = 2023))
# Make data frame with differnce
contrast_accel <- data.frame(years = '2023-1988', 
                             accel = (ppd_23 - ppd_88)[,1])
# Get 95% percentiles for contrast
quantile(contrast_accel$accel, probs = c(0.025, 0.5, 0.975))

# Get posterior predictions for interaction effect in LRS ~ afr*year of birth model
# ************* Need to figure out how to get the expected effect of the interaction in
# different years
# 
# This shows conditional effects of the interaction
foo <- conditional_effects(lrs_fr_mod)$`FirstRepro_sc:Born_sc`



c_eff_lh <- function(mod, sc_bdates, afr, p) {
  
  # Conditional effects from model with specified probability
  ce <- conditional_effects(lrs_fr_mod, prob = p, conditions = sc_bdates)$`FirstRepro_sc:Born_sc` %>%
    # Scale the conditional effects between -1 and 1
    mutate(across(c(estimate__:upper__), 
                  list(sc = function(x) as.vector(scales::rescale(x, to = c(0, 1)))))) %>%
    # Unscale birth dates and ages at first reproduction
    mutate(Born = round(Born_sc * sd(lh_pop_dat$Born) + mean(lh_pop_dat$Born)),
           FirstRepro = round(FirstRepro_sc * sd(lh_pop_dat$FirstRepro) + mean(lh_pop_dat$FirstRepro))) %>%
    # Filter out specified ages at first reproduction
    filter(FirstRepro %in% afr) %>%
    # Summarize effects by birth date and age at first repro
    group_by(Born, FirstRepro) %>%
    summarize(LRS = mean(estimate___sc), lower = mean(upper___sc), upper = mean(lower___sc)) %>%
    # Add column for probability
    mutate(probs = p)
  
  return(ce)
  
}

lrs_afr_coneffs <- c_eff_lh(mod = lrs_fr_mod, sc_bdates = c(-0.45, 0.45, 1.34), afr = c(4, 15), p = 0.95)

# # New AFR data sequence to predict
# fr_seq <- seq(from = min(lh_pop_dat_trunc$FirstRepro_sc), 
#               to = max(lh_pop_dat_trunc$FirstRepro_sc), 
#               length.out = 300)
# 
# # Define value of Born_sc variable in 1980, 1990, and 2000
# b_80 <- lh_pop_dat_trunc %>%
#   filter(Born == 1980) %>%
#   head(1) %>%
#   pull(Born_sc)
# b_90 <- lh_pop_dat_trunc %>%
#   filter(Born == 1990) %>%
#   head(1) %>%
#   pull(Born_sc)
# b_00 <- lh_pop_dat_trunc %>%
#   filter(Born == 2000) %>%
#   head(1) %>%
#   pull(Born_sc)
# # Create new data for predicting LRS ~ AFR for different years of birth
# # yob = 1980
# new_data_80 <- data.frame(FirstRepro_sc = fr_seq, Born_sc = b_80)
# # yob = 1990
# new_data_90 <- data.frame(FirstRepro_sc = fr_seq, Born_sc = b_90)
# # yob = 2000
# new_data_00 <- data.frame(FirstRepro_sc = fr_seq, Born_sc = b_00)
# 
# pp_data <- data.frame(FirstRepro = fr_seq, 
#                       LRS80 = colMeans(posterior_predict(lrs_fr_mod, new_data_80)),
#                       LRS90 = colMeans(posterior_predict(lrs_fr_mod, new_data_90)),
#                       LRS00 = colMeans(posterior_predict(lrs_fr_mod, new_data_00)))
# # Plot lines for interaction
# ggplot() +
#   geom_jitter(data = lh_pop_dat_trunc, aes(x = FirstRepro_sc, y = LRS), alpha = 0.3, colour = 'black') +
#   geom_line(data = pp_data, aes(x = FirstRepro, y = LRS80), colour = '#e9b91c') +
#   geom_line(data = pp_data, aes(x = FirstRepro, y = LRS90), colour = '#ce7b12') +
#   geom_line(data = pp_data, aes(x = FirstRepro, y = LRS00), colour = '#ae1324') +
#   theme(panel.background = element_rect(colour = 'black', fill = 'white'))

# Save posterior draws
saveRDS(temp_draws, 'models/temp_draws.rds')
saveRDS(ice_free_draws, 'models/ice_free_draws.rds')
saveRDS(accel_fr_draws, 'models/accel_fr_draws.rds')
saveRDS(repro_yr_draws, 'models/repro_yr_draws.rds')
saveRDS(accel_yr_draws, 'models/accel_yr_draws.rds')
saveRDS(accel_gr_draws, 'models/accel_gr_draws.rds')
saveRDS(lrs_accel_draws, 'models/lrs_accel_draws.rds')
saveRDS(lrs_yr_draws, 'models/lrs_yr_draws.rds')
saveRDS(lrs_afr_coneffs, 'models/lrs_afr_coneffs.rds')

