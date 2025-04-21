
# 05 - glm Models ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Fit GLMs to test hypotheses in MS

#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)

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
  # Filter individuals born after 1996 (we might not have captured their full
  # reproductive lifespan)
  filter(Born <= 1996)

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
accel_born_mod <-  lm(AgeAccel_sc ~ Born_sc, data = epi_dat)

# 2 Model epigenetic acceleration ~ first repro ====
accel_fr_mod <- lm(AgeAccel_sc ~ FirstRepro_sc, data = lh_epi_dat)

# 3 Model lifetime reproductive success ~ age at first reproduction ====

# Negative binomial
lrs_fr_mod <- MASS::glm.nb(LRS ~ FirstRepro*Born, data = lh_pop_dat)
# Poisson
lrs_fr_mod <- glm(LRS ~ FirstRepro*Born, data = lh_pop_dat, family = poisson(link = 'log'))

# Make new data for predicted distribution
pred_dat <- data.frame(FirstRepro = rep(seq(from = min(lh_pop_dat$FirstRepro), 
                                            to = max(lh_pop_dat$FirstRepro), 
                                            length.out = 100), 3),
                      Born = c(rep(1965, times = 50), 
                               rep(1990, times = 50), 
                               rep(2000, times = 50)))
# Predict
preds <- predict(lrs_fr_mod, new_dat, type = 'response', se.fit = T)

# Add predictions and confidence intervals
pred_dat$LRS_pred <- preds$fit
pred_dat$LRS_lower <- preds$fit - preds$se.fit*2
pred_dat$LRS_upper <- preds$fit + preds$se.fit*2

# Plot
pred_dat %>%
  ggplot(aes(x = FirstRepro, y = LRS_pred, group = Born)) +
  geom_jitter(data = lh_pop_dat, aes(x = FirstRepro, y = LRS, colour = Born)) +
  geom_ribbon(aes(x = FirstRepro, ymin = LRS_lower, ymax = LRS_upper), alpha = 0.3) +
  geom_line()

# Save models
saveRDS(accel_born_mod, 'models/accel_born_mod.rds')
saveRDS(accel_fr_mod, 'models/accel_fr_mod.rds')
saveRDS(lrs_fr_mod, 'models/lrs_fr_mod.rds')
