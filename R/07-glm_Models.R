
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

# Scaled
accel_born_mod <-  lm(AgeAccel_sc ~ Born_sc, data = epi_dat)
# Unscaled for prediction
accel_born_mod_unsc <-  lm(AgeAccel ~ Born, data = epi_dat)

# Check how much faster a bear aged over its lifetime from 1960-2020
# Predicted age acceleration in 1965
pred_65 <- predict(accel_born_mod_unsc, data.frame(Born = 1965), se.fit = T)
# Predicted age acceleration in 2020
pred_20 <- predict(accel_born_mod_unsc, data.frame(Born = 2020), se.fit = T)
# Difference in mean age acceleration
delta_y <- pred_20$fit - pred_65$fit
# Standard error of difference
delta_se <- sqrt(pred_20$se.fit^2 + pred_65$se.fit^2)
# 95% confidence interval
delta_ci <- c(delta_y - 1.96 * delta_se, delta_y + 1.96 * delta_se)

# 2 Model epigenetic acceleration ~ first repro ====

accel_fr_mod <- lm(AgeAccel_sc ~ FirstRepro_sc, data = lh_epi_dat)

# 3 Model lifetime reproductive success ~ age at first reproduction ====

# Without scaled variables for plotting
lrs_fr_mod <- MASS::glm.nb(LRS ~ FirstRepro*Born, data = lh_pop_dat)
# Scaled
lrs_fr_mod_sc <- MASS::glm.nb(LRS ~ FirstRepro_sc*Born_sc, data = lh_pop_dat)

# 4 Save models ====
saveRDS(accel_born_mod, 'models/accel_born_mod_lm.rds')
saveRDS(accel_fr_mod, 'models/accel_fr_mod_lm.rds')
saveRDS(lrs_fr_mod, 'models/lrs_fr_mod_nb_glm.rds')
