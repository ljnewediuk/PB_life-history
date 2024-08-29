
# X6 - Checks for tissue bias ====

# Author: Levi Newediuk

# Check whether we can detect an age accel ~ birth year relationship with only
# single tissues

library(tidyverse)
library(brms)
library(bayesplot)
library(ggdist)

# 1 Load data ====

# Load age data
ages <- readRDS('output/PB_clock_ages.rds') %>%
  mutate(Born = floor(yr - Age)) %>%
  mutate(across(c(AgeAccel, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Filter blood samples only
ages_blood <- ages %>%
  filter(Spec == 'Blood')

# Filter blood samples only
ages_skin <- ages %>%
  filter(Spec == 'Skin')

# 2 Count number of samples ===

# How many bears have multiple samples for blood?
ages_blood %>% 
  group_by(BearID) %>%
  summarize(nSamples = n()) %>%
  group_by(nSamples) %>%
  summarize(count = n()) 

# # A tibble: 2 × 2
# nSamples count
# <int> <int>
#   1     31
#   2     4

# How many bears have multiple samples for skin?
ages_skin %>% 
  group_by(BearID) %>%
  summarize(nSamples = n()) %>%
  group_by(nSamples) %>%
  summarize(count = n()) 

# A tibble: 5 × 2
# nSamples count
# <int> <int>
#    1    44
#    2     5
#    4     1
#    9     3
#   10     1

# 3 Fit models with single tissues ====

# Model age acceleration ~ year for skin only
skin_mod <- brm(AgeAccel_sc ~ Born_sc + Sex + (Born_sc + Sex | BearID),
                 family = gaussian, data = ages_skin,
                 iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                 prior = prior(normal(0,1), class = b),
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 backend = 'cmdstanr')

# Model age acceleration ~ year for blood only
# Weak effect... but sample size is really small
blood_mod <- brm(AgeAccel_sc ~ Born_sc + Sex,
                 family = gaussian, data = ages_blood,
                 iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                 prior = prior(normal(0,1), class = b),
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 backend = 'cmdstanr')

# Probability of directional effect
bayestestR::p_direction(blood_mod)
bayestestR::p_direction(skin_mod)

# Blood....
# Parameter   |     pd
# --------------------
#    (Intercept) | 58.16%
#   Born_sc     | 59.79%
#  SexM        | 53.35%
# 
# Skin....
# Parameter   |     pd
# --------------------
#    (Intercept) | 81.72%
#   Born_sc     | 90.83%
#  SexM        | 89.36%

# 4 Check agreement between blood and skin samples ====

# Which individuals have a combination of blood and skin in the same year?
bears_by_year <- ages %>% 
  group_by(BearID, Age) %>%
  summarize(nSamples = n()) %>%
  filter(nSamples > 1) %>%
  pull(BearID)

# Plot relationship between blood and skin age estimates when same bear was
# aged at same time using both tissues
ages  %>%
  filter(BearID %in% bears_by_year) %>%
  arrange(BearID) %>%
  filter(! BearID %in% c('X12697', 'X16239', 'X19212')) %>%
  select(BearID, Spec, AgePredict) %>%
  pivot_wider(names_from = Spec, values_from = AgePredict) %>%
  ggplot(aes(x = Skin, y = Blood)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
  geom_smooth(method = 'lm', se  = F, colour = 'black') +
  geom_point(size = 2) +
  labs(x = 'Skin Age', y = 'Blood Age') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18))

# Save plot
ggsave('blood_skin_agreement.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# 5 Fit and plot age ~ birth year model with effect for tissue ====

# Fit model
tissue_mod <- brm(AgeAccel_sc ~ Born_sc + Spec + (Born_sc + Spec | BearID),
                  family = gaussian, data = ages,
                  iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                  prior = prior(normal(0,1), class = b),
                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                  backend = 'cmdstanr')

# Make new data for predicting
# Skin
nd <- expand_grid(BearID = NA,
                  Spec = c('Blood', 'Skin'),
                  Born_sc = seq(from = min(ages$Born_sc), 
                                to = max(ages$Born_sc),
                                by = 0.05))

# Extract fitted values
f_tissue <- fitted(tissue_mod,
                   newdata = nd,
                   probs = c(0.025, 0.975),
                   summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd)) %>%
  # Rename and unscale
  mutate(Born = Born_sc * sd(ages$Born) + mean(ages$Born),
         AgeAccel = value * sd(ages$AgeAccel) + mean(ages$AgeAccel)) %>%
  select(Spec, Born, AgeAccel) %>%
  # Split by tissue
  group_by(Spec) %>%
  group_split()

# Get means of posterior for lines
line_bl <- f_tissue[[1]] %>%
  group_by(Born) %>%
  summarize(AgeAccel = mean(AgeAccel))

line_sk <- f_tissue[[2]] %>%
  group_by(Born) %>%
  summarize(AgeAccel = mean(AgeAccel))

# Plot
ggplot(f_tissue[[1]], aes(x = Born, y = AgeAccel)) +
  stat_lineribbon(data = f_tissue[[1]],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = 'red') +
  stat_lineribbon(data = f_tissue[[2]],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = 'grey') +
  geom_line(data = line_bl, colour = 'red') +
  geom_line(data = line_sk, colour = '#808080') +
  scale_colour_manual(values = c('red', '#808080'), name = 'Tissue type') +
  geom_jitter(data = ages, aes(x = Born, y = AgeAccel, colour = Spec), size = 3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = c(0.8, 0.15),
        legend.text = element_text(size = 18, colour = 'black'),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.title = element_text(size = 18, colour = 'black'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Age acceleration (years)') + xlab('Year of birth')

# Save plot
ggsave('accel_born_tissue_plot.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

