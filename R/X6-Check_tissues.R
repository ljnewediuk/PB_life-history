
# X6 - Checks for tissue bias ====

# Author: Levi Newediuk

# Check whether we can detect an age accel ~ birth year relationship with only
# single tissues

library(tidyverse)
library(brms)
library(bayesplot)
library(ggdist)

# 1 Load data ====

# Load epigenetic age data
epi_dat <- readRDS('output/WH_combined_ages.rds') 

# 2 Organize data by sex and tissue type ====

# Summarize by individual and tissue type
tissue_dat <- epi_dat %>%
  group_by(BearID, Spec) %>%
  summarize(AgeAccel = mean(AgeAccel), Born = mean(Born)) %>%
  ungroup() %>%
  # Scale and centre variables
  mutate(across(c(AgeAccel, Born),
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Summarize by individual and keep sex as variable
sex_dat <- epi_dat %>%
  group_by(BearID) %>%
  summarize(Sex = unique(Sex), AgeAccel = mean(AgeAccel), Born = mean(Born)) %>%
  # Scale and centre variables
  mutate(across(c(AgeAccel, Born),
                list(sc = function(x) as.vector(scale(x, center = T)))))

# 3 Fit models with tissue and sex ====

# Model age acceleration ~ year for skin only
tissue_m <- brm(AgeAccel_sc ~ Born_sc*Spec,
                family = gaussian, data = tissue_dat,
                iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                prior = prior(normal(0,1), class = b),
                control = list(adapt_delta = 0.99, max_treedepth = 20),
                backend = 'cmdstanr')

# Model age acceleration ~ year for blood only
# Weak effect... but sample size is really small
sex_m <- brm(AgeAccel_sc ~ Born_sc*Sex,
                 family = gaussian, data = sex_dat,
                 iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                 prior = prior(normal(0,1), class = b),
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 backend = 'cmdstanr')

# 4 Check agreement between blood and skin samples ====

# Which individuals have a combination of blood and skin in the same year?
bears_by_year <- epi_dat %>% 
  group_by(BearID, Age) %>%
  summarize(nSamples = n()) %>%
  filter(nSamples > 1) %>%
  pull(BearID)

# Plot relationship between blood and skin age estimates when same bear was
# aged at same time using both tissues
epi_dat  %>%
  filter(BearID %in% bears_by_year) %>%
  arrange(BearID) %>%
  filter(! BearID %in% c('X12697', 'X16239', 'X19212', 'X19627')) %>%
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

# 5 Plot age ~ birth year model with effects for tissue ====

# Make new data for predicting
nd_tissue <- expand_grid(Spec = c('Blood', 'Skin'),
                  Born_sc = seq(from = min(tissue_dat$Born_sc), 
                                to = max(tissue_dat$Born_sc),
                                by = 0.05))

# Extract fitted values
f_tissue <- fitted(tissue_m,
                   newdata = nd_tissue,
                   probs = c(0.025, 0.975),
                   summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_tissue)) %>%
  # Rename and unscale
  mutate(Born = Born_sc * sd(tissue_dat$Born) + mean(tissue_dat$Born),
         AgeAccel = value * sd(tissue_dat$AgeAccel) + mean(tissue_dat$AgeAccel)) %>%
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
                  alpha = .1, size = 0, fill = '#e57368') +
  stat_lineribbon(data = f_tissue[[2]],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#8092a6') +
  geom_line(data = line_bl, colour = '#d62d20') +
  geom_line(data = line_sk, colour = '#536878') +
  scale_colour_manual(values = c('#d62d20', '#536878'), name = 'Tissue type') +
  geom_jitter(data = ages, aes(x = Born, y = AgeAccel, colour = Spec), size = 3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = 'none',
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Epigenetic age acceleration (years)') + xlab('Year of birth')

# Save plot
ggsave('accel_born_tissue_plot.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# 6 Plot age ~ birth year model with effects for sex ====

# Make new data for predicting
nd_sex <- expand_grid(Sex = c('F', 'M'),
                      Born_sc = seq(from = min(tissue_dat$Born_sc), 
                                    to = max(tissue_dat$Born_sc),
                                    by = 0.05))

# Extract fitted values
f_sex <- fitted(sex_m,
                   newdata = nd_sex,
                   probs = c(0.025, 0.975),
                   summary = F) %>%
  data.frame() %>%
  # Pivot
  pivot_longer(everything()) %>%
  bind_cols(expand_grid(draws = 1:20000, nd_sex)) %>%
  # Rename and unscale
  mutate(Born = Born_sc * sd(sex_dat$Born) + mean(sex_dat$Born),
         AgeAccel = value * sd(sex_dat$AgeAccel) + mean(sex_dat$AgeAccel)) %>%
  select(Sex, Born, AgeAccel) %>%
  # Split by tissue
  group_by(Sex) %>%
  group_split()

# Get means of posterior for lines
line_F <- f_sex[[1]] %>%
  group_by(Born) %>%
  summarize(AgeAccel = mean(AgeAccel))

line_M <- f_sex[[2]] %>%
  group_by(Born) %>%
  summarize(AgeAccel = mean(AgeAccel))

# Plot
ggplot(f_sex[[1]], aes(x = Born, y = AgeAccel)) +
  stat_lineribbon(data = f_sex[[1]],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#f78ab5') +
  stat_lineribbon(data = f_sex[[2]],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#85b5f7') +
  geom_line(data = line_F, colour = '#f15097') +
  geom_line(data = line_M, colour = '#5097f1') +
  scale_colour_manual(values = c('#f15097', '#5097f1'), name = 'Tissue type') +
  geom_jitter(data = ages, aes(x = Born, y = AgeAccel, colour = Spec), size = 3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = 'none',
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Epigenetic age acceleration (years)') + xlab('Year of birth')

# Save plot
ggsave('accel_born_sex_plot.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

