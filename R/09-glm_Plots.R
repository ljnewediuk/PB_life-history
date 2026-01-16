
# 07 - Plots ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Plots from manuscript (2-4) and posterior predictive checks in supplement (5)
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(cowplot)

# 1 Load data ====

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  filter(! is.infinite(FirstRepro) | ! is.infinite(LastRepro)) %>%
  # Filter individuals born after 1996 (we might not have captured their full
  # reproductive lifespan)
  filter(Born <= 1996) %>%
  # Scale and centre variables
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T)))))

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

# Load models in list
mod_list <- lapply(list.files('models', pattern = 'lm', full.names = T), readRDS)

# 2 Plot age acceleration on year of birth ====

ggplot(data = epi_dat, aes(x = Born, y = AgeAccel)) +
  stat_smooth(method = 'lm', linewidth = 0.5, fill = '#425d9c50', color = '#425d9c') +
  geom_point(colour = '#425d9c', size = 3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Year of birth') + ylab('Age acceleration (years)')

# Save plot
ggsave('accel_born_plot.tiff', plot = last_plot(), path = 'figures/main/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# Save as svg
ggsave('accel_born_plot.svg', plot = last_plot(), path = 'figures/presentation/', 
       device = 'svg', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# 3 Plot lifetime reproductive success on year of first reproduction ====

# Make new data for predicted distributions (to show relationship for different
# cohorts)
pred_dat <- data.frame(FirstRepro = rep(seq(from = min(lh_pop_dat$FirstRepro), 
                                            to = max(lh_pop_dat$FirstRepro), 
                                            length.out = 100), 3),
                       Born = c(rep(1965, times = 50), 
                                rep(1980, times = 50), 
                                rep(1995, times = 50)))
# Predict 
preds <- predict(mod_list[[3]], pred_dat, type = 'response', se.fit = T)

# Add predictions and confidence intervals
pred_dat$LRS_pred <- preds$fit
pred_dat$LRS_lower <- preds$fit - preds$se.fit*2
pred_dat$LRS_upper <- preds$fit + preds$se.fit*2

# Plot
pred_dat %>%
  ggplot(aes(x = FirstRepro, y = LRS_pred, group = Born)) +
  scale_colour_gradient(high = '#BDD4FF', low = '#193A82',
                        breaks = c(1965, 1980, 1995),
                        name = 'Year born') +
  geom_jitter(data = lh_pop_dat, aes(x = FirstRepro, y = LRS, colour = Born)) +
  geom_ribbon(data = pred_dat[pred_dat$Born == '1965',],
              aes(x = FirstRepro, ymin = LRS_lower, ymax = LRS_upper),
              alpha = 0.5, fill = '#193A82') +
  geom_ribbon(data = pred_dat[pred_dat$Born == '1980',],
              aes(x = FirstRepro, ymin = LRS_lower, ymax = LRS_upper),
              alpha = 0.5, fill = '#6B84C0') +
  geom_ribbon(data = pred_dat[pred_dat$Born == '1995',],
              aes(x = FirstRepro, ymin = LRS_lower, ymax = LRS_upper),
              alpha = 0.5, fill = '#BDD4FF') +
  geom_line(data = pred_dat[pred_dat$Born == '1995',], color = '#94A9CC') +
  geom_line(data = pred_dat[pred_dat$Born == '1980',], color = '#5268A0') +
  geom_line(data = pred_dat[pred_dat$Born == '1965',], color = '#102A5C')  +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        legend.position = c(0.8, 0.7),
        legend.text = element_text(size = 15, colour = 'black'),
        legend.title = element_text(size = 15, colour = 'black', vjust = 5),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Age at first reproduction') + ylab('Lifetime reproductive success')

# Save plot
ggsave('lrs_fr_plot.tiff', plot = last_plot(), path = 'figures/main/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# 4 Plot age acceleration ~ age at first reproduction ====

ggplot(data = lh_epi_dat, aes(x = FirstRepro, y = AgeAccel)) +
  stat_smooth(method = 'lm', linewidth = 0.5, fill = '#425d9c50', color = '#425d9c') +
  geom_point(colour = '#425d9c', size = 3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Age at first reproduction (years)') + ylab('Age acceleration (years)')

# Save plot
ggsave('accel_fr_plot.tiff', plot = last_plot(), path = 'figures/main/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

