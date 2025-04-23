
# X6 - Checks for tissue bias ====

# Author: Levi Newediuk

# Check whether we can detect an age accel ~ birth year relationship with only
# single tissues

library(tidyverse)
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
tissue_m <- lm(AgeAccel_sc ~ Born_sc*Spec, data = tissue_dat)

# Model age acceleration ~ year for blood only
# Weak effect... but sample size is really small
sex_m <- lm(AgeAccel_sc ~ Born_sc*Sex, data = sex_dat)

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

ggplot(tissue_dat, aes(x = Born, y = AgeAccel)) +
  stat_smooth(method = 'lm', linewidth = 0.5, aes(fill = Spec, color = Spec)) +
  geom_point(aes(colour = Spec), size = 3) +
  scale_fill_manual(values = c('#d62d20', '#536878')) +
  scale_color_manual(values = c('#d62d20', '#536878')) +
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

ggplot(sex_dat, aes(x = Born, y = AgeAccel)) +
  stat_smooth(method = 'lm', linewidth = 0.5, aes(fill = Sex, color = Sex)) +
  geom_point(aes(colour = Sex), size = 3) +
  scale_fill_manual(values = c('#f15097', '#5097f1')) +
  scale_color_manual(values = c('#f15097', '#5097f1')) +
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

