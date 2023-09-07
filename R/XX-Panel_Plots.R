
library(tidyverse)
library(tidybayes)
library(cowplot)

# Load model draws
accel_fr_draws <- readRDS('models/accel_fr_draws.rds')
accel_growthr_draws <- readRDS('models/accel_growthrate_draws.rds')

# Load life history/epigenetic data
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  group_by(BearID) %>%
  mutate(MAccel = mean(AgeAccel))

# Join epigenetic aging and growth rates data
epi_growth_rates_dat <- readRDS('input/PB_clock_ages2.rds') %>%
  left_join(readRDS('models/size_age_id_draws.rds')) %>%
  group_by(BearID) %>%
  mutate(MAccel = mean(AgeAccel)) %>%
  na.omit()

# Plot epigenetic acceleration on age at first repro
age_afr_plot <- ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred, colour = Sex, fill = Sex), size = 0.5, .width = .90) +
  geom_point(data = lh_epi_dat, aes(y = MAccel, x = FirstRepro, fill = Sex), pch = 21, col = 'black', size = 3) +
  scale_colour_manual(values = c('#ffd700', '#ff4040')) +
  scale_fill_manual(values = c('#ffd70050', '#ff404050')) +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
        plot.margin = unit(c(1, 1, 1, 1), 'cm'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(colour = 'black', face = 'bold', hjust = 0, vjust = 2, size = 22),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3)) +
  ylim(-13,7) +
  xlab('Age at first reproduction') + ylab('') 

# Plot age acceleration on growth rates
age_growth_plot <- ggplot(epi_growth_rates_dat, aes(y = AgeAccel, x = mean_BearSlope)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  stat_lineribbon(data = accel_growthr_draws, aes(x = mean_BearSlope, y = .epred, colour = Sex, fill = Sex), size = 0.5, .width = .90) +
  geom_point(pch = 21, size = 3, aes(fill = Sex), colour = 'black') +
  scale_colour_manual(values = c('#ffd700', '#ff4040')) +
  scale_fill_manual(values = c('#ffd70050', '#ff404050')) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 1, 1, 1), 'cm'),
        axis.text.x = element_text(size = 18, colour = 'black', vjust = 0.5),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = 'none') +
  ylim(-13,7) +
  xlab('Growth rate (cm/year)') + ylab('Epigenetic age acceleration (years)')

# Plot  panels
panel_age_plot <- plot_grid(age_growth_plot, age_afr_plot, labels = c('A', 'B'), label_size = 18, align = 'hv')

# Save
ggsave('panel_age_plot.tiff', last_plot(), path = 'figures/', device = 'tiff', dpi = 300, width = 28, height = 15, units = 'cm')
