
library(tidyverse)
library(tidybayes)
library(cowplot)

# Load population and epigenetic sample life history data
lh_dat <- readRDS('output/lh_info_pop.rds')

# Load model draws
repro_yr_draws <- readRDS('models/repro_yr_draws.rds')
accel_fr_draws <- readRDS('models/accel_fr_draws.rds')

# Load life history/epigenetic data
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  group_by(BearID) %>%
  mutate(MAccel = mean(AgeAccel))

# Load epigenetic data and attach birth dates
epi_dat <- read.csv('input/bear_capture_info.csv') %>%
  select(BearCode, Born) %>%
  rename('BearID' = BearCode) %>%
  distinct() %>%
  right_join(readRDS('input/PB_clock_ages2.rds'))

# Plot age at first repro on birth date
age_born_plot <- ggplot(data = lh_dat, aes(x = Born, y = FirstRepro)) + 
  geom_point(pch = 21, size = 1, fill = '#a32372', col = '#3c0c2a') +
  stat_lineribbon(data = repro_yr_draws, 
                  aes(x = Born, y = .epred), alpha = 0.5, .width = .90,
                  colour = '#a32372', fill = '#a32372') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
        plot.margin = unit(c(0.5, 1, 1, 1), 'cm'),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5)) +
  ylab('Age at first reproduction') + xlab('Year of birth')

# Plot epigenetic acceleration on age at first repro
age_afr_plot <- ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  geom_point(data = lh_epi_dat, aes(y = MAccel, x = FirstRepro), pch = 21, fill = '#a32372', col = '#3c0c2a', size = 3) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred), alpha = 0.5, .width = .90, colour = '#a32372', fill = '#a32372') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
        plot.margin = unit(c(1, 1, 1, 1), 'cm'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(colour = 'black', face = 'bold', hjust = 0, vjust = 2, size = 22),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5)) +
  xlab('Age at first reproduction') + ylab('Epigenetic age acceleration (years)') 

# Plot  panels
panel_age_plot <- plot_grid(age_born_plot, age_afr_plot, labels = c('A', 'B'), label_size = 18, align = 'hv')

# Save
ggsave('panel_age_plot.tiff', last_plot(), path = 'figures/', device = 'tiff', dpi = 300, width = 28, height = 15, units = 'cm')

