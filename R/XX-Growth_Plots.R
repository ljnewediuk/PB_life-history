
library(tidyverse)

# Load predicted growth data
size_dat <- readRDS('output/cohort_growth_data.rds') %>%
  filter(Age %in% c(5, 10)) %>%
  filter(Sex == 'F' & Age < 10 | Sex == 'M' & Age > 5)

# Load estimated growth rates
growth_rates <- readRDS('models/size_age_id_draws.rds')

# Load epigenetic age data
age_dat <- readRDS('input/PB_clock_ages2.rds')

# Load growth model draws
accel_growthr_draws <- readRDS('models/accel_growthrate_draws.rds')

# Plot predicted size of mature bears by year 
ggplot(size_dat, aes(x = Cohort, y = PredGrowth, group = Sex)) +
  geom_line(aes(colour = Sex), size = 1) + 
  geom_point(aes(fill = Sex), pch = 21, size = 3) +
  scale_colour_manual(values = c('#e72131', '#3a4c97')) +
  scale_fill_manual(values = c('#e72131', '#3a4c97')) +
  ylab('Body length at maturity (cm)') + xlab('Cohort') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 1, 1, 1), 'cm'),
        axis.text.x = element_text(size = 18, colour = 'black', angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = 'none')

# Save
ggsave('size_cohort.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 19, units = 'cm')

# Join epigenetic aging and growth rates
epi_growth_rates_dat <- age_dat %>%
  left_join(growth_rates) %>%
  group_by(BearID) %>%
  mutate(MAccel = mean(AgeAccel))

# Plot age acceleration on growth rates
ggplot(epi_growth_rates_dat, aes(x = MAccel, y = mean_BearSlope)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  geom_point(pch = 21, size = 3, fill = '#a32372', colour = '#3c0c2a') +
  stat_lineribbon(data = accel_growthr_draws, aes(x = mean_BearSlope, y = .epred), colour = '#a32372', fill = '#a32372', alpha = 0.5, .width = .90) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 1, 1, 1), 'cm'),
        axis.text.x = element_text(size = 18, colour = 'black', vjust = 0.5),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = 'none') +
  xlab('Growth rate (cm/year)') + ylab('Epigenetic age acceleration (years)')


# Save
ggsave('age_growthr_fig.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 17, units = 'cm')
