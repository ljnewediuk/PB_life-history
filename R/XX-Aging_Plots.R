
library(tidyverse)
library(tidybayes)

# Load model draws
accel_fr_draws <- readRDS('models/accel_fr_draws.rds') %>%
  mutate(Label = factor(Sex, levels = c('F', 'M'), labels = c('A', 'B')))

# Load life history/epigenetic data
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  group_by(BearID) %>%
  mutate(MAccel = mean(AgeAccel),
         Label = factor(Sex, levels = c('F', 'M'), labels = c('A', 'B')))

lh_dat <- readRDS('output/lh_info_pop.rds')
# Load growth data
growth_dat <- readRDS('output/cohort_growth_residuals.rds') %>% 
  mutate(Cohort = ifelse(Born %in% 1965:1969, '1965-1969', '9999'),
         Cohort = ifelse(Born %in% 1970:1974, '1970-1974', Cohort),
         Cohort = ifelse(Born %in% 1975:1979, '1975-1979', Cohort),
         Cohort = ifelse(Born %in% 1980:1984, '1980-1984', Cohort),
         Cohort = ifelse(Born %in% 1985:1989, '1985-1989', Cohort),
         Cohort = ifelse(Born %in% 1990:1994, '1990-1994', Cohort),
         Cohort = ifelse(Born %in% 1995:1999, '1995-1999', Cohort),
         Cohort = ifelse(Born %in% 2000:2004, '2000-2004', Cohort),
         Cohort = ifelse(Born %in% 2005:2009, '2005-2009', Cohort),
         Cohort = ifelse(Born %in% 2005:2020, '2005-2020', Cohort)) %>%
  group_by(Age, Sex) %>%
  mutate(MSize = mean(GrowthParm),
         sdSize = sd(GrowthParm),
         ZSize = (GrowthParm - MSize)/sdSize) %>%
  ungroup() %>%
  select(BearID, Age, ResidGrowth, ZSize, GrowthParm, Cohort, Sex)  %>%
  left_join(lh_dat) %>%
  group_by(BearID, Cohort, Born, Sex) %>%
  mutate(MFR = mean(FirstRepro), MG = mean(ResidGrowth)) 

# Plot mean growth by age at first reproduction
ggplot(growth_dat, aes(x = MFR,  y = MG, group = Cohort, colour = Cohort)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = T) +
  xlab('Age at first reproduction') +  ylab('Residual growth') +
  facet_wrap(~ Sex, scales = 'free')

# Plot epigenetic acceleration on age at first repro
ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  scale_colour_manual(values = c('#ff8007', '#0028ff')) +
  scale_fill_manual(values = c('#ff800770', '#0028ff70')) +
  geom_point(data = lh_epi_dat, aes(y = MAccel, x = FirstRepro, fill = Label, col = Label), pch = 21, size = 5) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred, fill = Label, colour = Label), alpha = 0.05, .width = .90) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred, fill = Label, colour = Label), alpha = 0.2, .width = .60) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred, fill = Label, colour = Label), alpha = 0.4, .width = .30) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.margin = unit(c(1, 1, 1, 1), 'cm'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(colour = 'black', face = 'bold', hjust = 0, vjust = 2, size = 22),
        axis.line.x = element_line(linewidth = 0.6, colour = 'black'),
        axis.line.y = element_line(linewidth = 0.6, colour = 'black'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5)) +
  xlab('Age at first reproduction') + ylab('Mean epigenetic age acceleration (years)') +
  facet_wrap(~ Label, scales = 'free')

# Save
ggsave('age_accel_fig.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 18, width = 28, units = 'cm')
