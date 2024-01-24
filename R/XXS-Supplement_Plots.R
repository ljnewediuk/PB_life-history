
library(tidyverse)
library(tidybayes)
library(brms)

# List of models to load
mods <- list.files('models/', pattern = 'mod')

# Load each model and assign name
for(n in 1:length(mods)) {
  mod <- readRDS(paste0('models/', mods[n]))
  assign(str_extract(mods[n], '[^.]+'), mod)
}

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds') %>%
  mutate(across(Born:LRS, list(sc = function(x) as.vector(scale(x, center = T)))))
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  mutate(across(c(Born:LRS, AgeAccel), list(sc = function(x) as.vector(scale(x, center = T)))))
# Limit life history data to only reliable sampling dates
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)
# Load aging data
epi_dat <- readRDS('input/PB_clock_ages2.rds') %>%
  rename('Year' = yr) %>%
  mutate(Born = Year - round(Age),
         across(c(Year, AgeAccel, Age, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Load ice breakup dates and surface temperature data
# Ice data
ice_dat <- read.csv('input/WH_ice_breakup_dates.csv') %>%
  mutate(IceFree = jday_freezeup - jday_breakup) %>%
  rename('Year' = yr) %>%
  mutate(across(c(Year, IceFree), list(sc = function(x) as.vector(scale(x, center = T))))) %>%
  filter(Year %in% c(1988:2016))
# Temperature data
temp_dat <- read.table('input/global_temp_data.txt', header = T) %>%
  rename('Temp' = No_Smoothing) %>%
  mutate(across(c(Year, Temp), list(sc = function(x) as.vector(scale(x, center = T))))) %>%
  filter(Year %in% c(1988:2016))

# Plot surface temperature ~ year
temp_dat %>%
  add_predicted_draws(temp_mod) %>%
  ggplot(aes(x = Year_sc, y = Temp_sc)) +
  stat_lineribbon(aes(y = .prediction), .width = 0.95, alpha = 0.25, fill = 'grey') +
  geom_point(colour = 'black') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Global air surface temperature °C (scaled)') +
  xlab('Year (scaled)')

# Save
ggsave('SS_temp_year.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 15, units = 'cm', bg = 'white')

# Plot ice-free days  ~ year
ice_dat %>%
  add_predicted_draws(ice_free_mod) %>%
  ggplot(aes(x = IceFree_sc, y = Temp_sc)) +
  stat_lineribbon(aes(y = .prediction), .width = 0.95, alpha = 0.25, fill = 'grey') +
  geom_point(colour = 'black') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Ice free days/year (scaled)') +
  xlab('Year (scaled)')

# Save
ggsave('SS_ice_year.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 15, units = 'cm', bg = 'white')

# Plot age acceleration ~ birth year coloured by sex
epi_dat %>%
  add_predicted_draws(accel_born_mod) %>%
  ggplot(aes(x = AgeAccel_sc, y = Born_sc)) +
  stat_lineribbon(aes(y = .prediction), .width = 0.95, alpha = 0.25, fill = 'grey') +
  geom_point(colour = 'black') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Age acceleration (scaled)') +
  xlab('Birth year (scaled)')

# Save
ggsave('SS_accel_born.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 15, units = 'cm', bg = 'white')

# Plot age acceleration ~ age at first repro coloured by sex
lh_epi_dat %>%
  group_by(Sex) %>%
  add_predicted_draws(accel_fr_mod) %>%
  ggplot(aes(x = AgeAccel_sc, y = Born_sc)) +
  stat_lineribbon(aes(y = .prediction), .width = 0.95, alpha = 0.25, fill = 'grey') +
  geom_point(colour = 'black') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Age acceleration (scaled)') +
  xlab('Birth year (scaled)')

# Save
ggsave('SS_accel_born.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 15, units = 'cm', bg = 'white')

# Plot age at first reproduction ~ birth year

# Violin plots of age at first reproduction by birth year (1980-2000)
ggplot(lh_pop_dat_trunc, aes(x = factor(Born), y = FirstRepro)) +
  geom_violin(fill = '#82d0d2', colour = 'black') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Age at first reproduction') +
  xlab('Birth year')
  
# Save
ggsave('SS_fr_violins.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 40, units = 'cm', bg = 'white')
