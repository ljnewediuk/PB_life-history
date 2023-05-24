
library(tidyverse)
library(tidybayes)

# Load model draws
accel_fr_draws <- readRDS('models/accel_fr_draws.rds') %>%
  mutate(Label = factor(Sex, levels = c('F', 'M'), labels = c('A', 'B')))
repro_yr_draws <- readRDS('models/repro_yr_draws.rds')
accel_yr_draws <- readRDS('models/accel_yr_draws.rds')

# Load population and epigenetic sample life history data
lh_dat <- readRDS('output/lh_info_pop.rds')
lh_epi_dat <- readRDS('output/lh_info_epi.rds') %>%
  group_by(BearID) %>%
  mutate(MAccel = mean(AgeAccel),
         Label = factor(Sex, levels = c('F', 'M'), labels = c('A', 'B')))

# Load epigenetic data and attach birth dates
epi_dat <- read.csv('input/bear_capture_info.csv') %>%
  select(BearCode, Born) %>%
  rename('BearID' = BearCode) %>%
  distinct() %>%
  right_join(readRDS('input/PB_clock_ages2.rds'))

# Load ice break-up dates
ice_dat <- read.csv('input/WH_ice_breakup_dates.csv')

# Load temperature data
# (temp is global surface temperature compared to long-run average 1951-1980)
temp_dat <- read.table('input/global_temp_data.txt', header = T) %>%
  filter(Year >= min(lh_dat$Born) & Year <= max(lh_dat$Born))

# Calculate cumulative ice-free days over the average N years to first reproduction
afr <- 5
# OR age of senescence 
afr <- 20

cum_ice_free <- data.frame()
for(i in min(ice_dat$yr) : (max(ice_dat$yr) - afr)) {
  ice_free_d <- ice_dat %>%
    filter(yr %in% i:(i + afr)) %>%
    mutate(ice_free = jday_freezeup - jday_breakup) %>%
    pull(ice_free) %>%
    sum() %>% 
    as.numeric()
  ice_row <- data.frame(Born = i,
                        IceFreeDays = ice_free_d)
  cum_ice_free <- rbind(cum_ice_free, ice_row)
}

# Add cumulative ice-free days between birth and year of sample to epigenetic data
samples_ice <- data.frame()
for(i in 1: nrow(epi_dat)) {
  
  ice_free_d <- ice_dat %>%
    filter(yr %in% epi_dat[i,]$Born: epi_dat[i,]$yr) %>%
    mutate(ice_free = jday_freezeup - jday_breakup) %>%
    pull(ice_free) %>%
    sum() %>%
    as.numeric()
  ice_row <- data.frame(SampleID = epi_dat[i,]$SampleID,
                        IceFreeDays = ice_free_d,
                        IceFreeDaysYr = ice_free_d/ceiling(epi_dat[i,]$Age))
  samples_ice <- rbind(samples_ice, ice_row)
}


# Scale temperature and ice-free days
cum_ice_free_sc <- cum_ice_free %>%
  mutate(ScaledIceFree = scales::rescale(IceFreeDays)*3)
temp_dat_sc <- temp_dat %>%
  mutate(ScaledTemp = scales::rescale(No_Smoothing)*5)

# Plot age at first repro on birth date
ggplot(data = lh_dat, aes(x = Born, y = FirstRepro)) + 
  geom_point(aes(col = Sex, fill = Sex), alpha = 0.5, pch = 21, size = 2) + 
  stat_lineribbon(data = repro_yr_draws, aes(x = Born, y = .epred, fill = Sex, colour = Sex), alpha = 0.05, .width = .90) +
  stat_lineribbon(data = repro_yr_draws, aes(x = Born, y = .epred, fill = Sex, colour = Sex), alpha = 0.2, .width = .60) +
  stat_lineribbon(data = repro_yr_draws, aes(x = Born, y = .epred, fill = Sex, colour = Sex), alpha = 0.4, .width = .30) +
  scale_colour_manual(values = c('#2e21a1', '#2154a1')) +
  scale_fill_manual(values = c('#2e21a1', '#2154a1')) +
  geom_line(data = cum_ice_free_sc, aes(x = Born, y = ScaledIceFree), colour = '#209999') + 
  geom_area(data = cum_ice_free_sc, aes(x = Born, y = ScaledIceFree), 
            fill = '#afeeee', alpha = 0.5) + 
  geom_line(data = temp_dat_sc, aes(x = Year, y = ScaledTemp), col = '#cc181e') +
  geom_segment(aes(x = 1978, xend = 1982, y = 0, yend = 0), colour = '#209999') +
  annotate('text', x = 1975, y = 0, label = '784 d', size = 5, colour = '#209999') +
  geom_segment(aes(x = 1988, xend = 1992, y = 1.2, yend = 1.2), colour = '#209999') +
  annotate('text', x = 1985, y = 1.2, label = '873 d', size = 5, colour = '#209999') +
  geom_segment(aes(x = 2007, xend = 2011, y = 3, yend = 3), colour = '#209999') +
  annotate('text', x = 2003.7, y = 3, label = '1008 d', size = 5, colour = '#209999') +
  scale_y_continuous(sec.axis = sec_axis(~ . * 0.16, breaks = c(0, 0.5, 1), 
                                         name = 'Surface temp. (°C above long-term average)')) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.margin = unit(c(0.5, 1, 1, 1), 'cm'),
        axis.line.x = element_line(linewidth = 0.6, colour = 'black'),
        axis.line.y.left = element_line(linewidth = 0.6, colour = 'black'),
        axis.line.y.right = element_line(linewidth = 0.6, colour = 'black'),
        axis.text.y.left = element_text(size = 18, colour = 'black'),
        axis.text.y.right = element_text(size = 18, colour = '#cc181e'),
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.ticks.y.right = element_line(colour = '#cc181e'),
        axis.title.y.left = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.y.right = element_text(size = 18, colour = '#cc181e', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5)) +
  ylab('Age at first reproduction') + xlab('Year of birth')

# Save
ggsave('lhs_year_fig.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 18, width = 22, units = 'cm')

# Plot epigenetic acceleration on age at first repro
ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  scale_colour_manual(values = c('#2e21a1', '#2154a1')) +
  scale_fill_manual(values = c('#2e21a1', '#2154a1')) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred, fill = Label, colour = Label), alpha = 0.05, .width = .90) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred, fill = Label, colour = Label), alpha = 0.2, .width = .60) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred, fill = Label, colour = Label), alpha = 0.4, .width = .30) +
  geom_point(data = lh_epi_dat, aes(y = MAccel, x = FirstRepro, fill = Label), pch = 21, size = 5) +
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

# Add ice-free days to epigenetic data
accel_dat <- epi_dat %>%
  # left_join(cum_ice_free) %>%
  left_join(ice_dat) %>%
  mutate(IceFreeSeasonLength = jday_freezeup - jday_breakup,
         IceFreeSeasonLength_sc = scales::rescale(IceFreeSeasonLength)*5 - 10)

# Get repeat bears
rep_bears <- accel_dat %>%
  group_by(BearID) %>%
  mutate(N_samples = n()) %>%
  filter(N_samples > 2) %>%
  arrange(IceFreeDays) %>%
  mutate(Age30 = Born + 30,
         BearID_f = factor(BearID, 
                           levels = c('X12224', 'X12008', 'X12697', 'X12606', 'X19212'),
                           labels = c('1991', '1993', '1995', '1996', '1998')))

# Plot epigenetic aging on year
ggplot(data = accel_dat, aes(x = yr, y = AgeAccel)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black') +
  stat_lineribbon(data = accel_yr_draws, aes(x = yr, y = .epred), colour = 'black', fill = 'black', alpha = 0.05, .width = .90) +
  stat_lineribbon(data = accel_yr_draws, aes(x = yr, y = .epred), colour = 'black', fill = 'black', alpha = 0.2, .width = .60) +
  stat_lineribbon(data = accel_yr_draws, aes(x = yr, y = .epred), colour = 'black', fill = 'black', alpha = 0.4, .width = .30) +
  # geom_smooth(method = 'lm', colour = 'black') +
  geom_point(pch = 21, size = 5, fill = 'grey') +
  geom_line(aes(x = yr, y = IceFreeSeasonLength_sc), colour = '#209999') +
  scale_y_continuous(sec.axis = sec_axis(~ .,
                                         breaks = c(-10, -7.5, -5),
                                         labels = c(130, 160, 190),
                                         name = 'Length of ice-free season (days)')) +
  theme(panel.background = element_rect(fill = 'white', colour = NA),
        plot.margin = unit(c(0.5, 1, 1, 1), 'cm'),
        axis.line.x = element_line(linewidth = 0.6, colour = 'black'),
        axis.line.y = element_line(linewidth = 0.6, colour = 'black'),
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.text.y.left = element_text(size = 18, colour = 'black'),
        axis.text.y.right = element_text(size = 18, colour = '#209999'),
        axis.title.y.left = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.y.right = element_text(size = 18, colour = '#209999', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5),
        axis.ticks.y.right = element_line(colour = '#209999'),
        legend.key = element_rect(fill = 'white', colour = 'white'),
        legend.text = element_text(size = 15, colour = 'black'),
        legend.title = element_text(size = 18, colour = 'black'),
        legend.position = c(.3, .2)) +
  ylab('Epigenetic age acceleration (years)') + xlab('Year of sample')

# Save
ggsave('age_year_fig.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 18, width = 22, units = 'cm')




