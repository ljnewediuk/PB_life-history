
library(tidyverse)
library(ggdist)
library(cowplot)

# *** TO DO:
# 2) Get draws for conditional effects from accel ~ year model for orchard plot

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds')
lh_epi_dat <- readRDS('output/lh_info_epi.rds')
# Limit life history data to only reliable sampling dates
lh_epi_dat_trunc <- lh_epi_dat %>%
  filter(Born %in% 1980:2000)
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)

# List of draws to load
# Orchard plot
orchard_draws <- list.files('models/', pattern = 'draws') %>%
  setdiff(c('lrs_yr_draws.rds', 'lrs_accel_draws.rds',
            'size_age_id_draws_M.rds', 'size_age_id_draws_F.rds'))
# Panel plot
panel_draws <- c('lrs_yr_draws.rds', 'lrs_accel_draws.rds')

# Load panel model and assign name
for(n in 1:length(panel_draws)) {
  mod <- readRDS(paste0('models/', panel_draws[n]))
  assign(str_extract(panel_draws[n], '[^.]+'), mod)
}

# Load draws for orchard plot and combine into table
orchard_combined <- data.frame()
for(n in 1:length(orchard_draws)) {
  draw <- readRDS(paste0('models/', orchard_draws[n])) %>%
    mutate(model = str_extract(orchard_draws[n], '[^.]+')) 
  orchard_combined <- plyr::rbind.fill(orchard_combined, draw)
}

# Load conditional effects for lrs ~ afr
lrs_afr_draws <- readRDS('models/lrs_afr_coneffs.rds') %>%
  ungroup() %>%
  rename('slope' = LRS, '.lower_0.95' = lower, '.upper_0.95' = upper) %>%
  mutate(sex = NA, .lower_0.9 = NA, .upper_0.9 = NA, 
         model = paste('born', Born, 'afr', FirstRepro, sep = '_')) %>%
  select(slope, sex, .lower_0.9, .upper_0.9, .lower_0.95, .upper_0.95, model)

# Factor models for plotting
orchard_plots <- orchard_combined %>%
  # Remove growth models
  filter(! model == 'accel_gr_draws') %>%
  mutate(model = ifelse(! is.na(sex), paste0(model, '_', sex), model)) %>%
  # Add lrs ~ afr conditional effects
  rbind(lrs_afr_draws) %>%
  mutate(model = factor(model, 
                        levels = c('born_1995_afr_15', 'born_1995_afr_4',
                                   'born_1989_afr_15', 'born_1989_afr_4',
                                   'born_1983_afr_15', 'born_1983_afr_4',
                                   'accel_fr_draws_F', 'accel_fr_draws_M', 
                                   'repro_yr_draws_F', 'repro_yr_draws_M',
                                   'accel_yr_draws', 'ice_free_draws',
                                   'temp_draws')))

# Define list of labels
orchard_labs <- c('', '', '', '', '', 'E', '', 'D', '', 'C', 'B', '', 'A')
# Define colours
orchard_colours <- c('#f25022', '#ffb900', '#f25022', '#ffb900', '#f25022', 
                     '#ffb900', '#ff9af1', '#7fbfff',  '#ff9af1', '#7fbfff', 
                     '#cbbeb5', '#82d0d2', '#d80000')

# Orchard plot
# This is good but slopes are all on different scales... Will need to scale
# variables in models
ggplot() + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  theme(axis.title.y = element_blank()) +
  geom_vline(xintercept = 2.5, linetype = 'dashed') +
  geom_vline(xintercept = 4.5, linetype = 'dashed') +
  geom_vline(xintercept = 6.5, linetype = 'solid') +
  geom_vline(xintercept = 8.5, linetype = 'solid') +
  geom_vline(xintercept = 10.5, linetype = 'solid') +
  geom_vline(xintercept = 11.5, linetype = 'solid') +
  annotate(geom = 'text', x = 12.5, y = -0.4, label = 'climate ~ time', size = 5) +
  annotate(geom = 'text', x = 11, y = -0.45, label = 'age accel. ~ time', size = 5) +
  annotate(geom = 'text', x = 9.5, y = 0.35, label = 'afr ~ time', size = 5) +
  annotate(geom = 'text', x = 7.5, y = 0.45, label = 'age accel. ~ afr', size = 5) +
  annotate(geom = 'text', x = 5.5, y = -0.5, label = 'lifetime \n reproductive \n success', size = 5) +
  annotate(geom = 'text', x = 5.5, y = 0.5, label = '1985', size = 5) +
  annotate(geom = 'text', x = 3.5, y = 0.5, label = '1989', size = 5) +
  annotate(geom = 'text', x = 1.5, y = 0.5, label = '1995', size = 5) +
  geom_linerange(data = orchard_plots, 
                 aes(x = model, ymin = .lower_0.95, ymax = .upper_0.95, colour = model),
                 lwd = 1.5) +
  geom_linerange(data = orchard_plots, 
                 aes(x = model, ymin = .lower_0.9, ymax = .upper_0.9, colour = model),
                 lwd = 2.5) +
  geom_pointrange(data = orchard_plots,
                  aes(x = model, y = slope, ymin = .lower_0.9, ymax = .upper_0.9, colour = model),
                  lwd = 1, shape = 21, fill = 'white', stroke = 2) +
  coord_flip() + ylab('slope') +
  scale_colour_manual(values = orchard_colours) + 
  scale_x_discrete(labels = orchard_labs, expand = c(0,0.3)) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        plot.margin = unit(c(0.5, 0.5, 1, 0.5), 'cm'),
        axis.text.y = element_text(size = 18, colour = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, vjust = -3),
        axis.ticks.y = element_blank(),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),
        legend.position = 'none') +
  ylab('Slope estimate') 

# Save
ggsave('orchard_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 20, width = 12, units = 'cm')

# Panel plot for lifetime reproductive success
# LRS ~ age acceleration
lrs_age_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_jitter(data = lh_epi_dat_trunc, aes(y = LRS, x = AgeAccel), pch = 20, fill = '#cbbeb5', alpha = 0.1, size = 2) +
  stat_lineribbon(data = lrs_accel_draws, aes(x = AgeAccel, y = .epred), linewidth = 0.75, colour = '#82d0d2', fill = '#82d0d270', .width = .90) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),) +
  ylim(0, 16) +
  xlab('Age acceleration (years)')
# LRS ~ time
lrs_time_plot <- ggplot() +
  geom_jitter(data = lh_pop_dat_trunc, aes(y = LRS, x = Born), pch = 20, fill = '#cbbeb5', alpha = 0.1, size = 2) +
  stat_lineribbon(data = lrs_yr_draws, aes(x = Born, y = .epred), linewidth = 0.75, colour = '#82d0d2', fill = '#82d0d270', .width = .90) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        axis.text.y = element_text(size = 18, colour = 'black'), 
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 4),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.75, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylim(0, 16) +
  xlab('Year of birth') + ylab('Lifetime reproductive \n success')

# Plot panels
panel_lrs_plot <- plot_grid(lrs_time_plot, NULL, lrs_age_plot,
                            labels = c('A', '', 'B'), rel_widths = c(1, 0, 1), label_size = 18, align = 'hv', label_x = .1, nrow = 1)

# Save
ggsave('lrs_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 10, width = 23, units = 'cm')

