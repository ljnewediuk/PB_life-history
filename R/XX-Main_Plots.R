
library(tidyverse)
library(ggdist)
library(cowplot)

# *** TO DO:
# 2) Get draws for conditional effects from accel ~ year model for orchard plot

# Load life-history data
lh_pop_dat <- readRDS('output/lh_info_pop.rds')
lh_epi_dat <- readRDS('output/lh_info_epi.rds')
# Limit life history data to only reliable sampling dates
lh_pop_dat_trunc <- lh_pop_dat %>%
  filter(Born %in% 1980:2000)

# List of draws to load
# Orchard plot
orchard_draws <- list.files('models/', pattern = 'draws') %>%
  setdiff(c('lrs_yr_draws.rds', 'lrs_accel_draws.rds',
            'size_age_id_draws_M.rds', 'size_age_id_draws_F.rds'))
# Panel plot
panel_draws <- c('lrs_yr_draws.rds', 'lrs_accel_draws.rds', 'lrs_afr_coneffs.rds')

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

# Factor models for plotting
orchard_plots <- orchard_combined %>%
  # Remove growth models
  filter(! model %in% c('accel_gr_draws', 'accel_yr_draws',
                        'ice_free_draws', 'temp_draws')) %>%
  mutate(model = ifelse(! is.na(sex), paste0(model, '_', sex), model)) %>%
  # Factor to order models
  mutate(model = factor(model, 
                        levels = c('accel_fr_draws_F', 'accel_fr_draws_M', 
                                   'repro_yr_draws_F', 'repro_yr_draws_M',
                                   'accel_born_draws_F', 'accel_born_draws_M')))

# Define list of labels
orchard_labs <- c('', 'C', '', 'B', '', 'A')
# Define colours
orchard_colours <- c('#fadf00', '#001BFA',  '#fadf00', '#001BFA', 
                     '#fadf00', '#001BFA')

# Orchard plot
# This is good but slopes are all on different scales... Will need to scale
# variables in models
ggplot() + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  theme(axis.title.y = element_blank()) +
  geom_vline(xintercept = 2.5, linetype = 'solid') +
  geom_vline(xintercept = 4.5, linetype = 'solid') +
  annotate(geom = 'text', x = 5.5, y = 1.2, label = 'Epigenetic aging accelerated \n over time in both sexes', size = 4.5, hjust = 0) +
  annotate(geom = 'text', x = 3.5, y = 1.2, label = 'Bears born later reproduced \n earlier in both sexes', size = 4.5, hjust = 0) +
  annotate(geom = 'text', x = 1.5, y = 1.2, label = 'Bears that first reproduced \n when older aged slower \n in both sexes', size = 4.5, hjust = 0) +
  geom_linerange(data = orchard_plots, 
                 aes(x = model, ymin = .lower_0.95, ymax = .upper_0.95, colour = model),
                 lwd = 2.5) +
  geom_linerange(data = orchard_plots, 
                 aes(x = model, ymin = .lower_0.9, ymax = .upper_0.9, colour = model),
                 lwd = 3.5) +
  geom_pointrange(data = orchard_plots,
                  aes(x = model, y = slope, ymin = .lower_0.9, ymax = .upper_0.9, colour = model),
                  lwd = 1, shape = 21, fill = 'white', stroke = 3) +
  coord_flip(ylim = c(-1,1), clip = 'off') + ylab('slope') +
  scale_colour_manual(values = orchard_colours) + 
  scale_x_discrete(labels = orchard_labs, expand = c(0,0.3)) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = .5),
        plot.margin = unit(c(0.5, 7, 1, 0.5), 'cm'),
        axis.text.y = element_text(size = 18, colour = 'black', face = 'bold'), 
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, vjust = -3),
        axis.ticks.y = element_blank(),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),
        legend.position = 'none') +
  ylab('Slope estimate') 

# Save
ggsave('orchard_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 14, width = 18, units = 'cm')

# Panel plot for lifetime reproductive success
# LRS ~ age acceleration
lrs_age_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_jitter(data = lh_epi_dat, aes(y = LRS, x = AgeAccel), pch = 20, fill = '#cbbeb5', alpha = 0.1, size = 2) +
  stat_lineribbon(data = lrs_accel_draws, aes(x = AgeAccel, y = .epred), linewidth = 0.75, colour = '#82d0d2', fill = '#82d0d270', .width = .90) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0, 0.75, 0.5), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),) +
  ylim(0, 16) +
  xlab('Age acceleration (years)')
# LRS ~ time
lrs_time_plot <- ggplot() +
  geom_jitter(data = lh_pop_dat_trunc, aes(y = LRS, x = Born), pch = 20, fill = '#cbbeb5', alpha = 0.1, size = 2) +
  stat_lineribbon(data = lrs_yr_draws, aes(x = Born, y = .epred), linewidth = 0.75, colour = '#82d0d2', fill = '#82d0d270', .width = .90) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25), 
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.75, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylim(0, 16) +
  xlab('Birth year') + ylab('Lifetime reproductive \n success')
# LRS ~ first reproduction conditional effects
lrs_coneff_plot <- ggplot(lrs_afr_coneffs, aes(x = FirstRepro_sc, y = estimate__)) +
  scale_x_continuous(breaks = c(-1.38, -0.0725, 1.13), labels = c(5, 10, 15)) +
  geom_line(aes(colour = factor(Born))) + 
  geom_ribbon(aes(x = FirstRepro_sc, ymin = lower__, ymax = upper__, fill = factor(Born)), alpha = 0.5, colour = NA) +
  scale_colour_manual(name = 'Birth year', values = c('#c0e7e8', '#82d0d2', '#68A6A8')) +
  scale_fill_manual(name = 'Birth year', values = c('#c0e7e870', '#82d0d270', '#68A6A870')) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = c(0.8, 0.7),
        legend.text = element_text(size = 18, colour = 'black'),
        legend.title = element_text(size = 18, colour = 'black'),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(0.25, 0, 0.75, 0.5), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Age at first reproduction')

# Make y axis label
Ylab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Lifetime reproductive success', size = 7,angle = 90) + theme_void()

# Plot panels
panel_lrs_plot <- plot_grid(lrs_coneff_plot, NULL, lrs_age_plot,
                            labels = c('A', '', 'B'), rel_widths = c(1, 0.3, 1), label_size = 18, label_x = -.15, nrow = 1, align = 'hv')

# Plot big figure
plot_grid(Ylab, NULL, panel_lrs_plot, nrow = 1, rel_widths = c(0.1, 0.1, 1.2))

# Save
ggsave('lrs_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 12, width = 25, units = 'cm', bg = 'white')

