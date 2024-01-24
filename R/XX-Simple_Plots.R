
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)
library(ggdist)
library(cowplot)

# Model epigenetic acceleration ~ first repro
accel_fr_mod <- brm(AgeAccel ~ FirstRepro + (1 | BearID),
                    data = lh_epi_dat, family = gaussian, 
                    iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                    prior = prior(normal(0,1), class = b),
                    control = list(adapt_delta = 0.99, max_treedepth = 18),
                    backend = 'cmdstanr')

# Get model draws for plotting
accel_fr_draws <- lh_epi_dat |>
  data_grid(FirstRepro = seq_range(FirstRepro, n = 1000), BearID = NA)  |>
  add_epred_draws(object = accel_fr_mod, ndraws = 1000)

# Plot regression
ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro, y = .epred), linewidth = 0.75, colour = 'black', fill = 'lightgrey', .width = .90, alpha = 0.5) +
  geom_jitter(data = lh_epi_dat, aes(y = AgeAccel, x = FirstRepro), size = 2) +
  scale_y_continuous(breaks = c(-10, 0, 10)) +
  scale_x_continuous(breaks = c(0, 10, 20))  +
  theme(panel.background = element_rect(colour = 'black', 
                                        fill = 'white', 
                                        linewidth = 1),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 14, 
                                   colour = 'black'),
        axis.title.x = element_text(size = 16, 
                                    colour = 'black', 
                                    vjust = -5),
        axis.title.y = element_text(size = 16, 
                                    colour = 'black', 
                                    vjust = 5)) +
  ylab('Age acceleration (years)') + 
  xlab('Age at first reproduction') +
  ylim(-10, 10) +
  xlim(0, 20)

# Save
ggsave('simple_accel_afr.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 10, width = 13, units = 'cm', bg = 'white')

