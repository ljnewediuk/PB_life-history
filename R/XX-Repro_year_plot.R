
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)
library(bayestestR)
library(performance)

# Load models
repro_yr_mod <- readRDS('models/repro_yr_mod.rds')

# Load life history data
lh_pop_dat_trunc <- readRDS('output/lh_info_pop.rds') %>%
  mutate(across(Born:LRS, 
                list(sc = function(x) as.vector(scale(x, center = T))))) %>%
  filter(Born %in% 1980:2000)

# Predict model draws by sex 
# Male
repro_yr_draws_M <- lh_pop_dat_trunc |>
  data_grid(Born_sc = seq_range(Born_sc, n = 1000), 
            BearID = NA, 
            Sex = 'M') |>
  add_epred_draws(object = repro_yr_mod, ndraws = 1000)
# Female
repro_yr_draws_F <- lh_pop_dat_trunc |>
  data_grid(Born_sc = seq_range(Born_sc, n = 1000), 
            BearID = NA, 
            Sex = 'F') |>
  add_epred_draws(object = repro_yr_mod, ndraws = 1000)
# Bind together
repro_yr_draws <- repro_yr_draws_M %>%
  bind_rows(repro_yr_draws_F)

# Plot age at first repro over time
ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = '#E7E6E6') +
  scale_colour_manual(values = c('#ffd700', '#3399ff')) +
  scale_fill_manual(values = c('#ffd70070', '#3399ff70')) +
  scale_y_continuous(breaks = c(-2.16, -0.07, 2.02), labels = c('2', '10', '18'), expand = expansion(mult = 0.3)) +
  scale_x_continuous(breaks = c(-0.45, 0.45, 1.35), labels = c(1980, 1990, 2000)) +
  stat_lineribbon(data = repro_yr_draws, aes(x = Born_sc, y = .epred, fill = Sex, colour = Sex), .width = .90) +
  geom_jitter(data = lh_pop_dat_trunc, aes(x = Born_sc, y = FirstRepro_sc, col = Sex), size = 3, alpha = 0.5) +
  theme(legend.position = 'none',
        plot.background = element_rect(fill = '#222A35', colour = NA),
        panel.background = element_rect(fill = '#222A35', 
                                        colour = '#E7E6E6',
                                        linewidth = 1),
        plot.margin = unit(c(1, 1, 1, 1), 'cm'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = '#E7E6E6'),
        axis.text.y = element_text(size = 18, colour = '#E7E6E6'),
        axis.text.x = element_text(size = 18, colour = '#E7E6E6', angle = 45, vjust = 0.45),
        axis.title.y = element_text(size = 18, colour = '#E7E6E6', vjust = 5),
        axis.title.x = element_text(size = 18, colour = '#E7E6E6', vjust = -5)) +
  xlab('Year of birth') + ylab('Age at first reproduction') +
  facet_wrap(~ Sex) 

# Save
ggsave('repro_yr_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 10, width = 15, units = 'cm')
