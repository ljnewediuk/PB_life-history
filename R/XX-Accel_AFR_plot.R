
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)
library(bayestestR)
library(performance)

# Load models
accel_fr_mod <- readRDS('models/accel_fr_mod.rds')

# Load life history/epigenetic data
lh_epi_dat <- readRDS('output/lh_info_epi.rds')  %>%
  # Scale variables
  mutate(across(c(Born:LRS, AgeAccel), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Predict model draws by sex 
# Male
accel_fr_draws_M <- lh_epi_dat |>
  data_grid(FirstRepro_sc = seq_range(FirstRepro_sc, n = 1000), 
            BearID = NA, 
            Sex = 'M') |>
  add_epred_draws(object = accel_fr_mod, ndraws = 1000)
# Female
accel_fr_draws_F <- lh_epi_dat |>
  data_grid(FirstRepro_sc = seq_range(FirstRepro_sc, n = 1000), 
            BearID = NA, 
            Sex = 'F') |>
  add_epred_draws(object = accel_fr_mod, ndraws = 1000)
# Bind together
accel_fr_draws <- accel_fr_draws_M %>%
  bind_rows(accel_fr_draws_F)

# Plot epigenetic acceleration on age at first repro
ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = '#E7E6E6') +
  scale_colour_manual(values = c('#ffd700', '#3399ff')) +
  scale_fill_manual(values = c('#ffd70070', '#3399ff70')) +
  geom_point(data = lh_epi_dat, aes(y = AgeAccel_sc, x = FirstRepro_sc, col = Sex), size = 3) +
  stat_lineribbon(data = accel_fr_draws, aes(x = FirstRepro_sc, y = .epred, fill = Sex, colour = Sex), .width = .90) + 
  scale_x_continuous(breaks = c(-2, 0, 2), labels = c(2, 10, 18)) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2), labels = c(-12, -6, 0, 6)) +
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
        axis.text = element_text(size = 18, colour = '#E7E6E6'),
        axis.title.y = element_text(size = 18, colour = '#E7E6E6', vjust = 5),
        axis.title.x = element_text(size = 18, colour = '#E7E6E6', vjust = -5)) +
  xlab('Age at first reproduction') + ylab('Epigenetic age \n acceleration (years)') +
  facet_wrap(~ Sex)

# Save
ggsave('accel_afr_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 10, width = 15, units = 'cm')
