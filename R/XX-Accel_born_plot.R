
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)
library(bayestestR)
library(performance)

# Load models
accel_born_mod <- readRDS('models/accel_born_mod.rds')

# Load epigenetic aging data
epi_dat <- readRDS('input/PB_clock_ages2.rds') %>%
  rename('Year' = yr) %>%
  mutate(Born = Year - round(Age),
         across(c(Year, AgeAccel, Age, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Predict model draws by sex 
# Male
accel_born_draws_M <- epi_dat |>
  data_grid(Born_sc = seq_range(Born_sc, n = 1000), 
            BearID = NA, 
            Sex = 'M') |>
  add_epred_draws(object = accel_born_mod, ndraws = 1000)
# Female
accel_born_draws_F <- epi_dat |>
  data_grid(Born_sc = seq_range(Born_sc, n = 1000), 
            BearID = NA, 
            Sex = 'F') |>
  add_epred_draws(object = accel_born_mod, ndraws = 1000)
# Bind together
accel_born_draws <- accel_born_draws_M %>%
  bind_rows(accel_born_draws_F)

# Plot age at first repro over time
ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = '#E7E6E6') +
  scale_colour_manual(values = c('#ffd700', '#3399ff')) +
  scale_fill_manual(values = c('#ffd70070', '#3399ff70')) +
  scale_y_continuous(breaks = c(-4.2, -2.1, -0.002, 2.2), labels = c(-12, -6, 0, 6)) +
  scale_x_continuous(breaks = c(-1.85, -0.34, 1.2), labels = c(1970, 1985, 2000)) +
  stat_lineribbon(data = accel_born_draws, aes(x = Born_sc, y = .epred, fill = Sex, colour = Sex), .width = .90) +
  geom_jitter(data = epi_dat, aes(x = Born_sc, y = AgeAccel_sc, col = Sex), size = 3, alpha = 0.5) +
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
  xlab('Year of birth') + ylab('Epigenetic age \n acceleration (years)') +
  facet_wrap(~ Sex) 

# Save
ggsave('accel_born_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 10, width = 15, units = 'cm')
