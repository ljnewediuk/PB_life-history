
# 07 - Plots ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Plots from manuscript (2-4) and posterior predictive checks in supplement (5)
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(cowplot)
library(ggdist)
library(bayesplot)
library(brms)

# 1 Load data ====

# List of fitted effects to load
fes <- list.files('output/', pattern = 'f_effects')

# Load fitted effects and assign name
for(n in 1:length(fes)) {
  fes_n <- readRDS(paste0('output/', fes[n]))
  assign(str_extract(fes[n], '[^.]+'), fes_n)
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
epi_dat <- readRDS('output/PB_clock_ages.rds') %>%
  rename('Year' = yr) %>%
  mutate(Born = Year - round(Age),
         across(c(Year, AgeAccel, Age, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# 2 Plot age acceleration on year of birth ====

# Get mean of posterior for age accel of F & M in accel ~ born models
accel_born_mean <- f_effects_accel_born %>%
group_by(Born) %>%
  summarize(AgeAccel = mean(AgeAccel))
accel_born_mean_F <- f_effects_accel_born_F %>%
  group_by(Born) %>%
  summarize(AgeAccel = mean(AgeAccel))
accel_born_mean_M <- f_effects_accel_born_M %>%
  group_by(Born) %>%
  summarize(AgeAccel = mean(AgeAccel))

# Plot age acceleration ~ birth year models (female and male models separately)
ggplot(data = f_effects_accel_born_M, aes(x = Born, y = AgeAccel)) +
  stat_lineribbon(.width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#677daf') +
  stat_lineribbon(data = f_effects_accel_born_F, 
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#86b7b9') +
  geom_line(data = accel_born_mean_M, colour = '#425d9c') +
  geom_line(data = accel_born_mean_F, colour = '#68A6A8') +
  geom_point(data = epi_dat, aes(x = Born, y = AgeAccel, colour = Sex), size = 3) +
  scale_color_manual(values = c('#68A6A8', '#425d9c')) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = c(0.8, 0.15),
        legend.text = element_text(size = 18, colour = 'black'),
        legend.title = element_text(size = 18, colour = 'black'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Year of birth') + ylab('Age acceleration (years)')

# Save plot
ggsave('accel_born_plot.tiff', plot = last_plot(), path = 'figures/main/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# 3 Plot lifetime reproductive success on year of first reproduction ====

# Get mean of posterior for all AFR in LRS ~ first repro model in 1980, 1990, 2000
lrs_fr_means <- f_effects_lrs_fr %>%
  group_by(Born, FirstRepro) %>%
  summarize(value = mean(value))

# Plot LRS ~ first repro
ggplot(f_effects_lrs_fr, aes(x = FirstRepro, y = value)) +
  geom_jitter(data = lh_pop_dat_trunc, aes(x = FirstRepro, y = LRS, colour = Born_sc), size = 3) +
  stat_lineribbon(data = f_effects_lrs_fr[f_effects_lrs_fr$Born == 2000,],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#193A82') +
  geom_line(data = lrs_fr_means[lrs_fr_means$Born == 2000,], colour = '#112757') +
  stat_lineribbon(data = f_effects_lrs_fr[f_effects_lrs_fr$Born == 1990,],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#6B84C0') +
  geom_line(data = lrs_fr_means[lrs_fr_means$Born == 1990,], colour = '#4a67ae') +
  stat_lineribbon(data = f_effects_lrs_fr[f_effects_lrs_fr$Born == 1980,],
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#BDD4FF') +
  geom_line(data = lrs_fr_means[lrs_fr_means$Born == 1980,], colour = '#8ab3ff') +
  scale_colour_gradient(low = '#BDD4FF', high = '#193A82',
                        breaks = c(-0.45, 0.45, 1.34), 
                        labels = c(1980, 1990, 2000),
                        name = 'Year born') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = c(0.8, 0.7),
        legend.text = element_text(size = 18, colour = 'black'),
        legend.title = element_text(size = 18, colour = 'black'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Age at first reproduction') + ylab('Lifetime reproductive success')

# Save plot
ggsave('lrs_fr_plot.tiff', plot = last_plot(), path = 'figures/main/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# 4 Plot age acceleration ~ age at first reproduction ====

# Get mean of posterior for age accel of F & M in accel ~ born model
accel_fr_mean_F <- f_effects_accel_fr_F %>%
  group_by(FirstRepro) %>%
  summarize(AgeAccel = mean(AgeAccel))
accel_fr_mean_M <- f_effects_accel_fr_M %>%
  group_by(FirstRepro) %>%
  summarize(AgeAccel = mean(AgeAccel))

# Plot age accel ~ first repro
ggplot(f_effects_accel_fr_M, aes(x = FirstRepro, y = AgeAccel)) +
  stat_lineribbon(data = f_effects_accel_fr_M,
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#6B84C0') +
  geom_line(data = accel_fr_mean_M, colour = '#425d9c') +
  stat_lineribbon(data = f_effects_accel_fr_F,
                  .width = seq(from = .03, to = .975, by = .03),
                  alpha = .1, size = 0, fill = '#86b7b9') +
  geom_line(data = accel_fr_mean_F, colour = '#68A6A8') +
  scale_colour_manual(values = c('#68A6A8', '#425d9c')) +
  geom_jitter(data = lh_epi_dat, aes(x = FirstRepro, y = AgeAccel, colour = Sex), size = 3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = c(0.2, 0.2),
        legend.text = element_text(size = 18, colour = 'black'),
        legend.title = element_text(size = 18, colour = 'black'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Age acceleration (years)') + xlab('Age at first reproduction')

# Save plot
ggsave('accel_fr_plot.tiff', plot = last_plot(), path = 'figures/extended/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')

# 5 Plot posterior predictive checks ====

mods <- list.files('models/', pattern = 'mod')

color_scheme_set('red')

# Load each model, plot pp_check, and assign plot name
for(n in 1:length(mods)) {
  mod <- readRDS(paste0('models/', mods[n]))
  
  y <- mod$data[,1]
  yrep <- posterior_predict(mod, ndraws = 100)
  mod_name <- str_extract(mods[n], '[^.]+')
  
  pp_plot <- ppc_dens_overlay(y, yrep) + 
    theme(axis.text.x = element_text(size = 20, family = 'sans'),
          legend.text = element_text(size = 20, family = 'sans'),
          legend.position = 'none',
          axis.line = element_line(colour = 'black', linewidth = 0.75),
          plot.tag = element_text(size = 20, family = 'sans'),
          plot.tag.position = 'top',
          plot.margin = unit(rep(0.5, times = 4), 'cm')) +
    labs(tag = mod_name)
  
  assign(paste0(mod_name, '_plot'), pp_plot)
  
}

#  Make a list of plots
all_plots <- list(accel_born_mod_plot, accel_fr_mod_plot, lrs_fr_mod_plot)

# Plots in panels
plot_grid(plotlist = all_plots, nrow = 1)

# Save plots
ggsave('pp_checks.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/extended', dpi = 300, height = 12, width = 25, units = 'cm', bg = 'white')

