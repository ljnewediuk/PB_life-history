
library(bayesplot)
library(tidyverse)
library(brms)
library(cowplot)

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

#  List plots
all_plots <- list(accel_born_mod_plot, accel_fr_mod_plot, lrs_fr_mod_plot)
  
# Plot grid
plot_grid(plotlist = all_plots, nrow = 1)

# Save
ggsave('pp_checks.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 12, width = 25, units = 'cm', bg = 'white')


