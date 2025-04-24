
# X8 - Estimate black bear ages ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Estimate ages for black bear samples from Czajka et al. 2024 Molecular Ecology
#Resources https://doi.org/10.1111/1755-0998.13956
#===============================================================================


#------------------------------------------------------------------------------
# load packages
library(tidyverse)

# 1 Load data and source functions ====

# Source functions to normalize the betas and age new samples
source('functions/NormalizeBetas.R')
source('functions/AgeSamples.R')

# Load PB clock
PB_clock <- readRDS('output/PB_clock.rds')

# 2 Normalize the betas ====
normBetas('_black_bears')

# 3 Estimate ages ====
bb_ages <- ageNew(batch_no = '_black_bears', clock = PB_clock, 
                  failed_s = NULL, correct_ages = F, calc_accel = F)

# Calculate accuracy ====

# MAE
median(abs(bb_ages$AgePredict - bb_ages$Age))
# Pearson's correlation
as.numeric(cor.test(bb_ages$AgePredict, bb_ages$Age)$estimate)

# 5 Plot black bear ages ====

bb_ages %>%
  ggplot(aes(x = Age, y = AgePredict)) + 
  geom_point(size = 3, color = '#420420') +
  geom_smooth(method = 'lm', se = F, color = '#420420') +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Chronological age (years)') + ylab('Epigenetic age (years)') +
  xlim(-5, 35) + ylim(-5, 35)

# Save plot
ggsave('black_bear_check.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')
