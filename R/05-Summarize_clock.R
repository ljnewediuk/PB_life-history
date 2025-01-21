

# 05 - Summarize clock results ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Calculate clock performance metrics (median absolute error and correlation),
#make table summarizing test samples, and plot the clock
#===============================================================================


#------------------------------------------------------------------------------
# load packages
library(tidyverse)
library(cowplot)

# 1 Load data and source functions ====

# Load sample specs
sample_specs <- data.frame()
for(B in c(1:3, 9)) {
  S <- readRDS(paste0('input/batch', B, '_samples.rds')) %>% mutate(Batch = B)
  sample_specs <- bind_rows(sample_specs, S)
}

# Load ages
epi_ages <- readRDS('output/WH_combined_ages.rds') %>%
  # Factor tissue types so skin appears first in plots
  mutate(Spec = factor(Spec, levels = c('Skin', 'Blood')))

# Load IDs of samples used to make clock versus testing
PB_clock_IDs <- readRDS('output/PB_clock_IDs.rds')

# Load failed QC samples
failedQC <- readRDS('output/failed_QC_samples.rds')

# 2 Make table of training/testing bears for supplement ====

supp_table <- sample_specs %>%
  # Distinguish sample by training, validation, or additional samples from after
  # the clock was constructed
  mutate(SampleType = case_when(sampleId %in% PB_clock_IDs$train ~ 'Training',
                                sampleId %in% PB_clock_IDs$test ~ 'Validation',
                                Batch == 9 ~ 'Additional')) %>%
  # Omit NAs and samples that failed QC
  na.omit() %>%
  filter(! sampleId %in% failedQC) %>%
  # Fix typo in tissue type
  mutate(Spec = str_replace(Spec, '_Bloo', 'Blood')) %>%
  # Clean up columns
  rename('ID' = id, 'SampleTissue' =  Spec, 
         'DateSampled' = YMD, 'Age' = age, 'Sex' = sex) %>%
  select(ID, DateSampled, SampleType, SampleTissue, Age, Sex)

# 3 Calculate clock median absolute error and correlation ====

# MAE
median(abs(epi_ages$AgePredict - epi_ages$Age))
# Pearson's correlation
as.numeric(cor.test(epi_ages$AgePredict, epi_ages$Age)$estimate)

# 5 Plot clock ====

# Function to plot clock panel
plotClock <- function(dat, panelP = F) {
  
  P <- ggplot(dat) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    geom_smooth(aes(x = Age, y = AgePredict), 
                colour = 'black', linewidth = 1.5, method = 'lm', se = F) +
    geom_point(aes(x = Age, y = AgePredict, colour = Spec), size = 5) +
    scale_colour_manual(values = c('#536878', '#d62d20')) +
    ylab('Epigenetic age (years)') + xlab('Chronological age (years)') +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), 'cm'),
          panel.background = element_rect(fill = 'white', colour = 'black',
                                          linewidth = 1),
          panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),
          axis.title = element_blank(),
          axis.text = element_text(colour = 'black', size = 18),
          legend.position = 'none')
  
  if(isTRUE(panelP)) {
    P <- P + 
      lims(x = c(0, 31), y = c(0, 31)) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 2), 'cm'))
  }
  
  return(P)
}

# Pull out repeat bears (at least four samples over lifetime) and plot panels
epi_repeats <- epi_ages %>%
  group_by(BearID) %>%
  mutate(N = n()) %>%
  filter(N >= 4)

# Plot panel for each repeat bear and list
ID_plots <- lapply(group_split(group_by(epi_repeats, BearID)), 
                    function(x) plotClock(x, panelP = T))
ID_panels <- plot_grid(plotlist = ID_plots, labels = LETTERS[2:6], 
          label_size = 22, ncol = 1, align = 'v', label_x = -.005)

# Plot all samples in testing data
clock_plot <- plotClock(epi_ages)
clock_panel <- plot_grid(clock_plot, labels = 'A', label_size = 22, ncol = 1)

# Plot panels
plot_grid(clock_panel, ID_panels, ncol = 2, rel_widths = c(1, 0.5))

# Save plot
ggsave('clock_panel.tiff', plot = last_plot(), path = 'figures/main/', 
       device = 'tiff', dpi = 300, height = 22, width = 33, units = 'cm', bg = 'white')

# Save table for supplement
write.csv(supp_table, 'output/supplementary_bear_data.csv', row.names = F)
