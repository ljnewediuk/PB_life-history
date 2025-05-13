
library(tidyverse)
library(EnvStats)
library(cowplot)

# Load capture data and select bear IDs and aging technique
capt_dat <- read.csv('input/bear_capture_info.csv') %>%
  rename(BearID = BearCode) %>%
  select(BearID, AQual)

# Load aging data
epi_dat <- readRDS('output/WH_combined_ages.rds') %>%
  # Summarize by individual
  group_by(BearID) %>%
  summarize(AgeAccel = mean(AgeAccel), Born = mean(Born)) %>%
  # Scale and centre variables
  mutate(across(c(AgeAccel, Born),
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Run Rosner's test to check for outliers, get row numbers of outliers
outlier_rows <- rosnerTest(epi_dat$AgeAccel)$all.stats$Obs.Num

# Exclude suspected outliers from the dataset and re-run the model
sub_epi_dat <- epi_dat[-outlier_rows ,]

# Original model
accel_born_mod <-  lm(AgeAccel_sc ~ Born_sc, data = epi_dat)

# Without suspected outliers 
accel_born_mod_no <-  lm(AgeAccel_sc ~ Born_sc, data = sub_epi_dat)

# Subset only the outliers and join aging technique
epi_outliers <- epi_dat[outlier_rows ,] %>%
  left_join(capt_dat) %>%
  select(BearID, AQual) %>%
  distinct()

# Plot the relationship with the outliers highlighted, and excluded in a different panel

# Panel A - outliers highlighted
p_outliers <- ggplot(data = epi_dat, aes(x = Born, y = AgeAccel)) +
  stat_smooth(method = 'lm', linewidth = 0.5, fill = '#425d9c50', color = '#425d9c') +
  geom_point(colour = '#425d9c', size = 3) +
  geom_point(data = epi_outliers, colour = 'red') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylim(-18, 12) +
  xlab('Year of birth') + ylab('Age acceleration (years)') 

# Panel B - outliers removed
p_no_outliers <- ggplot(data = sub_epi_dat, aes(x = Born, y = AgeAccel)) +
  stat_smooth(method = 'lm', linewidth = 0.5, fill = '#425d9c50', color = '#425d9c') +
  geom_point(colour = '#425d9c', size = 3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylim(-18, 12)

# Make panel plot

# Plot panels
out_panel <- plot_grid(p_outliers, p_no_outliers, 
                         ncol = 2, labels = c('A', 'B'), label_size = 22)

# X and Y labels for ewas panel plot
Ylab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Epigenetic age acceleration (years)', 
                             size = 7, angle = 90) + theme_void()
Xlab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Year of birth', 
                             size = 7, hjust = 0.4) + theme_void()

# Add axis labels
out_panel_y <- plot_grid(Ylab, out_panel, rel_widths = c(0.1, 1))
plot_grid(out_panel_y, Xlab, rel_heights = c(1, 0.08), ncol = 1)

# Save plot
ggsave('outliers_plot.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 27, units = 'cm', bg = 'white')
