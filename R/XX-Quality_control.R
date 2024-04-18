
# Check for relationship between detection p-values and age of sample

library(tidyverse)

# Load chip sample sheets
sample_sheet <- readRDS('output/updated_sample_sheet_PB_array1.rds') %>%
  rbind(readRDS('output/updated_sample_sheet_PB_array2.rds')) %>%
  rbind(readRDS('output/updated_sample_sheet_PB_array3.rds')) %>%
  dplyr::rename('SampleID' = Sample_Name) %>%
  mutate('YrSampled' = as.numeric(substr(SampleID, 8, 11))) %>%
  select(chip.ID.loc, YrSampled)

# Load detection p-values
QC_checks <- readRDS('output/detection_p_batch1.rds') %>%
  rbind(readRDS('output/detection_p_batch2.rds')) %>%
  rbind(readRDS('output/detection_p_batch3.rds'))

# Join data
QC_dat <- QC_checks %>%
  left_join(sample_sheet)

# Plot relationship between detection p-value and year of sample
ggplot() +
  geom_point(data = QC_dat, aes(x = YrSampled, y = detection_p), colour = 'black', size = 2) +
  ylab('Detection p-value') + xlab('Sample year') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18))

# Also do so for non-outliers
ggplot() +
  geom_point(data = QC_dat, aes(x = YrSampled, y = detection_p), colour = 'black', size = 2) +
  ylab('Detection p-value') + xlab('Sample year') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18)) +
  ylim(0, 0.1)
