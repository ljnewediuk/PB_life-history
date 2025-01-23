
# X5 - Supplemental quality control ====

# Author: Levi Newediuk

# Check for relationship between detection p-values and age of sample

library(tidyverse)
library(vegan)

# Load data
sample_sheet <- data.frame()
QC_checks <- data.frame()
rmv_samps <- data.frame()
t_betas <- data.frame()
for(B in c(1:3, 9)) {
  # Load sample sheets
  sheet_B <- readRDS(paste0('output/updated_sample_sheet_PB_array', B, '.rds')) %>%
    mutate(Batch = B)
  # Load transposed betas for PCA
  t_betas_B <- readRDS(paste0('output/tbetas_PB_array', B, '.rds')) %>%
    rownames_to_column('chip.ID.loc') %>%
    left_join(select(sheet_B, Sample_Name, chip.ID.loc)) %>%
    relocate(Sample_Name, .after = chip.ID.loc)
  # Load detection p-values
  QC_B <- readRDS(paste0('output/detection_p_batch', B, '.rds')) %>%
    left_join(sheet_B) %>%
    mutate(Batch = B)
  # Make list of potential samples to remove if they are more than twice the 
  # next smallest detection p value
  rmv_B <- data.frame()
  for(i in 2:nrow(QC_B)-1) {
    if(QC_B[i,][[2]]/QC_B[i+1,][[2]] > 2) {
      rmv_B <- rbind(rmv_B, QC_B[1:i,])
    }
  }
  # Add sample IDs
  rmv_B <- rmv_B %>%
    select(chip.ID.loc:Sample_Name) %>%
    mutate(Batch = B)
  # Bind with larger data frames
  sample_sheet <- rbind(sample_sheet, sheet_B)
  QC_checks <- rbind(QC_checks, QC_B)
  rmv_samps <- rbind(rmv_samps, rmv_B)
  t_betas <- rbind(t_betas, t_betas_B)
}

# Plot relationship between detection p-value and year of sample
QC_checks %>%
  mutate(YrSampled = as.numeric(substr(Sample_Name, 8, 11)),
         Label = ifelse(detection_p > 0.1, Sample_Name, '')) %>%
  ggplot(aes(x = YrSampled, y = detection_p)) +
  geom_point(colour = 'black', size = 2) +
  ylab('Detection p-value') + xlab('Sample year') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18))

# Also do so for non-outliers
QC_checks %>%
  mutate('YrSampled' = as.numeric(substr(Sample_Name, 8, 11))) %>%
ggplot() +
  geom_point(aes(x = YrSampled, y = detection_p), colour = 'black', size = 2) +
  ylab('Detection p-value') + xlab('Sample year') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18)) +
  ylim(0, 0.1)

# Save plots
ggsave('detection_p_vals.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/supplementary', dpi = 300, height = 12, width = 18, units = 'cm', bg = 'white')

# Run PCAs to visually check samples that cluster away from others
betas_ord <- t_betas %>%
  select(! Sample_Name) %>%
  column_to_rownames('chip.ID.loc') %>%
  rda(scale = T)

# Pull out PCA 1 and 2
PCA_pts <- betas_ord$CA$u[, c(1:2)] %>%
  as.data.frame() %>%
  rownames_to_column('chip.ID.loc') %>%
  # Join with sample metadata
  right_join(select(t_betas, chip.ID.loc, Sample_Name)) 

# Plot PCA
PCA_pts %>%
  mutate(Label = ifelse(PC1 > 0.05, Sample_Name, '')) %>%
ggplot(aes(x = PC1, y = PC2)) + 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(size = 2) +
  geom_text(aes(label = Label), position = position_jitter(), hjust = -.1) +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.position = c(0.8,0.6),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, colour = 'black')) +
  ylab('PC2') +
  xlab('PC1')

# To rmv_samps, add any additional samples with detection p-values > 0.1 and
# samples that cluster away from others in PCA
rmv_samp_IDs <- QC_checks %>%
  select(chip.ID.loc, detection_p, Sample_Name, Batch) %>%
  filter(detection_p > 0.1 | Sample_Name %in% c('X11418_1990-09-20_Blood', 
                                                'X09304_2001-09-06_Blood')) %>%
  rbind(rmv_samps) %>%
  distinct() %>%
  pull(Sample_Name)

# Save final list of samples to remove
saveRDS(rmv_samp_IDs, 'output/failed_QC_samples.rds')

