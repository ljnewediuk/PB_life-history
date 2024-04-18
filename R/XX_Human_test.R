
# Test polar bear clock on human samples
# 107 human blood samples from "Methylation studies in human blood N94 N68", 
# accession # GSE184221

# Samples were run on Horvath Mammal 40 array, so should be able to age them
# using the polar bear clock

library(tidyverse)
library(sesame)
library(minfi)

# Load sample info
sample_sheet <- read.table('input/human_sample_info.txt', header = T) %>%
  mutate(chip.ID.loc = paste(Sample_Name, chip.ID, Stripe, sep = '_'),
         # Add basenames (i.e., file paths for iscan files)
         Basename = paste0('iscans/human_test/', Sample_Name, '_', chip.ID, '_', Stripe))

# Load PB clock
PB_clock <- readRDS('output/PB_clock.rds')

# Create an RGChannelSet object containing raw red green channel data from .idat
RGset <- minfi::read.metharray.exp(base = NULL, targets = sample_sheet, recursive = T) 

# Annotate the RGset object with probe coordinates
RGset@annotation <- c(array = 'HorvathMammalMethylChip40', annotation = "test.unknown")

# Get data frame of raw beta values
raw_betas_minfi <- as_tibble(minfi::getBeta(RGset), rownames = "CGid")

# Normalize betas
Mset <- minfi::preprocessNoob(RGset)

# Get data frame of normalized beta values 
n_betas <- as_tibble(minfi::getBeta(Mset), rownames = "CGid")

# 5 Transpose normalized betas
n_betas_t <- n_betas %>%
  select(! CGid) %>%
  t() %>%
  as.data.frame()

# Add back CGid by renaming rest of columns to CG sites
colnames(n_betas_t) <- n_betas$CGid

# Prep betas and select only sites in clock
meth_betas <- n_betas_t %>%
  rownames_to_column('chip.ID.loc') %>%
  select(c(chip.ID.loc, PB_clock[2:126,]$cg))

# Predict age using polar bear clock
age_preds <- meth_betas %>%
  rowwise() %>%
  # Add intercept and multiply meth betas by Cpg coefficients
  mutate(AgePredict = PB_clock[1,]$beta + 
           sum(c_across(cg01127764:cg25430089) * PB_clock[2:126,]$beta)) %>%
  # Join with sample info
  select(chip.ID.loc, AgePredict) %>%
  left_join(sample_sheet)

# Correlation
age_corr <- cor.test(y = age_preds$AgePredict, x = age_preds$Age)
corr_label <- paste0('correlation = ', round(age_corr$estimate, 2))

# MAE
age_mae <- median(abs(age_preds$AgePredict - age_preds$Age))

# Plot
ggplot() +
  geom_point(data = age_preds, aes(x = Age, y = AgePredict), colour = 'black', size = 2) +
  annotate(geom = 'text', x = 78, y = -9, label = corr_label, size = 6) +
  ylab('Epigenetic age') + xlab('Chronological age') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.key = element_rect(colour = 'white'),
        legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, colour = 'black'))

# Save plot
ggsave('human_test_plot.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 14, units = 'cm', bg = 'white')
