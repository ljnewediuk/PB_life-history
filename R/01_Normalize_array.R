
# 01 - Normalize array ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Normalize methylation data
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(sesame)
library(minfi)

# Get normalized betas using Minfi, transpose betas, then save 
# the transposed betas and sample sheet with updated columns for aging. 
# Batch number needs to be specified, corresponding to both the sample sheet 
# for the batch and the folder containing the idat files for the batch. This 
# code is adapted from some version of code from the Clock Foundation and 
# A. Shafer.

# Specify batch number
batch_no <- 1

# 1 Install packages if not already installed ====

# SeSaMe and minfi packages from Bioconductor:
bioc_packs <- c('minfi', 'sesame')
# Install if one or either is not
if(length(setdiff(bioc_packs, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(bioc_packs, rownames(installed.packages())))
}

# Then install sesame and minfi developmental versions:
dev_packs <- c('HorvathMammalMethylChip40manifest', 'HorvathMammalMethylChip40anno.test.unknown')
# Install if one or either is not
if(length(setdiff(dev_packs, rownames(installed.packages()))) > 0) {
  install.packages(paste0(dev_packs, '_0.2.2.tar.gz'), type = 'source', repos = NULL)
}

# 2 Load data ====

# Load sample sheet
sample_sheet_file_name <- paste0('input/PB_array', batch_no, '_sample_sheet', batch_no, '.rds')

# Load manifest file
manifest_file_name_sesame <- 'input/HorvathMammal40.CanonicalManifest.3.2019.sesame.csv'

# 3 Update sample sheet ====

# List idat file names to remove any chip positions without data
chip.IDs <- substr(list.files(paste0('iscans/batch', batch_no, '/'), recursive = T, full.names = F), 1, 12)
chip.positions <- substr(list.files(paste0('iscans/batch', batch_no, '/'), recursive = T, full.names = F), 27, 32)

# Remove any files without both red and green (i.e., two files)
run_samples <- data.frame(chip.ID.loc = paste(chip.IDs, chip.positions, sep = '_')) %>%
  group_by(chip.ID.loc) %>%
  summarize(n_samples = n()) %>%
  filter(! n_samples == 1)

# Sample sheet needs a column 'Basename' pointing to the basename of a two-colour
# .idat file (i.e., either _Red.idat or _Grn.idat). Need to add the col Basename
# with the file path for each sample in the corresponding row. Then use 
# 'read.metharray.exp'to find the corresponding files using the sample sheet.
sample_sheet <- readRDS(sample_sheet_file_name) %>%
  mutate(chip.ID.loc = paste(chip.ID, stripe, sep = '_'),
    # Add basenames (i.e., file paths for iscan files)
    Basename = paste0('iscans/batch', 
                      batch_no, '/', chip.ID, '/', chip.ID, '_', stripe)) %>%
  # Filter only iscans with both red and green iscan files
  filter(chip.ID.loc %in% run_samples$chip.ID.loc)

# 4 Normalize betas ====

# Create an RGChannelSet object containing raw red green channel data from .idat
RGset <- minfi::read.metharray.exp(base = NULL, targets = sample_sheet, recursive = T) 

# Annotate the RGset object with probe coordinates. 
# This line is currently just a place holder as the annotation is empty, 
# but needs to be run anyway.
RGset@annotation <- c(array = 'HorvathMammalMethylChip40', annotation = "test.unknown")

# Calling getBeta on the RGset object will return a data frame of raw beta values 
# for each CG site.
raw_betas_minfi <- as_tibble(minfi::getBeta(RGset), rownames = "CGid")

# Calling preprocessNoob on RGset will return a MethylSet object, which contains 
# normalized beta values along with a few other things
Mset <- minfi::preprocessNoob(RGset)

# Calling getBeta on Mset will return a data frame of normalized beta values for 
# each CG Site
n_betas <- as_tibble(minfi::getBeta(Mset), rownames = "CGid")

# 5 Check quality of methylation data ====

# Get QC for methylset
QCset <- getQC(Mset)

# Get median unmethylated and methylated intensity (log2)
QCdf <- data.frame(mMed = QCset$mMed, uMed = QCset$uMed)

# Plot with line for "bad" sample cutoff (bad samples below line)
ggplot(QCdf, aes(x = mMed, y = uMed)) + 
  geom_point() + 
  geom_segment(x = 8, xend = 13, y = 13, yend = 8, linetype = 'dashed') + 
  xlim(8,14) + ylim(8,14)

# Do this for all pops, then plot good/bad samples for each

# 5 Transpose normalized betas for clock-fitting ====

# Transpose matrix after removing CGid column
n_betas_t <- n_betas %>%
  select(! CGid) %>%
  t() %>%
  as.data.frame()

# Add back CGid by renaming rest of columns to CG sites
colnames(n_betas_t) <- n_betas$CGid

# 6 Save for clock ====

# Save betas
saveRDS(n_betas_t, paste0('output/tbetas_PB_array', batch_no, '.rds'))
# Save normalized betas
saveRDS(n_betas, paste0('output/nbetas_PB_array', batch_no, '.rds'))
# Save sample sheet
saveRDS(sample_sheet, paste0('output/updated_sample_sheet_PB_array', batch_no, '.rds'))
