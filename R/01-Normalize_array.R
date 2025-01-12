
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
# The function also runs a quality-control check for samples with high detection
# p-values on average, indicating poorer quality.
# Batch number needs to be specified, corresponding to both the sample sheet 
# for the batch and the folder containing the idat files for the batch. This 
# code is adapted from some version of code from the Clock Foundation and 
# A. Shafer.

# 1 Source the function to normalize the betas ====
source('functions/NormalizeBetas.R')

# 2 Install packages if not already installed ====

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

# 3 Normalize the betas ====

for(B in c(1:3, 9)) {
  normBetas(B)
}


