
# 04 - Estimate new ages ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Estimate new ages for batch 9 containing more western Hudson Bay samples, 
#n = 96, and combine with the original test samples
#===============================================================================


#------------------------------------------------------------------------------
# load packages
library(tidyverse)

# 1 Load data and source functions ====

# Source function for aging new samples
source('functions/AgeSamples.R')

# Load PB clock
PB_clock <- readRDS('output/PB_clock.rds')

# Load WHB ages
first_batch_ages <- readRDS('output/PB_clock_ages.rds') %>%
  # Adjust columns to match spatial data
  mutate(Born = yr - floor(Age),
         YMD = as.Date(substr(SampleID, 8, 17))) %>%
  rename(Sample_Name = SampleID) %>%
  select(AgePredict, Age, AgeAccel, Sample_Name, YMD, 
         BearID, Born, Spec, Sex)

# Load samples that failed QC
failed_QC <- readRDS('output/failed_QC_samples.rds')

# 2 Age new samples in batch 9 ====

WH_ages <- ageNew(batch_no = 9, clock = PB_clock, failed_s = failed_QC) %>%
  select(! Batch) %>%
  rbind(first_batch_ages) 

# 3 Save combined ages ====

saveRDS(WH_ages, 'output/WH_combined_ages.rds')

