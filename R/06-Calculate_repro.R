
# 06 - Calculate reproduction ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Calculate lifetime reproductive success and life history data
#===============================================================================

# NOTE: bear_capture_info.csv and bearPED.csv are embargoed files due to 
# restrictions on data sharing. The script will break without these files, but
# lh_info_pop.rds and lh_info_epi.rds are available for running subsequent scripts.

#------------------------------------------------------------------------------
#load packages
library(tidyverse)

# Calculate some metrics related to life history and reproduction from pedigree.

# 1 Load data and source functions ====

# Load epigenetic age estimates
epi_age <- readRDS('output/WH_combined_ages.rds')

# Load pedigree (embargoed file)
b_ped <- read.csv('input/bearPED.csv') %>%
  rename('Offspring' = id)

# Load capture info for birth years
captures <- read.csv('input/bear_capture_info.csv') 

# Load function to calculate life history metrics
source('functions/GetLifeHistoryMetrics.R')

# 2 Organize data ====

# Get age info
age_info <- captures %>%
  rename('BearID' = BearCode) %>%
  select(BearID, Born) %>%
  distinct()

# Join with pedigree
offspring_info <- captures %>%
  rename('Offspring' = BearCode) %>%
  select(Offspring, Born) %>%
  right_join(b_ped) 

# 4 Calculate age at first reproduction and lifetime reproductive success ====

# Get reproductive info for the entire population
whole_population <- firstReproLRS(cubdat = offspring_info, 
                                  ageinfo = age_info, wholepop = T)

# And bears with age samples only
epi_population <- firstReproLRS(cubdat = offspring_info, epidat = epi_age, 
                                ageinfo = age_info, wholepop = F)

# 5 Save data ====
saveRDS(whole_population, 'output/lh_info_pop.rds')
saveRDS(epi_population, 'output/lh_info_epi.rds')
