
# 05 - Calculate reproduction ====

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
# Some things of general interest (interval between births, estimated numbers
# of lost litters). Other things are necessary for analyses (age at first repro,
# lifetime repro success)

# Removed any sires/dams born before 1980 or after 2010 to 
# ensure they can be captured within the sampling window. The population was
# not comprehensively sampled prior to 1980, so more reason to suspect that 
# any bears born earlier than 1980 may not have had their first reproduction
# captured. Similarly, bears born after 2000 may not have had their first
# reproduction captured because there haven't been enough sampling years yet.

# 1 Load data ====

# Load epigenetic age estimates
epi_age <- readRDS('output/WH_combined_ages.rds')

# Load pedigree (embargoed file)
b_ped <- read.csv('input/bearPED.csv') %>%
  rename('Offspring' = id)

# Load capture info for birth years
captures <- read.csv('input/bear_capture_info.csv') 

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

# 4 Calculate age at first reproduction ====

# Function to get age at first repro for either entire population or just
# bears with epigenetic age data
firstRepro <- function(cubdat, epidat, ageinfo, wholepop = T) {
  # Rename dams and sires as BearID for grouping
  sires <- cubdat %>%
    rename('BearID' = sire) %>%
    select(BearID, Offspring) %>%
    filter(! is.na(BearID)) %>%
    distinct()
  dams <- cubdat %>%
    rename('BearID' = dam) %>%
    select(BearID, Offspring) %>%
    filter(! is.na(BearID)) %>%
    distinct()
  # Subset just bears with epigenetic age if wholepop == F
  if(wholepop == F) {
    sires <- sires %>%
      filter(BearID %in% epidat$BearID) %>%
      left_join(select(offspring_info, Offspring, Born)) %>%
      distinct()
    dams <- dams %>%
      filter(BearID %in% epidat$BearID) %>%
      left_join(select(offspring_info, Offspring, Born)) %>%
      distinct()
  }
  # Group by and get reproductive info for males and females
  mf_repro_dat <- data.frame()
  for(i in c('dams', 'sires')) {
    # Get male/female data
    dat <- get(i)
    # Earliest and latest reproductions and LRS by bear
    repro_dat <- dat %>%
      group_by(BearID) %>%
      summarize(Earliest = min(Born, na.rm = T),
                Latest = max(Born, na.rm = T),
                LRS = n()) %>%
      # Add birth info for bears of interest
      left_join(epidat) %>%
      # Calculate first and last reproduction based on year of birth
      mutate(FirstRepro = Earliest - Born, LastRepro = Latest - Born) %>%
      select(BearID, AgeAccel, Born, FirstRepro, LastRepro, LRS)
    # Filter
    # Either offspring had no birth date (results in Inf when calculating 
    # min/max birth), or the bear had no ID (causes all bears with is.na(ID) to 
    # be linked to the same offspring)
    # Also filter out nonsensical ages at reproduction (before birth or after
    # age 35)
    repro_dat <- repro_dat %>%
      # Infinites and NAs
      filter(! is.infinite(FirstRepro)) %>%
      filter(! is.na(FirstRepro)) %>%
      filter(! FirstRepro < 0 & ! LastRepro > 35)
    # Add column for sex
    if(i == 'sires') {
      repro_dat <- repro_dat %>%
        mutate(Sex = 'M')
    } else {
      repro_dat <- repro_dat %>%
        mutate(Sex = 'F')
    }
    # Combine data
    mf_repro_dat <- rbind(mf_repro_dat, repro_dat)
    
  }
  return(mf_repro_dat)
}


# Get reproductive info for the entire population
whole_population <- firstRepro(cubdat = offspring_info, 
                               epidat = epi_age, ageinfo = age_info, 
                               wholepop = T)

# And bears with age samples only
epi_population <- firstRepro(cubdat = offspring_info, 
                             epidat = epi_age, ageinfo = age_info, wholepop = F)


