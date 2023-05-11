
library(tidyverse)


# Load sample specs (i.e., sample type)
samp_specs <- readRDS('input/batch1_samples.rds') %>%
  rbind(readRDS('input/batch2_samples.rds')) %>%
  rbind(readRDS('input/batch3_samples.rds')) %>%
  select(sampleId, Spec, age, sex) %>%
  dplyr::rename('SampleID' = sampleId, 'Age' = age, 'Sex' = sex) %>%
  mutate(SampleYear = as.numeric(substr(SampleID, 8, 11))) %>%
  mutate(BirthYear = SampleYear - Age) %>% 
  mutate(BearID = substr(SampleID, 1, 6))

# Load ice data and calculate ice-free days/year
WH_ice <- read.csv('input/WH_ice_breakup_dates.csv') %>%
  mutate(ice_free = jday_freezeup - jday_breakup)

# Function to calculate cumulative ice-free days
ice_days <- function(id, sample_data, ice_data) {
  
  dat <- sample_data %>%
    filter(BearID == id)
  
  # Calculate cumulative ice-free days over the age of the individual
  avg_cum_ice_free <- data.frame()
  for(i in 1:nrow(dat)) {
    
    ice_year <- ice_data %>%
      filter(yr %in% seq(dat[i ,]$BirthYear, dat[i, ]$SampleYear))
    
    cum_ice_free <- data.frame(SampleID = dat[i ,]$SampleID,
      ice_free = mean(ice_year$ice_free))
    
    avg_cum_ice_free <- rbind(avg_cum_ice_free, cum_ice_free)
    
  }
  
  return(avg_cum_ice_free)
  
}

# Get ice-free days by sample

ice_free_data <- data.frame() 
for(j in unique(samp_specs$BearID)) {
  ice_free <- ice_days(id = j, sample_data = samp_specs, ice_data = WH_ice)
  ice_free_data <- rbind(ice_free_data, ice_free)
}

# Save
saveRDS(ice_free_data, 'output/cum_ice_free.rds')

