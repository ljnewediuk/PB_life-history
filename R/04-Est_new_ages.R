
library(tidyverse)

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

# Age new samples in batch 9
WH_ages <- ageNew(batch_no = 9, clock = PB_clock, failed_s = failed_QC) %>%
  select(! Batch) %>%
  rbind(first_batch_ages) 

# Get mean age accel for individuals
mean_ages <- WH_ages %>% 
  group_by(BearID) %>%
  summarize(AgeAccel = mean(AgeAccel), Born = unique(Born))

# Quick linear model
summary(lm(AgeAccel ~ Born, data = mean_ages))

# Quick plot
mean_ages %>%
  ggplot(aes(x = Born, y = AgeAccel)) +
  geom_point() + 
  geom_smooth(method = 'lm')

