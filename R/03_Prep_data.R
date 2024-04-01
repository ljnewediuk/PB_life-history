# 03 - Prep data ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Prep data for polar bear clock
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)

# We have an initial sample of 288 bears. We are removing specific individuals
# that might bias the clock (failed/poor quality samples, related individuals) 
# and sub-setting our initial probe search space to those correlated with age
# and uncorrelated with sex.

# 1 Load data ====

# Load related bears
sibs <- readRDS('input/full_sibs.rds')

# Load list of probes to include in clock
probes <- readRDS('output/clock_Cpgs.rds')

# Load chip sample sheets
sample_sheet <- readRDS('output/updated_sample_sheet_PB_array1.rds') %>%
  rbind(readRDS('output/updated_sample_sheet_PB_array2.rds')) %>%
  rbind(readRDS('output/updated_sample_sheet_PB_array3.rds')) %>%
  select(Sample_Name, chip.ID.loc) %>%
  dplyr::rename('SampleID' = Sample_Name) %>%
  distinct()

# Load sample specs (i.e., sample type)
samp_specs <- readRDS('input/batch1_samples.rds') %>%
  rbind(readRDS('input/batch2_samples.rds')) %>%
  rbind(readRDS('input/batch3_samples.rds')) %>%
  select(sampleId, Spec, age, sex, Born) %>%
  dplyr::rename('SampleID' = sampleId, 'Age' = age, 'Sex' = sex) %>%
  left_join(sample_sheet, relationship = 'many-to-many') %>% distinct()

# Check which bears have multiple samples (need to exclude these too)
multi_samps <- sample_sheet %>%
  mutate(BearID = substr(SampleID, 1, 6)) %>%
  group_by(BearID) %>%
  summarize(n_samps = n()) %>%
  filter(n_samps > 1) %>%
  pull(BearID)

# 2 Calculate age to 0.25 yr ====

# Calculate the actual age of the bear based on date of sample - age, assigning
# a birth date of January 1 of year, then get quarterly age by 3 mos.
birth_dates <- samp_specs %>%
  mutate(birth_year = as.numeric(substr(SampleID, 8, 11)) - Age) %>%
  mutate(sample_date = as.Date(substr(SampleID, 8, 17)),
         birth_date = as.Date(paste0(birth_year, '-01-01')),
         AgeDays = as.numeric(difftime(sample_date, birth_date, units = 'days')),
         CorrectedAge = AgeDays/365,
         CorrectedAgeQuarterly = ifelse(CorrectedAge - Age <= 0.25, Age + 0.25, 0),
         CorrectedAgeQuarterly = ifelse(CorrectedAge - Age > 0.25 & 
                                          CorrectedAge - Age <= 0.5, Age + 0.5, 
                                        CorrectedAgeQuarterly),
         CorrectedAgeQuarterly = ifelse(CorrectedAge - Age > 0.5 & 
                                          CorrectedAge - Age <= 0.75, Age + 0.75, 
                                        CorrectedAgeQuarterly),
         CorrectedAgeQuarterly = ifelse(CorrectedAge - Age > 0.75 & 
                                          CorrectedAge - Age <= 0.99, Age + 1, 
                                        CorrectedAgeQuarterly)) %>%
  select(! c(Age, birth_year:CorrectedAge)) %>%
  dplyr::rename('Age' = CorrectedAgeQuarterly)

# 3 Load sample matrices ====

# Combine matrices fro first three arrays and exclude outliers
meth_betas <- readRDS('output/tbetas_PB_array1.rds') %>%
  rbind(readRDS('output/tbetas_PB_array2.rds')) %>%
  rbind(readRDS('output/tbetas_PB_array3.rds')) %>%
  rownames_to_column('chip.ID.loc') %>%
  left_join(birth_dates) %>%
  # Create column with ID to exclude related bears
  mutate(BearID = substr(SampleID, 1, 6)) %>%
  # Fix one misspelled blood sample
  mutate(Spec = ifelse(Spec == '_Bloo', 'Blood', Spec)) %>%
  # Make column for year sampled
  mutate(YrSampled = str_sub(SampleID, 8, 11)) %>%
  # Move ID columns to front
  relocate(c(SampleID, BearID, YrSampled, Age, Sex, Spec, Born, chip.ID.loc)) %>%
  # Exclude outlier samples (failed scans/wierd samples)
  # ßs cluster away from others in array 1: 
  #   - X10228_1998_08-31_Blood
  #   - X10776_1997-09-06_Blood
  # ßs cluster away from others in array 2:
  #   - X03367_1997-09-22_Blood
  # ßs cluster away from others in array 2:
  #   - X12697_2008-09-16_Skin
  #   - X12697_1997_09_22_Blood
  #   - X09365_1988-09-29_Blood
  #   - X03292_2001-09-15_Blood
filter(! SampleID %in% c('X09304_2001-09-06_Blood', 'X09407_1988-09-07_Blood',
                         'X12606_1997-08-28_Blood', 'X10776_1997-09-06_Blood', 
                         'X10228_1998-08-31_Blood', 'X03292_2001-09-15_Blood',
                         'X12697_2008-09-16_Skin', 'X09365_1988-09-29_Blood',
                         'X12697_1997-09-22_Blood')) %>%
  # Select only probes we want to include based on EWAS
  select(c(SampleID:chip.ID.loc, all_of(probes)))

# 4 Make testing and training data ====

# Training betas
meth_betas_train <- meth_betas %>%
  # Exclude bears that are full sibs
  filter(! BearID %in% sibs) %>%
  # Exclude bears that have multiple samples
  filter(! BearID %in% multi_samps)

# Testing betas
meth_betas_test <- meth_betas %>%
  # Exclude bears in test set
  filter(! BearID %in% meth_betas_train$BearID)

# 4 - Save data for clock ====

saveRDS(meth_betas, 'output/prepped_all_data.rds')
saveRDS(meth_betas_train, 'output/prepped_train_data.rds')
saveRDS(meth_betas_test, 'output/prepped_test_data.rds')
