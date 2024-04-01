
# 04 - Polar bear clock ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Fit polar bear clock
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)

#get clock function
source('R/00_Clock_Function.R')

# We fit a glmnet model to the remaining unbiased
# samples (n = 144) and used the remaining (excluding the poor quality samples) 
# for our biological questions.

# 1 Load prepped data ====

test_betas <- readRDS('output/prepped_test_data.rds')
train_betas <- readRDS('output/prepped_train_data.rds')
all_betas <- readRDS('output/prepped_all_data.rds')

# 2 Fit clock with all data ====

full_clock <- fit_clock(train = train_betas, test = test_betas, betas = all_betas)

# 3 Fit clock using only bears born before 1985 ====

# Filter training and testing data
train_early <- train_betas %>%
  filter(! YrSampled > 1995 & ! Born > 1985)
test_early <- test_betas %>%
  filter(! Born <= 1985 & ! YrSampled <= 1995)

# Fit clock
early_clock <- fit_clock(train = train_early, test = test_early, betas = all_betas)

# 4 Fit clock using only adult bears (> 5 years) ====

# Filter training and testing data
train_adult <- train_betas %>%
  # Exclude any bears < 5 yrs
  filter(! Age < 5)
test_adult <- test_betas %>%
  # Exclude any bears < 5 yrs
  filter(! Age < 5)

# Fit clock
mature_clock <- fit_clock(train = train_adult, test = test_adult, betas = all_betas)

# 5 Save clock with all data ====

# Clock model
saveRDS(full_clock$`Clock model`, 'output/PB_clock_mod.rds')
saveRDS(full_clock$`Clock CpGs`, 'output/PB_clock.rds')

# Save predicted ages
saveRDS(full_clock$Predictions, 'output/PB_clock_ages.rds')

