
# 03 - Polar bear clock ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Fit polar bear clock and save clock and predictions
#===============================================================================


#------------------------------------------------------------------------------
# load packages
library(tidyverse)
library(glmnet)

# 1 Load data and source functions ====

# Load samples that failed QC
failed_QC <- readRDS('output/failed_QC_samples.rds')

# Load related bears
sibs <- readRDS('input/full_sibs.rds')

# Load list of probes to include in clock
probes <- readRDS('output/clock_Cpgs.rds')

# Source function for cleaning the betas
source('functions/CleanBetas.R')

# 2 Clean data ==== 

# Clean betas in batches 1-3 and separate into training and testing (this 
# function also removes samples that did not pass QC, excludes siblings from
# the clock training data, and subsets probes to those identified in the EWAS)
meth_dat <- cleanBetas(batches = 1:3, failed_s = failed_QC, 
                       excl_oth = sibs, keep_p = probes, sep_train = T)

# Specify training and testing data
meth_betas_train <- meth_dat$train
meth_betas_test <- meth_dat$test

# List IDs of training and testing samples
train_bears <- meth_betas_train %>%
  pull(SampleID)
test_bears <- meth_betas_test %>%
  pull(SampleID)

# Add sample IDs to list
PB_clock_IDs <- list(train = train_bears, test = test_bears)

# Get matrix of betas for training data
meth_betas_train_m <- meth_betas_train %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! SampleID:Spec) %>%
  # Convert to matrix
  as.matrix()

# Get matrix of betas for test data
meth_betas_test_m <- meth_betas_test %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! SampleID:Spec) %>%
  # Convert to matrix
  as.matrix()

# 5 Check ages ====

# Add ages for training and testing
age_df <- meth_betas_train %>%
  mutate(AgePredict = 0) %>%
  select(SampleID:chip.ID.loc, AgePredict)

# Get betas and ages
betasLoop <- meth_betas_train_m
ageLoop <- as.numeric(age_df[age_df$SampleID %in% meth_betas_train$SampleID ,]$Age)

# Make sure ages for training match training samples
meth_betas_train$chip.ID.loc == rownames(meth_betas_train_m)

set.seed(2)

# 6 Fit clock and predict on training data ====

# Glmnet model (training betas ~ ages)
cvfit <- cv.glmnet(betasLoop, ageLoop, nfolds = 10, alpha = .5)

# Add predictions as column to ages in training data
age_preds <- meth_betas_test %>%
  select(SampleID:Spec) %>%
  # Predict model
  mutate(AgePredict = as.numeric(predict(cvfit, newx = meth_betas_test_m, 
                                         type = "response", s = "lambda.min")))
# Add residuals
age_preds <-  age_preds %>%
  mutate(AgeAccel = lm(age_preds$AgePredict ~ age_preds$Age)$residuals,
         # Add year
         yr = as.numeric(substr(SampleID, 8, 11)),
         Born = yr - floor(Age))

# 9 Clock coefficients in data frame ====

# Clock coefs
PB_clock <- as.matrix(coef(cvfit, s = 'lambda.min')) %>%
  as.data.frame() %>%
  rownames_to_column('cg') %>%
  dplyr::rename('beta' = s1) %>%
  filter(beta != 0)

# 10 Save clock, predicted ages, and training/testing IDs ====

# Clock cpgs as .rds, .csv
saveRDS(PB_clock, 'output/PB_clock.rds')
write.csv(PB_clock, 'output/PB_clock.csv')

# Save predicted ages
saveRDS(age_preds, 'output/PB_clock_ages.rds')

# Save clock sample IDs
saveRDS(PB_clock_IDs, 'output/PB_clock_IDs.rds')
