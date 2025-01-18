
library(tidyverse)
library(glmnet)
library(brms)

# 1 Load data and source functions ====

# Load samples that failed QC
failed_QC <- readRDS('output/failed_QC_samples.rds')

# Load related bears
sibs <- readRDS('input/full_sibs.rds')

# Load list of probes to include in clock
probes <- readRDS('output/clock_Cpgs.rds')

# Source function for cleaning the betas
source('functions/CleanBetas.R')

# Clean betas in batches 1-3 and separate into training and testing (this 
# function also removes samples that did not pass QC, excludes siblings from
# the clock training data, and subsets probes to those identified in the EWAS)
meth_dat <- cleanBetas(batches = 1:3, failed_s = failed_QC, 
                       keep_p = probes, sep_train = F)

# Set iterations
it <- 1
# Run while loop cross-validation
while(it <= 50) {

  print(it)
  
  # Specify training and testing data
  meth_betas_train <- meth_dat %>%
    group_by(Age, Spec, Sex) %>%
    sample_n(1) %>%
    filter(! BearID %in% sibs)
  
  meth_betas_test <- meth_dat %>%
    filter(! SampleID %in% meth_betas_train$SampleID)
  
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
  
  # Get betas and ages
  betasLoop <- meth_betas_train_m
  ageLoop <- meth_betas_train$Age
  
  # Make sure ages for training match training samples
  meth_betas_train$chip.ID.loc == rownames(meth_betas_train_m)
  
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
  
  # MAE
  MAE <- median(abs(age_preds$AgePredict - age_preds$Age))
  # Pearson's correlation
  corr <- as.numeric(cor.test(age_preds$AgePredict, age_preds$Age)$estimate)
  
  # Fit model age acceleration ~ year of birth
  m <-  brm(AgeAccel ~ Born,
            data = age_preds, family = gaussian, 
            iter = 10000, warmup = 5000, chains = 4, cores = 4, 
            prior = prior(normal(0,1), class = b),
            control = list(adapt_delta = 0.99, max_treedepth = 20),
            backend = 'cmdstanr')
  
  # Pull out posterior draws
  b <- as.data.frame(m) %>%
    pull(b_Born)
  
  # Add metrics and posterior draws to growing objects
  if(it == 1) {
    b_draws <- b
    mets <- data.frame(iteration = 1, MAE, corr, N = nrow(meth_betas_train))
  } else {
    b_draws <- c(b_draws, b)
    mets <- mets %>%
      bind_rows(data.frame(iteration = it, MAE, corr))
  }
  
  # Set new iterations
  it <- it + 1
  
}

saveRDS(b_draws, 'output/bootstrap_age_birth_draws.rds')
saveRDS(mets, 'output/bootstrap_clock_metrics.rds')


