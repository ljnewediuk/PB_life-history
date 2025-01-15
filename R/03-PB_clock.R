
# 03 - Polar bear clock ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Fit polar bear clock
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

# List IDs of training and testing bears
train_bears <- meth_betas_train %>%
  pull(BearID)
test_bears <- meth_betas_test %>%
  pull(BearID)

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
         Born = yr - floor(Age)) %>%
  group_by(BearID) %>%
  summarize(AgeAccel = mean(AgeAccel), Born = unique(Born))

summary(lm(AgeAccel ~ Born, data = age_preds))

# 7 Plot clock ====

# Calculate median absolute error, correlation, and label for plot
age_mae <- median(abs(age_preds$AgePredict - age_preds$Age))
age_corr <- as.numeric(cor.test(age_preds$AgePredict, age_preds$Age)$estimate)
mae_label <- paste0('mae = ', round(age_mae, 1))
corr_label <- paste0('correlation = ', round(age_corr, 2))

ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_smooth(data = age_preds, aes(x = Age, y = AgePredict), 
              colour = 'black', linewidth = 1.5, method = 'lm', se = F) +
  geom_point(data = age_preds, aes(x = Age, y = AgePredict, colour = Spec), size = 2) +
  annotate(geom = 'text', x = 23, y = 0, vjust = 2, label = mae_label, size = 6) +
  annotate(geom = 'text', x = 23, y = 0, vjust = 0, label = corr_label, size = 6) +
  scale_colour_manual(values = c('#d62d20', '#536878')) +
  ylab('Epigenetic age') + xlab('Chronological age') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, colour = 'black'))

# 8 - Plot age acceleration ~ birth year ====

age_preds %>%
  mutate(Born = yr - floor(Age)) %>%
  # filter(! SampleID == 'X10994_1988-09-23_Blood') %>%
  ggplot() +
  geom_smooth(aes(x = Born, y = AgeAccel), 
              method = 'lm', colour = 'black') +
  scale_colour_manual(values = c('#d62d20', '#536878')) +
  geom_point(aes(x = Born, y = AgeAccel, colour = Spec)) +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.position = 'none')

# 9 Get clock coefficients/bears in data frames ====

# Clock coefs
PB_clock <- as.matrix(coef(cvfit, s = 'lambda.min')) %>%
  as.data.frame() %>%
  rownames_to_column('cg') %>%
  dplyr::rename('beta' = s1) %>%
  filter(beta != 0)

# Make table of training/testing bears for supplement
supp_table <- samp_specs %>%
  mutate(ID = substr(SampleID, 1, 6),
         DateSampled = substr(SampleID, 8, 17),
         Testing = ifelse(ID %in% test_bears, 'Yes', 'No'),
         Training = ifelse(ID %in% train_bears, 'Yes', 'No')) %>%
  rename('SampleType' =  Spec) %>%
  select(ID, DateSampled, SampleType, Age, Sex, Testing, Training)

# 10 Save ====

# Clock cpgs as .rds, .csv
saveRDS(PB_clock, 'output/PB_clock.rds')
write.csv(PB_clock, 'output/PB_clock.csv')

# Save predicted ages
saveRDS(age_preds, 'output/PB_clock_ages.rds')

# Save table for supplement
write.csv(supp_table, 'output/supplementary_bear_data.csv')
