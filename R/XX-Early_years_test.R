
# Supplement - Early bears clock ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
# NOTE: This is only a validation of our original clock. The clock produced here
# has issues that make it unreliable.

# We are including bears born < 1986 to train the clock and bears born > 1986 to
# predict. This is a good idea but potentially problematic because there were
# fewer younger bears sampled in the early years and bears sampled in later years
# are primarily younger. We tried to avoid the age bias problem by sampling
# evenly across years, which means we have an even spread of ages in our data,
# but this is not the case when we use all the early bears to predict later bears.
# All of the early bears are mature (older than 5 yrs with most in their 20s), 
# and most later bears are  younger (under 5). The clock under-predicts the ages 
# of younger bears, creating a spurious and relatively strong positive
# association between age and age acceleration not present in our clock. Because
# of this spurious association, most of the later-born bears (which are also 
# primarily younger) show a spurious pattern of deceleration, so that the positive
# association we found between birth year and age acceleration using our much more
# accurate clock appears negative using the early-years clock.
#===============================================================================


#------------------------------------------------------------------------------
# load packages
library(tidyverse)
library(glmnet)

# Load related bears
sibs <- readRDS('input/full_sibs.rds')

# Birth year info
born_info <- read.csv('input/bear_capture_info.csv') %>%
  select(BearCode, Born) %>%
  rename('BearID' = BearCode) %>%
  distinct()

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
  select(sampleId, Spec, age, sex) %>%
  dplyr::rename('SampleID' = sampleId, 'Age' = age, 'Sex' = sex) %>%
  left_join(sample_sheet) %>% distinct()

# Check which bears have multiple samples (need to exclude these too)
multi_samps <- sample_sheet %>%
  mutate(BearID = substr(SampleID, 1, 6),
         YrSampled = substr(SampleID, 8, 11),
         Spec = substr(SampleID, 19, 20)) %>%
  left_join(born_info) %>%
  group_by(BearID) %>%
  mutate(n_samps = n()) %>%
  filter(n_samps > 1)

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
  # Join with birth year data
  left_join(born_info) %>%
  # Move ID columns to front
  relocate(c(SampleID, BearID, Age, Sex, Born, Spec, chip.ID.loc)) %>%
  # Fix typo
  mutate(Spec = ifelse(Spec == '_Bloo', 'Blood', Spec)) %>%
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

# Find earliest samples from bears with multiple
earliest_multi_samps <- multi_samps %>%
  filter(Spec == 'Sk' & YrSampled == min(YrSampled)) %>%
  pull(SampleID) %>%
  unique()

# These are all the longitudinal samples
all_multi_samps <- multi_samps %>%
  pull(SampleID) %>%
  unique()

# These are all the longitudinal samples excluding the earliest samples that we
# can use to train the clock
to_exclude <- all_multi_samps[! all_multi_samps %in% earliest_multi_samps]

# Training betas for early years
meth_betas_train <- meth_betas %>%
  # Exclude bears that are full sibs
  filter(! BearID %in% sibs) %>%
  # Exclude bears that have multiple samples (except earliest)
  filter(! SampleID %in% to_exclude) %>%
  # Exclude bears born after '85 and skin samples
  filter(Born <= 1986) %>%
  mutate(type = 'training \n (pre-85)')

# Testing betas for PB clock
meth_betas_test <- meth_betas %>%
  # Exclude bears that have multiple samples (except earliest)
  filter(! SampleID %in% to_exclude) %>%
  # Exclude bears born after '85 and skin samples
  filter(Born > 1986) %>%
  mutate(type = 'testing \n (post-85)')

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
age_df <- meth_betas %>%
  mutate(AgePredict = 0) %>%
  select(SampleID:chip.ID.loc, AgePredict)

# Get betas and ages
betasLoop <- meth_betas_train_m
ageLoop <- as.numeric(age_df[age_df$SampleID %in% meth_betas_train$SampleID ,]$Age)

# Make sure ages for training match training samples
age_df[age_df$SampleID %in% meth_betas_train$SampleID ,]$chip.ID.loc == rownames(meth_betas_train_m)

set.seed(2)

# 6 Fit clock and predict on training data ====

# Glmnet model (training betas ~ ages)
cvfit <- cv.glmnet(betasLoop, ageLoop, nfolds = 10, alpha = .5)

# Add predictions as column to ages in training data
age_preds <- age_df %>%
  filter(SampleID %in% meth_betas_test$SampleID) %>%
  # Predict model
  mutate(AgePredict = as.numeric(predict(cvfit, newx = meth_betas_test_m, 
                                         type = "response", s = "lambda.min")))
# Add residuals
age_preds <-  age_preds %>%
  mutate(AgeAccel = lm(age_preds$AgePredict ~ age_preds$Age)$residuals,
         # Add year
         yr = as.numeric(substr(SampleID, 8, 11))) %>%
  # Remove chip column
  select(! chip.ID.loc)

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
  annotate(geom = 'text', x = 23, y = 0, vjust = -2, label = mae_label, size = 6) +
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
  ggplot(aes(x = Born, y = AgeAccel, colour = Spec)) +
  scale_colour_manual(values = c('#d62d20', '#536878')) +
  geom_point() +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.position = 'none')

# Plot density plot of training/testing data from early clock
meth_betas_train %>%
  mutate(type = 'training') %>%
  bind_rows(meth_betas_test) %>%
  mutate(type = ifelse(is.na(type), 'testing', type)) %>%
  ggplot(aes(x = Age, group = type, colour = type, fill = type)) + 
  geom_density(alpha = 0.5) +
  scale_colour_manual(values = c('#f7ad0d', '#6b879b')) +
  scale_fill_manual(values = c('#f7ad0d', '#6b879b')) +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.key = element_rect(colour = 'white'),
        legend.position = c(0.85, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, colour = 'black'))

# Plot density plot of training/testing from our clock:
bind_rows(train_early, test_early) %>%
  ggplot(aes(x = Age, group = type, colour = type, fill = type)) + 
  geom_density(alpha = 0.5) +
  scale_colour_manual(values = c('#f7ad0d', '#6b879b')) +
  scale_fill_manual(values = c('#f7ad0d', '#6b879b')) +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.key = element_rect(colour = 'white'),
        legend.position = c(0.45, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, colour = 'black'))
