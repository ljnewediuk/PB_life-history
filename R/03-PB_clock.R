
# 04 - Polar bear clock ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Fit polar bear clock
#===============================================================================


#------------------------------------------------------------------------------
# load packages
library(tidyverse)
library(glmnet)

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
  select(sampleId, Spec, age, sex) %>%
  dplyr::rename('SampleID' = sampleId, 'Age' = age, 'Sex' = sex) %>%
  mutate(Spec = ifelse(Spec == '_Bloo', 'Blood', Spec)) %>%
  left_join(sample_sheet) %>% distinct()

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
  # Move ID columns to front
  relocate(c(SampleID, BearID, Age, Sex, Spec, chip.ID.loc)) %>%
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
