
# X4 - Supplemental Mature bears clock ====

# Author: Levi Newediuk

# Fit clock with only sexually mature bears in case we are capturing different
# processes when we include young bears undergoing development. It seems the 
# predictions are slightly poorer, suggesting the aging process we're capturing 
# is likely similar across ages.

# load packages
library(tidyverse)
library(glmnet)
library(ggdist)

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

# Clean betas in batch 9
meth_b9 <- cleanBetas(batches = 9, failed_s = failed_QC, keep_p = probes, 
                      sep_train = F)

# 4 Make testing and training data ====

# Training betas
meth_betas_train <- meth_dat$train %>%
  # Filter only sexually mature bears
  filter(Age >= 5)

# Testing betas
meth_betas_test <- meth_dat$test %>%
  # Add batch 9
  rbind(meth_b9) %>%
  # Filter only sexually mature bears
  filter(Age >= 5)

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

# Add ages for training and testing
age_df <- meth_betas_train %>%
  mutate(AgePredict = 0) %>%
  select(SampleID:chip.ID.loc, AgePredict)

# Get betas and ages
betasLoop <- meth_betas_train_m
ageLoop <- as.numeric(age_df[age_df$SampleID %in% meth_betas_train$SampleID ,]$Age)

# 6 Fit clock and predict on training data ====

# Glmnet model (training betas ~ ages)
cvfit <- cv.glmnet(betasLoop, ageLoop, nfolds = 10, alpha = .5)

# Add predictions as column to ages in testing data
age_preds <- meth_betas_test %>%
  select(SampleID: Spec) %>%
  # Predict model
  mutate(AgePredict = as.numeric(predict(cvfit, newx = meth_betas_test_m, 
                                         type = "response", s = "lambda.min")))

# Get only age preds in the validation set
age_vals <- age_preds %>%
  filter(BearID %in% meth_dat$test$BearID)

# MAE
median(abs(age_vals$AgePredict - age_vals$Age))
# Pearson's correlation
as.numeric(cor.test(age_vals$AgePredict, age_vals$Age)$estimate)

# Add residuals
age_accel <-  age_preds %>%
  mutate(AgeAccel = lm(age_preds$AgePredict ~ age_preds$Age)$residuals,
         # Add year
         yr = as.numeric(substr(SampleID, 8, 11)),
         Born = yr - floor(Age))

# 7 Plot clock ====

# Calculate median absolute error, correlation, and label for plot
clock_plot <- ggplot(data = age_vals, aes(x = Age, y = AgePredict, colour = Spec)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_smooth(method = 'lm', se = F, linewidth = 1.5, colour = 'black') +
  geom_point(size = 3) +
  scale_colour_manual(values = c('#d62d20', '#536878')) +
  ylab('Epigenetic age') + xlab('Chronological age') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.position = 'none',
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylab('Epigenetic age (years)') + xlab('Chronological age (years)')

# 8 - Fit model ====

model_dat <- age_accel %>%
  group_by(BearID) %>% 
  summarize(Born = mean(Born), AgeAccel = mean(AgeAccel)) %>%
  ungroup() %>%
  mutate(across(c(AgeAccel, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

accel_mod  <- lm(AgeAccel_sc ~ Born_sc, data = model_dat)

# 8 - Plot age acceleration ~ birth year model and facet ====

# Get mean of posterior for age accel of F & M in accel ~ born models
accel_plot <- age_accel %>%
  ggplot(aes(x = Born, y = AgeAccel)) +
  stat_smooth(method = 'lm', linewidth = 0.5, fill = '#425d9c50', color = '#425d9c') +
  geom_point(colour = '#425d9c', size = 3) +
  theme(plot.margin = unit(c(1, 0.5, 1, 2), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        legend.position = 'none',
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylim(-19, 8) +
  ylab('Epigenetic age acceleration (years)') + xlab('Year of birth')

# Plot clock and age accel ~ birth year as panels
plot_grid(clock_plot, accel_plot, ncol = 2, align = 'h',
          labels = c('A', 'B'), label_size = 20)

# 10 - Save plot ====

ggsave('mature_clock.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/supplementary', dpi = 300, height = 12, width = 27, units = 'cm', bg = 'white')

