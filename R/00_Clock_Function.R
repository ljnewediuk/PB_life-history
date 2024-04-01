
# Polar bear clock function ====

# Author: Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#Function to fit polar bear clock and return a list including a data.frame with
#age predictions and age acceleration for test samples ('Predictions'), a 
#data.frame with clock CpGs and betas ('Clock CpGs'), the glmnet model ('Clock
#model'), a ggplot of epigenetic ~ chronological age ('Clock plot'), a
#data.frame with identifying information from testing and training samples 
#('Table'), the mean absolute error ('MAE'), and the correlation between 
#epigenetic and chronological age ('Correlation').
#===============================================================================


#------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(glmnet)

fit_clock <- function(train, test, betas) {
  
  # List IDs of training and testing bears
  train_bears <- train %>%
    pull(BearID)
  test_bears <- test %>%
    pull(BearID)
  
  # Get matrix of betas for training data
  meth_betas_train_m <- train %>%
    # Make chip positions rownames
    column_to_rownames('chip.ID.loc') %>%
    # Remove extra cols
    select(! SampleID:Spec) %>%
    # Convert to matrix
    as.matrix()
  
  # Get matrix of betas for test data
  meth_betas_test_m <- test %>%
    # Make chip positions rownames
    column_to_rownames('chip.ID.loc') %>%
    # Remove extra cols
    select(! SampleID:Spec) %>%
    # Convert to matrix
    as.matrix()
  
  # 5 Check ages ====
  
  # Add ages for training and testing
  age_df <- betas %>%
    mutate(AgePredict = 0) %>%
    select(SampleID:chip.ID.loc, AgePredict)
  
  # Get betas and ages
  betas_vec <- meth_betas_train_m
  age_vec <- as.numeric(age_df[age_df$SampleID %in% train$SampleID ,]$Age)
  
  # Make sure ages for training match training samples
  age_df[age_df$SampleID %in% train$SampleID ,]$chip.ID.loc == rownames(meth_betas_train_m)
  
  # 6 Fit clock and predict on training data ====
  
  # Glmnet model (training betas ~ ages)
  cvfit <- cv.glmnet(betas_vec, age_vec, nfolds = 10, alpha = .5)
  
  # Add predictions as column to ages in training data
  age_preds <- age_df %>%
    filter(SampleID %in% test$SampleID) %>%
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
  
  # Calculate mean absolute error, correlation, and label for plot
  age_mae <- sum(abs(age_preds$AgePredict - age_preds$Age))/nrow(age_preds)
  age_corr <- as.numeric(cor.test(age_preds$AgePredict, age_preds$Age)$estimate)
  mae_label <- paste0('mae = ', round(age_mae, 1))
  
  # 7 Plot clock ====
  
  clock_plot <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    geom_point(data = age_preds, aes(x = Age, y = AgePredict, colour = Spec), size = 2) +
    annotate(geom = 'text', x = 25, y = 0, label = mae_label, size = 6) +
    scale_colour_manual(values = c('#d62d20', '#536878')) +
    ylab('Epigenetic age') + xlab('Chronological age') +
    theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(),
          axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
          axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
          axis.text = element_text(colour = 'black', size = 18),
          legend.key = element_rect(colour = 'white'),
          legend.position = c(0.2, 0.8),
          legend.title = element_blank(),
          legend.text = element_text(size = 18, colour = 'black'))
  
  # 8 Get clock coefficients/bears in data frames ====
  
  # Clock coefs
  PB_clock <- as.matrix(coef(cvfit, s = 'lambda.min')) %>%
    as.data.frame() %>%
    rownames_to_column('cg') %>%
    dplyr::rename('beta' = s1) %>%
    filter(beta != 0)
  
  # Make table of training/testing bears for supplement
  supp_table <- all_betas %>%
    mutate(ID = substr(SampleID, 1, 6),
           DateSampled = substr(SampleID, 8, 17),
           Testing = ifelse(ID %in% test_bears, 'Yes', 'No'),
           Training = ifelse(ID %in% train_bears, 'Yes', 'No')) %>%
    rename('SampleType' =  Spec) %>%
    select(ID, DateSampled, SampleType, Age, Sex, Testing, Training)
  
  # Make list for output
  out <- list('Predictions' = age_preds, 
              'Clock CpGs' = PB_clock,
              'Clock model' = cvfit,
              'Clock plot' = clock_plot,
              'Table' = supp_table,
              'MAE' = age_mae,
              'Correlation' = age_corr)
  
  # Return output list
  return(out)
  
}
