
library(tidyverse)
library(FSAdata)
library(FSA)

# Load growth data
growth_dat <- read.csv('input/growth_data.csv') %>%
  mutate(AXG = as.numeric(AXG),
         SLen = as.numeric(SLen)) %>%
  # Remove error/no age data
  filter(! Age > 35, ! Born <= 1950, AXG < 280, AXG > 20, Sex %in% c('F', 'M'))

# Set function for typical vBF model
( vb <- vbFuns(param="Typical") )

# Function to model growth based on subset of years (yrs), sexes (sex), 
# growth parameters (parm), and data (growth_dat)
growth_fun <- function(yrs, sex, parm, obj_type = 'resid', growth_dat) {
  
  # Subset data
  growth_subset <- growth_dat %>%
    filter(Born %in% yrs, Sex == sex) %>%
    select(BearID, Age, Sex, Born, all_of(parm)) %>%
    na.omit()
  # Rename growth parameter to something general for model
  colnames(growth_subset)[5] <- 'GrowthParm'
  
  # Get reasonable starting values for the optimization algorithm
  ( f.starts <- vbStarts(GrowthParm ~ Age, data = growth_subset) )
  # Estimate growth parameters from data (Growth ~ parameters)
  f.fit <- nls(GrowthParm ~ vb(Age, Linf, K, t0), data = growth_subset, start = f.starts)
  
  # Get predicted growth by age
  pred_growth <- predict(f.fit, data.frame(Age = 0:30))
  
  # Either get model predictions with confidence intervals...
  # **** DOESN'T WORK... SOMETHING IN car PACKAGE INTERFERING WTIH FUNCTION ****
  if(obj_type == 'predict') {
    # Bootstrap function
    boot_fun <- function(x) predict(x, data.frame(Age = 0:30))
    # Get confidence intervals from model
    pred_growth_lower <- car::Boot(f.fit, f = boot_fun) %>% 
      confint() %>% as.data.frame() %>% pull(`2.5 %`) %>% as.numeric()
    pred_growth_upper <- car::Boot(f.fit, f = boot_fun) %>% 
      confint() %>% as.data.frame() %>% pull(`97.5 %`) %>% as.numeric()
    # Make data frame for output
    growth_out <- data.frame(Age = 0:30,
                             Cohort = paste(min(yrs), max(yrs), sep = '-'),
                             GrowthName = parm,
                             PredGrowth = pred_growth,
                             PredGrowthLower = pred_growth_lower,
                             PredGrowthUpper = pred_growth_upper)
  } 

  # Or residuals...
  if(obj_type == 'resid') {
    # Get growth residuals
    growth_out <- growth_subset %>%
      mutate(Cohort = paste(min(yrs), max(yrs), sep = '-'),
             GrowthName = parm,
             ResidGrowth = as.numeric(resid(f.fit)),
             Born = factor(Born))
  }
  
  # Return output
  return(growth_out)
  
}

# Model growth by cohort and sex, saving predictions and residuals in data frames
# List of cohorts
cohort_list <- list(1965: 1969, 1970:1974, 1975:1979, 1980:1984, 1985:1989,
                    1990:1994, 1995:1999, 2000:2004, 2005:2009, 2010:2014, 2015:2020)

growth_preds_ch <- data.frame()
for(i in cohort_list) {
  
  for(j in c('M', 'F')) {
    
    # Model
    cohort_growth <- growth_fun(yrs = i, sex = j, parm = 'AXG', 
                                obj_type = 'predict', growth_dat = growth_dat)
    # Predictions
    growth_pred <- cohort_growth %>%
      mutate(Sex = j) 
    # Bind together
    growth_preds_ch <- rbind(growth_preds_ch, growth_pred)
    
  }
}

growth_resids_ch <- data.frame()
for(i in cohort_list) {
  
  for(j in c('M', 'F')) {
    
    # Model
    cohort_growth <- growth_fun(yrs = i, sex = j, parm = 'AXG', 
                                obj_type = 'resid', growth_dat = growth_dat)
    # Bind together
    growth_resids_ch <- rbind(growth_resids_ch, cohort_growth)
    
  }
}

growth_resids_sum <- growth_resids_ch %>%
  group_by(BearID, Cohort) %>%
  summarize(MResid = mean(ResidGrowth))

ggplot(growth_preds_ch, aes(x = Age, y = PredGrowth, colour = Cohort, fill = Cohort)) +
  geom_ribbon(aes(x = Age, ymin = PredGrowthLower, ymax = PredGrowthUpper), alpha = 0.5) +
  geom_line() + facet_wrap(~ Sex)

growth_pop <- growth_fun(yrs = 1965:2018, sex = 'M', parm = 'AXG', growth_dat = growth_dat) %>%
  rbind(growth_fun(yrs = 1965:2018, sex = 'F', parm = 'AXG', growth_dat = growth_dat)) %>%
  group_by(BearID) %>%
  summarize(MRGrowthPop = mean(ResidGrowth))

lh_dat <- readRDS('output/lh_info_pop.rds') %>%
  mutate(Born = factor(Born))

temp_dat <- read.table('input/global_temp_data.txt', header = T) %>%
  mutate(Born = factor(Year)) %>%
  rename('Temp' = Lowess.5.) %>%
  select(Born, Temp)

growth_repro_ch <- growth_resids_ch %>%
  left_join(lh_dat) %>%
  left_join(temp_dat) %>%
  group_by(Cohort) %>%
  mutate(MFirstRepro = mean(FirstRepro, na.rm = T),
         sdFirstRepro = sd(FirstRepro, na.rm = T),
         ZFirstRepro = (FirstRepro - MFirstRepro)/sdFirstRepro) %>%
  # filter(Age %in% 0:5) %>%
  na.omit() %>%
  group_by(BearID, Sex, Born, ZFirstRepro, FirstRepro) %>%
    summarize(MRGrowth = mean(ResidGrowth)) %>%
  # filter(! is.na(LRS)) %>%
  mutate(NBorn = as.numeric(Born)) %>%
  left_join(growth_pop) %>%
  filter(! MRGrowth > 40) 

ggplot(growth_repro_ch, aes(x = MRGrowthPop, y = FirstRepro, colour = ch))  +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~ Sex) 
  
epi_growth <- readRDS('input/PB_clock_ages2.rds') %>%
  group_by(BearID) %>%
  summarize(MAccel = mean(AgeAccel)) %>%
  left_join(growth_repro_ch) %>%
  left_join(lh_dat) %>%
  mutate(ch = ifelse(Born %in% 1965:1974, '65s', NA),
         ch = ifelse(Born %in% 1975:1984, '75s', ch),
         ch = ifelse(Born %in% 1985:1994, '85s', ch),
         ch = ifelse(Born %in% 1995:2004, '95s', ch),
         ch = ifelse(Born %in% 2005:2014, '00s', ch),
         ch = ifelse(Born %in% 2015:2022, '10s', ch)) %>%
  na.omit()

ggplot(epi_growth, aes(y = MAccel, x = MRGrowth)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~ Sex)

ggplot(epi_growth, aes(y = MAccel, x = FirstRepro)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~ Sex)

ggplot(epi_growth, aes(y = MRGrowth, x = FirstRepro)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~ Sex)


# Save growth residuals
saveRDS(growth_resids_ch, 'output/cohort_growth_residuals.rds')
# Save models
saveRDS(growth_preds_ch, 'output/cohort_growth_data.rds')


