
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
    # boot_fun <- function(x) predict(x, data.frame(Age = 0:30))
    # # Get confidence intervals from model
    # pred_growth_lower <- car::Boot(f.fit, f = boot_fun) %>% 
    #   confint() %>% as.data.frame() %>% pull(`2.5 %`) %>% as.numeric()
    # pred_growth_upper <- car::Boot(f.fit, f = boot_fun) %>% 
    #   confint() %>% as.data.frame() %>% pull(`97.5 %`) %>% as.numeric()
    # Make data frame for output
    growth_out <- data.frame(Age = 0:30,
                             Cohort = paste(min(yrs), max(yrs), sep = '-'),
                             GrowthName = parm,
                             PredGrowth = pred_growth
                             # PredGrowthLower = pred_growth_lower,
                             # PredGrowthUpper = pred_growth_upper
                             )
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

# Model growth (length) by cohort and sex, saving predictions and residuals in data frames
# List of cohorts
cohort_list <- list(1965: 1969, 1970:1974, 1975:1979, 1980:1984, 1985:1989,
                    1990:1994, 1995:1999, 2000:2004, 2005:2009, 2010:2014, 2015:2020)

growth_preds_ch <- data.frame()
for(i in cohort_list) {
  
  for(j in c('M', 'F')) {
    
    # Model
    cohort_growth <- growth_fun(yrs = i, sex = j, parm = 'SLen', 
                                obj_type = 'predict', growth_dat = growth_dat)
    # Predictions
    growth_pred <- cohort_growth %>%
      mutate(Sex = j) 
    # Bind together
    growth_preds_ch <- rbind(growth_preds_ch, growth_pred)
    
  }
}

# Residual growth in axillary girth width
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

ggplot(growth_preds_ch, aes(x = Age, y = PredGrowth, colour = Cohort)) +
  scale_colour_manual(values = c('#d53d9b', '#cd2c8f', '#b82881', '#a32372', '#8e1e63',
                                 '#791a55', '#641546', '#5a123f', '#501038', '#460e31', '#200716')) +
  # geom_ribbon(aes(x = Age, ymin = PredGrowthLower, ymax = PredGrowthUpper), alpha = 0.5) +
  geom_line(size = 1) + facet_wrap(~ Sex) +
  theme(legend.position = c(.5, .5),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
        plot.margin = unit(c(0.5, 1, 1, 1), 'cm'),
        axis.text.y.left = element_text(size = 18, colour = 'black'),
        axis.text.y.right = element_text(size = 18, colour = '#cc181e'),
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.ticks.y.right = element_line(colour = '#cc181e'),
        axis.title.y.left = element_text(size = 18, colour = 'black', vjust = 5),
        axis.title.y.right = element_text(size = 18, colour = '#cc181e', vjust = 5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, colour = 'black'))

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

growth_rate_ch <- data.frame()
for(i in unique(growth_resids_ch$BearID)) {
  
  growth_rate_id <- growth_resids_ch %>%
    filter(BearID == i)
  
  if(nrow(growth_rate_id) < 2) next
  
  min_age <- min(growth_rate_id$Age)
  max_age <- max(growth_rate_id$Age)
  
  if(min_age > 5) next
  
  growth_early <- growth_rate_id %>% filter(Age <= 5) %>%
    summarize(MEarly = mean(ResidGrowth, na.rm = T)) %>% 
    pull(MEarly)
  
  if(max_age < 20) {
    
    growth_rate_row <- data.frame(BearID = i, 
                                  GrowthEarly = growth_early,
                                  GrowthRateYear = NA)
    
  } else {
    
    growth_d <- growth_rate_id[growth_rate_id$Age == max_age ,]$ResidGrowth -
      growth_rate_id[growth_rate_id$Age == min_age ,]$ResidGrowth
    growth_rt  <- growth_d/(max_age - min_age)
    
    growth_rate_row <- data.frame(BearID = i, 
                                  GrowthEarly = growth_early,
                                  GrowthRateYear = growth_rt)
    
  }
    
  growth_rate_ch <- rbind(growth_rate_ch, growth_rate_row)
  
}

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
  left_join(growth_rate_ch) %>%
  filter(! MRGrowth > 40) 

growth_repro_M <- growth_repro_ch %>%
  filter(Sex == 'M')

# Fit model to check how age at first repro has changed with body size for males over time
frg_mod <- lm(ZFirstRepro~MRGrowthPop*NBorn, data = growth_repro_M)
# Predict model on new data
frg_new_dat <- expand.grid(MRGrowthPop = seq(from = min(growth_repro_M$MRGrowthPop),
                                             to = max(growth_repro_M$MRGrowthPop), length.out = 15),
                           NBorn = seq(from = 1, to = 40, by = 1)) %>%
  # Filter only first and last years
  filter(NBorn %in% c(1, 40))
frg_pred <- predict(frg_mod, frg_new_dat, interval = 'prediction') %>%
  as.data.frame() %>% cbind(frg_new_dat)

# Shift from the faster-growing males reproducing earlier in the early years (probably
# early-reproducers also grew faster and were of better quality overall) to the 
# SLOWER-GROWING MALES REPRODUCING EARLIER IN LATER YEARS, SUGGESTING A TRADE-OFF
# WHERE MALES GROW TO A SMALLER SIZE IF THEY REPRODUCE EARLIER. Two overall strategies,
# and since males are getting smaller and reproducing earlier on average over time,
# suggests possible selection for investing more in reproduction early on.
ggplot(frg_pred, aes(x = MRGrowthPop, y = fit)) +
  geom_ribbon(aes(x = MRGrowthPop, ymin = lwr, ymax = upr, fill = factor(NBorn)), colour = NA, alpha = 0.5) +
  geom_line(aes(colour = factor(NBorn))) 

ggplot(growth_repro_ch, aes(x = MRGrowthPop, y = ZFirstRepro, group = Born, colour = Born))  +
  geom_point(colour = 'black') +
  geom_smooth(method = 'lm', se = F) +
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
  filter(! is.na(Sex))

ggplot(epi_growth, aes(y = MAccel, x = MRGrowth)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~ Sex) +
  ylab('Age acceleration') + xlab('Residual axillary girth')

ggplot(epi_growth, aes(y = MAccel, x = FirstRepro)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~ Sex) +
  ylab('Age acceleration') + xlab('Age at first reproduction')

ggplot(epi_growth, aes(y = MRGrowth, x = FirstRepro)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~ Sex)


# Save growth residuals
saveRDS(growth_resids_ch, 'output/cohort_growth_residuals.rds')
# Save models
saveRDS(growth_preds_ch, 'output/cohort_growth_data.rds')


