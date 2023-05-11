
library(marked) 
library(tidyverse)

# Organize bear samples by cohort

# Load bear epigenetic ages
epi_ages <- readRDS('input/PB_clock_ages.rds') %>%
  mutate(Date = substr(SampleID, 8, 17))

# Load capture info and join with ages
captures <- read.csv('input/bear_capture_info.csv') %>%
  rename('BearID' = BearCode) %>%
  right_join(epi_ages) %>% 
  mutate(CaptYear = as.numeric(substr(Date, 1, 4))) %>%
  select(BearID, Spec, Born, Sex, CaptYear, VisitAge, AgeAccel) %>%
  # Fix bear with wrong date
  mutate(Born = ifelse(BearID == 'X19212', 1998, Born),
         VisitAge = ifelse(BearID == 'X19212', CaptYear - Born, VisitAge))

# Capture history for individuals up to maximum of 35 years 
# All individuals born ≤ 1986 would reach their maximum lifespan of 
# 35 years as of 2021
captures_year <- captures %>%
  filter(Born <= 1991)

# yr <- 1990
# 
# captures_year <- captures %>%
#   filter(Born == yr)

#### TURN INTO FUNCTION
capture_histories <- data.frame()
for(i in unique(captures_year$BearID)) {
  
  id_BC <- captures_year %>%
    filter(BearID == i)
  
  yr <- unique(id_BC$Born)
  
  ch <- id_BC %>%
    right_join(expand.grid(CaptYear = yr:(yr+30))) %>%
    mutate(ch = ifelse(is.na(VisitAge), 0, 1)) %>%
    select(BearID, CaptYear, ch) %>%
    distinct() %>%
    arrange(CaptYear) %>%
    pull(ch)
  
  accel <- id_BC %>%
    summarize(AgeAccel = mean(AgeAccel)) %>%
    pull(AgeAccel)
  # accel_skin <- id_BC %>%
  #   filter(Spec %in% c('Skin')) %>%
  #   summarize(AgeResid = mean(AgeResid)) %>%
  #   pull(AgeResid)
  # 
  # accel_blood <- id_BC %>%
  #   filter(Spec %in% c('Blood')) %>%
  #   summarize(AgeResid = mean(AgeResid)) %>%
  #   pull(AgeResid)
  
  ch_row <- data.frame(BearID = unique(id_BC$BearID),
                       BearCohort = unique(id_BC$Born),
                       BearSex = unique(id_BC$Sex),
                       AgeAccel = accel,
                       ch = paste(ch, collapse = ''))
  
  capture_histories <- rbind(capture_histories, ch_row)
}

# Get quantiles of age acceleration for predicting
l_quant_accel <- quantile(capture_histories$AgeAccel, probs =  0.025)
u_quant_accel <- quantile(capture_histories$AgeAccel, probs =  0.975)

capt_proc <- process.data(capture_histories)
design.Phi <- list(static = c('AgeAccel', 'BearCohort', 'BearSex'), age.bins = c(0, 10, 20, 30))
design.parameters <- list(Phi = design.Phi)

capt_ddl <- make.design.data(capt_proc, parameters = design.parameters)

names(capt_ddl$Phi)

Phi.sfw <- list(formula = ~AgeAccel*age)
mod <- crm(capt_proc, capt_ddl, model.parameters = list(Phi = Phi.sfw), hessian = T)

new_epi_ages <- expand.grid(AgeAccel = c(l_quant_accel, u_quant_accel), BearCohort = c(1965, 1985, 2020), BearSex = c('F', 'M'))

reals <- predict(mod, newdata = new_epi_ages, parameter = 'Phi', se = T) %>%
  mutate(age_cat = ifelse(age %in% 0:4, 'Cub (0-4 y)', NA),
         age_cat = ifelse(age %in% 5:20, 'Adult (5-20 y)', age_cat),
         age_cat = ifelse(age %in% 21:29, 'Old (21-30 y)', age_cat)) %>%
  mutate(age_cat = factor(age_cat, levels = c('Cub (0-4 y)', 
                                              'Adult (5-20 y)', 
                                              'Old (21-30 y)'))) %>%
  group_by(AgeAccel, age_cat) %>%
  summarize(est = mean(estimate),
            lower = mean(lcl),
            upper = mean(ucl)) %>%
  mutate(AgeAccel = ifelse(AgeAccel < 0, 'Lower % accel.', 'Upper % accel.'))

ggplot(reals, aes(x = age_cat, y = est, group = AgeAccel, col = AgeAccel)) +
  scale_colour_brewer(palette = 'Set1') +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
        axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
        axis.text = element_text(colour = 'black', size = 18),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 18, colour = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, colour = 'black'),
        legend.position = c(.2, .2)) +
  ylab('P survival') +
  xlab('')


# Phi is the probability the individual is alive at time t+1 given it is alive
# at time t. Since the time intervals are ages, maybe we can include time as 
# a covariate

