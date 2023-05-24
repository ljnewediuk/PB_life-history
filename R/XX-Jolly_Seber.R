
library(marked) 
library(tidyverse)

# Organize bear samples by cohort

# Load bear epigenetic ages
epi_ages <- readRDS('input/PB_clock_ages.rds') %>%
  mutate(Date = substr(SampleID, 8, 17)) %>%
  group_by(BearID) %>%
  summarize(BearSex = unique(Sex),
            MAccel = mean(AgeAccel))

# Load growth data
growth_dat <- readRDS('output/male_growth_data.rds')$residuals %>%
  rbind(readRDS('output/female_growth_data.rds')$residuals) %>%
  group_by(BearID) %>%
  summarize(MeanRGrowth = mean(ResidGrowth))

# Load life history data
lh_dat <- readRDS('output/lh_info_pop.rds')

# Load capture info and join with ages
captures <- read.csv('input/bear_capture_info.csv') %>%
  rename('BearID' = BearCode) %>%
  mutate(CaptYear = as.numeric(substr(Date, 1, 4))) %>%
  select(BearID, VisitAge, CaptYear) %>%
  left_join(lh_dat) %>%
  left_join(growth_dat) %>%
  select(BearID, Born, Sex, CaptYear, VisitAge, FirstRepro, MeanRGrowth) %>%
  # Fix bear with wrong date
  mutate(Born = ifelse(BearID == 'X19212', 1998, Born),
         VisitAge = ifelse(BearID == 'X19212', CaptYear - Born, VisitAge))

# Capture history for individuals up to maximum of 25 years 
# All individuals born ≤ 1991 would reach their maximum lifespan of 
# 25 years as of 1996
captures_year <- captures %>%
  filter(Born <= 1996)

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
    right_join(expand.grid(CaptYear = yr:(yr+25))) %>%
    mutate(ch = ifelse(is.na(VisitAge), 0, 1)) %>%
    select(BearID, CaptYear, ch) %>%
    distinct() %>%
    arrange(CaptYear) %>%
    pull(ch)
  
  # accel <- id_BC %>%
  #   summarize(AgeAccel = mean(AgeAccel)) %>%
  #   pull(AgeAccel)
  
  
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
                       FirstRepro = unique(id_BC$FirstRepro),
                       MeanRGrowth = unique(id_BC$MeanRGrowth),
                       BearSex = unique(id_BC$Sex),
                       # AgeAccel = accel,
                       ch = paste(ch, collapse = ''))
  
  capture_histories <- rbind(capture_histories, ch_row)
}

captures_prepped <- capture_histories %>%
  filter(! ch == '00000000000000000000000000') %>%
  left_join(epi_ages) %>%
  filter(! is.na(MeanRGrowth)) %>%
  # filter(BearSex == 'M') %>%
  mutate(BearCohortGrp = ifelse(BearCohort %in% 1960:1969, '60s', NA),
         BearCohortGrp = ifelse(BearCohort %in% 1970:1979, '70s', BearCohortGrp),
         BearCohortGrp = ifelse(BearCohort %in% 1980:1989, '80s', BearCohortGrp),
         BearCohortGrp = ifelse(BearCohort %in% 1990:1996, '90s', BearCohortGrp))


# Get quantiles of growth
l_growth <- quantile(captures_prepped$MeanRGrowth, probs =  0.025)
u_growth <- quantile(captures_prepped$MeanRGrowth, probs =  0.975)

capt_proc <- process.data(captures_prepped)

design.Phi <- list(static = c('BearSex', 'BearCohortGrp', 'FirstRepro', 'MeanRGrowth'), 
                   age.bins = c(0, 10, 15, 20, 25))

design.p <- list(static = c('BearSex', 'BearCohortGrp', 'FirstRepro', 'MeanRGrowth'), 
                 age.bins = c(0, 10, 15, 20, 25))

design.parameters <- list(p = design.p,
                          Phi = design.Phi)

capt_ddl <- make.design.data(capt_proc, parameters = design.parameters)

Phi.f.repro <- list(formula = ~ FirstRepro*age + BearCohortGrp)
p.f.repro <- list(formula = ~ FirstRepro*age + BearCohortGrp)
Phi.growth <- list(formula = ~ MeanRGrowth*age + BearCohortGrp)
p.growth <- list(formula = ~ MeanRGrowth*age + BearCohortGrp)

mod <- crm(capt_proc, capt_ddl, 
           model.parameters = list(Phi = Phi.growth), 
           hessian = T)

new_dat <- expand.grid(FirstRepro = c(5, 10, 15, 20), 
                       MeanRGrowth = seq(from = l_growth, to = u_growth, length.out = 4),
                       BearSex = c('F', 'M'),
                       BearCohortGrp = c('60s', '70s', '80s', '90s'))

new_epi_ages <- expand.grid(MAccel = seq(from = l_quant_accel, to = u_quant_accel, length.out = 5), BearCohort = c(1965, 1985, 2020), BearSex = c('F', 'M'))

reals <- predict(mod, newdata = new_dat, parameter = 'Phi', se = T) %>%
  filter(BearCohortGrp %in% c('70s', '80s', '90s'))
  # mutate(age_cat = ifelse(age %in% 0:4, 'Cub (0-4 y)', NA),
  #        age_cat = ifelse(age %in% 5:20, 'Adult (5-20 y)', age_cat),
  #        age_cat = ifelse(age %in% 21:29, 'Old (21-30 y)', age_cat)) %>%
  # mutate(age_cat = factor(age_cat, levels = c('Cub (0-4 y)', 
  #                                             'Adult (5-20 y)', 
  #                                             'Old (21-30 y)'))) %>%
  group_by(AgeAccel, age_cat) %>%
  summarize(est = mean(estimate),
            lower = mean(lcl),
            upper = mean(ucl)) %>%
  mutate(AgeAccel = ifelse(AgeAccel < 0, 'Lower % accel.', 'Upper % accel.'))

ggplot(reals, aes(x = age, y = estimate, col = MeanRGrowth, group = MeanRGrowth)) +
  # scale_colour_brewer(palette = 'Set1') +
  geom_smooth() +
  facet_wrap(~ BearCohortGrp, ncol = 3)
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

