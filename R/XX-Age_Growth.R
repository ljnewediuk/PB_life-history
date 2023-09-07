
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)

# Try this at age 8 for males instead... see Atkinson et al. 1996 for ref that
# male bears stop growing at age 8.
growth_dat <- read.csv('input/growth_data.csv', header = T) %>%
  mutate(SLen = as.numeric(SLen)) %>%
  filter(! is.na(SLen)) %>%
  filter(Age %in% 0:8) %>%
  filter(Sex == 'M' | Age  %in% 0:5) %>%
  mutate(Age_sc = scale(Age, center = F)[, 1], SLen_sc = scale(SLen, center = F)[, 1])
# Separate out males and females
growth_dat_F <- growth_dat %>%
  filter(Sex == 'F')
growth_dat_M <- growth_dat %>%
  filter(Sex == 'M')

# growth_dat_F_sub <- growth_dat_F %>%
#   arrange(BearID) %>%
#   head(300)

# lh_dat <- readRDS('output/lh_info_pop.rds')
# 
# age_dat <- readRDS('input/PB_clock_ages2.rds')
# 
# growth_resids <- readRDS('output/cohort_growth_residuals.rds') %>%
#   filter(Age %in% 6:35) %>%
#   filter(Sex == 'F' | Age  %in% 9:35) %>%
#   select(BearID, ResidGrowth)

# Model size ~ age (polynomial model)
# Define function for model
poly_mod <- function(dat) {
  mod <- brm(SLen_sc ~ 1 + Age_sc + I(Age_sc^2) + (Age_sc + I(Age_sc^2) | BearID),
             data = dat, family = gaussian, 
             iter = 10000, warmup = 5000, chains = 4, cores = 4, 
             prior = prior(normal(0,1), class = b),
             control = list(adapt_delta = 0.99, max_treedepth = 18),
             backend = 'cmdstanr')
}

# Female bears
mod_F <- poly_mod(growth_dat_F)
# Male bears
mod_M <- poly_mod(growth_dat_M)

# Individual slopes
# Female bears
id_vars_F <- mod_F %>%
  spread_draws(b_Age_sc, r_BearID[BearID, term], ndraws = 1000) %>%
  filter(! term %in% 'Intercept') %>%
  pivot_wider(names_from = term, values_from = r_BearID)  %>%
  mutate(BearSlope = b_Age_sc + Age_sc + IAge_scE2) %>%
  group_by(BearID) %>% 
  summarize(mean_BearSlope = mean(BearSlope))

# Male bears
id_vars_M <- mod_M %>%
  spread_draws(b_Age_sc, r_BearID[BearID, term], ndraws = 1000) %>%
  filter(! term %in% 'Intercept') %>%
  pivot_wider(names_from = term, values_from = r_BearID)  %>%
  mutate(BearSlope = b_Age_sc + Age_sc + IAge_scE2) %>%
  group_by(BearID) %>% 
  summarize(mean_BearSlope = mean(BearSlope))

# # Population slopes
# fixed_vars <- mod %>%
#   spread_draws(b_Age) %>%
#   median_qi(b_Age) %>%
#   rename('slope' = b_Age) %>%
#   select(slope, .lower, .upper) %>%
#   mutate(var = 'Age',
#          response = 'Growth')

# Save model
saveRDS(mod_F, 'models/size_age_mod_F.rds')
saveRDS(mod_M, 'models/size_age_mod_M.rds')

# Save individual slopes
saveRDS(id_vars_F, 'models/size_age_id_draws_F.rds')
saveRDS(id_vars_M, 'models/size_age_id_draws_M.rds')




# Join life history and growth rates
lh_growth_dat <- lh_dat %>% 
  left_join(id_vars) %>%
  left_join(growth_resids) %>%
  mutate(Cohort = ifelse(Born %in% 1965:1969, '1965-1969', '9999'),
         Cohort = ifelse(Born %in% 1970:1974, '1970-1974', Cohort),
         Cohort = ifelse(Born %in% 1975:1979, '1975-1979', Cohort),
         Cohort = ifelse(Born %in% 1980:1984, '1980-1984', Cohort),
         Cohort = ifelse(Born %in% 1985:1989, '1985-1989', Cohort),
         Cohort = ifelse(Born %in% 1990:1994, '1990-1994', Cohort),
         Cohort = ifelse(Born %in% 1995:1999, '1995-1999', Cohort),
         Cohort = ifelse(Born %in% 2000:2004, '2000-2004', Cohort),
         Cohort = ifelse(Born %in% 2005:2009, '2005-2009', Cohort),
         Cohort = ifelse(Born %in% 2005:2020, '2005-2020', Cohort))

ggplot(lh_growth_dat, aes(y = mean_BearSlope, x = Born)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  ylab('Growth slope (length/yr)') + xlab('Year born') 

summary(lm(mean_BearSlope ~ Born, data = lh_growth_dat))

ggplot(lh_growth_dat, aes(y = FirstRepro, x = mean_BearSlope)) + 
  geom_point() + geom_smooth(method = 'lm') + 
  xlab('Growth slope (length/yr)') + ylab('Age at first reproduction') +
  facet_wrap(~ Sex)

summary(lm(FirstRepro ~ mean_BearSlope*Sex, data = lh_growth_dat))

# Join epigenetic aging and growth rates
epi_growth_dat <- age_dat %>%
  left_join(id_vars)

ggplot(epi_growth_dat, aes(y = AgeAccel, x = mean_BearSlope)) + 
  geom_point(aes(colour = Sex)) + geom_smooth(method = 'lm') + 
  xlab('growth slope') + ylab('age acceleration')

summary(lm(AgeAccel~mean_BearSlope, data  = epi_growth_dat))

