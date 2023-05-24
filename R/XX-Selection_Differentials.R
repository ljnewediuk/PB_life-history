
library(tidyverse)

# Load pedigree
b_ped <- read.csv('input/bearPED.csv')

# Load capture info for birth years
captures <- read.csv('input/bear_capture_info.csv') 

# Get male and female growth information 
growth_dat <- readRDS('output/cohort_growth_residuals.rds')
  # Females
F_growth <- growth_dat %>%
  filter(Sex == 'F') %>%
  # Rename for joining to offspring data
  rename('Dam' = BearID) %>%
  group_by(Dam) %>%
  # Summarize mean growth residual
  summarize(DamRGrowth = mean(ResidGrowth)) %>%
  # Get absolute value of the residual
  mutate(DamRGrowthAbs = abs(DamRGrowth))
# Males
M_growth <- growth_dat %>%
  filter(Sex == 'M') %>%
  rename('Sire' = BearID) %>%
  group_by(Sire) %>%
  summarize(SireRGrowth = mean(ResidGrowth)) %>%
  mutate(SireRGrowthAbs = abs(SireRGrowth))

capt_by_year_F <- captures %>%
  rename('Dam' = BearCode) %>%
  mutate(CaptYear = substr(Date, 1, 4)) %>%
  left_join(F_growth) %>%
  filter(! is.na(DamRGrowth)) 

capt_by_year_M <- captures %>%
  rename('Sire' = BearCode) %>%
  mutate(CaptYear = substr(Date, 1, 4)) %>%
  left_join(M_growth) %>%
  filter(! is.na(SireRGrowth))

# Get offspring birth dates and dam and sire for all cubs
offspring_info <- captures %>%
  rename('id' = BearCode) %>%
  select(id, Born) %>%
  right_join(b_ped) %>%
  rename('Offspring' = id, 'Sire' = sire, 'Dam' = dam) %>%
  # Join with growth data
  left_join(F_growth) %>%
  left_join(M_growth) %>%
  distinct()

# Function to calculate selection differential
sel_diff <- function(year, sex, mod_type = 'growth', abs = F) {
  
  # Subset to either male or female data
  if(sex == 'F') {
    capt_sex <- capt_by_year_F %>%
      rename('BearID' = Dam, 'RGrowth' = DamRGrowth, 'RAbs' = DamRGrowthAbs)
    birth_sex <- offspring_info %>%
      select(Dam, Born, DamRGrowth, DamRGrowthAbs, Offspring) %>%
      rename('BearID' = Dam, 'RGrowth' = DamRGrowth, 'RAbs' = DamRGrowthAbs)
  } else {
    capt_sex <- capt_by_year_M %>%
      rename('BearID' = Sire, 'RGrowth' = SireRGrowth, 'RAbs' = SireRGrowthAbs)
    birth_sex <- offspring_info %>%
      select(Sire, Born, SireRGrowth, SireRGrowthAbs, Offspring) %>%
      rename('BearID' = Sire, 'RGrowth' = SireRGrowth, 'RAbs' = SireRGrowthAbs)
  }
  
  if(abs == T) {
    capt_sex <- capt_sex %>%
      mutate(GrowthParam = RAbs)
    birth_sex <- birth_sex %>%
      mutate(GrowthParam = RAbs)
  } 
  if(abs == F) {
    capt_sex <- capt_sex %>%
      mutate(GrowthParam = RGrowth)
    birth_sex <- birth_sex %>%
      mutate(GrowthParam = RGrowth)
  }
  
  # Subset birth data to year
  birth_sub <- birth_sex %>%
    filter(Born == year) %>%
    # Summarize number of cubs born per bear
    group_by(BearID) %>%
    summarize(GrowthParam = unique(GrowthParam),
              NCubs = n()) %>%
    mutate(PCub = 1) %>%
    na.omit()
  
  # Subset capture data to year
  capt_sub <- capt_sex %>%
    filter(CaptYear == year,
           VisitAge %in% 5:20,
           ! BearID %in% birth_sub$BearID) %>%
    mutate(NCubs = 0, PCub = 0) %>%
    select(BearID, GrowthParam, NCubs, PCub) %>%
    rbind(birth_sub) %>%
    left_join(lh_info_pop)
  
  if(nrow(capt_sub) < 5) stop('No data')
  
  if(mod_type == 'growth') {
    
    poiss_mod <- glm(NCubs ~ GrowthParam, data = capt_sub, family = poisson)
    bin_mod <- glm(PCub ~ GrowthParam, data = capt_sub, family = binomial)
    
  }
  
  if(mod_type == 'afr') {
    
    capt_sub_afr <- capt_sub %>% filter(! is.na(FirstRepro))
  
    if(nrow(capt_sub_afr) < 5) stop('No data')
    
    poiss_mod <- glm(NCubs ~ FirstRepro, data = capt_sub_afr, family = poisson)
    bin_mod <- glm(PCub ~ FirstRepro, data = capt_sub_afr, family = binomial)
    
  }
  
  coef_list <- list(poisson = as.numeric(coef(poiss_mod)[2]), 
                    binomial = as.numeric(coef(bin_mod)[2]))
  
  return(coef_list)
  
}

sds_M <- data.frame()
for(y in 1967:2017) {
  sd_calc <- try(sel_diff(year = y, sex = 'M', abs = F, mod_type = 'afr'))
  if(! is.list(sd_calc)) next
  sd_row <- data.frame(Year = y, 
                       SelectionDiff_p = sd_calc$poisson,
                       SelectionDiff_b = sd_calc$binomial)
  sds_M <- rbind(sds_M, sd_row)
}

sds_F <- data.frame()
for(y in 1967:2017) {
  sd_calc <- try(sel_diff(year = y, sex = 'F', abs = T))
  if(! is.list(sd_calc)) next
  sd_row <- data.frame(Year = y, 
                       SelectionDiff_p = sd_calc$poisson,
                       SelectionDiff_b = sd_calc$binomial)
  sds_F <- rbind(sds_F, sd_row)
}


