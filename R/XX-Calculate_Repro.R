
library(tidyverse)

# Load epigenetic age estimates
epi_age <- readRDS('input/PB_clock_ages2.rds')

# Load pedigree
b_ped <- read.csv('input/bearPED.csv')

# Load capture info for birth years
captures <- read.csv('input/bear_capture_info.csv') 

age_info <- captures %>%
  rename('BearID' = BearCode) %>%
  select(BearID, Born) %>%
  distinct()

offspring_info <- captures %>%
  rename('id' = BearCode) %>%
  select(id, Born) %>%
  right_join(b_ped)

sires <- offspring_info %>%
  rename('BearID' = sire) %>%
  right_join(epi_age)

dams <- offspring_info %>%
  rename('BearID' = dam) %>%
  right_join(epi_age)

male_info <- sires %>%
  filter(! is.na(id)) %>%
  group_by(BearID) %>%
  select(BearID, Born, id) %>%
  distinct() %>%
  summarize(Earliest = min(Born, na.rm = T),
            Latest = max(Born, na.rm = T),
            LRS = n()) %>%
  left_join(age_info) %>%
  left_join(epi_age) %>%
  mutate(FirstRepro = Earliest - Born,
         LastRepro = Latest - Born) %>%
  select(BearID, FirstRepro, LastRepro, LRS) %>%
  distinct() %>%
  left_join(epi_age) %>%
  left_join(age_info) %>%
  mutate(NOffspringYr = ifelse(Born <= 1991, LRS/((Born + 30) - Born), LRS/(2021 - Born))) %>%
  filter(! is.infinite(FirstRepro))
  
  
female_info <- dams %>%
  filter(! is.na(id)) %>%
  group_by(BearID) %>%
  select(BearID, Born, id) %>%
  distinct() %>%
  summarize(Earliest = min(Born, na.rm = T),
            Latest = max(Born, na.rm = T),
            LRS = n()) %>%
  left_join(age_info) %>%
  mutate(FirstRepro = Earliest - Born,
         LastRepro = Latest - Born) %>%
  select(BearID, FirstRepro, LastRepro, LRS) %>%
  distinct() %>%
  left_join(epi_age) %>%
  left_join(age_info) %>%
  mutate(NOffspringYr = ifelse(Born <= 1991, LRS/((Born + 30) - Born), LRS/(2021 - Born))) %>%
  filter(! is.infinite(FirstRepro))

repro_dat <- female_info %>%
  rbind(male_info) %>%
  group_by(BearID, Spec) %>%
  mutate(MAccel = mean(AgeAccel)) %>%
  select(BearID, Sex, Spec, AgeAccel, MAccel, FirstRepro, LastRepro, LRS, NOffspringYr) %>%
  # Remove outlier
  filter(! BearID == 'X10994')

