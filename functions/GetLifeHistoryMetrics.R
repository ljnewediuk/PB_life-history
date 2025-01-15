
# Function to get lifetime reproductive success and age at first reproduction
# for downstream analysis

# Need to supply:
#   - The offspring data (offspring, their birth dates, and dam and sire)
#   - Epigenetic ages (not required if calculating for the whole population)
#   - Age info (bear IDs and birth dates)

firstReproLRS <- function(cubdat, epidat = NULL, ageinfo, wholepop = T) {
  # Rename dams and sires as BearID in each data frame for grouping, and filter
  # out any rows without dam or sire
  sires <- cubdat %>%
    rename('BearID' = sire) %>%
    select(BearID, Offspring) %>%
    filter(! is.na(BearID)) %>%
    distinct()
  dams <- cubdat %>%
    rename('BearID' = dam) %>%
    select(BearID, Offspring) %>%
    filter(! is.na(BearID)) %>%
    distinct()
  # Join with offspring info for birth dates of offspring
  sires <- sires %>%
    left_join(select(offspring_info, Offspring, Born)) %>%
    distinct()
  dams <- dams %>%
    # filter(BearID %in% epidat$BearID) %>%
    left_join(select(offspring_info, Offspring, Born)) %>%
    distinct()
  # Subset just bears with epigenetic age if wholepop == F
  if(isFALSE(wholepop)) {
    dams <- dams %>%
      filter(BearID %in% epidat$BearID)
    sires <- sires %>%
      filter(BearID %in% epidat$BearID)
  }
  # Group by and get reproductive info for males and females
  mf_repro_dat <- data.frame()
  for(i in c('dams', 'sires')) {
    # Get male/female data
    dat <- get(i)
    # Calculate earliest and latest reproductions and LRS by bear
    repro_dat <- dat %>%
      group_by(BearID) %>%
      summarize(Earliest = min(Born, na.rm = T),
                Latest = max(Born, na.rm = T),
                LRS = n()) %>%
      # Add birth info for bears of interest
      left_join(age_info) %>%
      # Calculate first and last reproduction based on year of birth
      mutate(FirstRepro = Earliest - Born, LastRepro = Latest - Born) %>%
      select(BearID, Born, FirstRepro, LastRepro, LRS)
    # Filter out nonsensical ages at reproduction (before birth or after
    # age 35)
    repro_dat <- repro_dat %>%
      filter(! FirstRepro < 0 & ! LastRepro > 35)
    # Add column for sex
    if(i == 'sires') {
      repro_dat <- repro_dat %>%
        mutate(Sex = 'M')
    } else {
      repro_dat <- repro_dat %>%
        mutate(Sex = 'F')
    }
    # Add epigenetic aging data if wholepop == F
    if(isFALSE(wholepop)) {
      repro_dat <- repro_dat %>%
        left_join(select(epidat, BearID, AgeAccel))
    }
    # Combine data
    mf_repro_dat <- rbind(mf_repro_dat, repro_dat)
    
  }
  return(mf_repro_dat)
}
