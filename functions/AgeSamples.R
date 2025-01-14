
# Function to age new samples using clock of choice

# Load sample sheet, sample specs, transposed betas, then apply the clock to
# predict new epigenetic ages and calculate residuals for age acceleration.

ageNew <- function(batch_no, clock, failed_s) {
  
  # Load sample sheets for spatial bears
  sample_sheet <- readRDS(paste0('output/updated_sample_sheet_PB_array', 
                                 batch_no, '.rds')) %>%
    select(Sample_Name, chip.ID.loc) 
  
  # Load batch info
  sample_specs <- readRDS(paste0('input/batch', batch_no, '_samples.rds')) %>%
    rename(Sample_Name = sampleId) %>%
    right_join(sample_sheet) %>%
    distinct() 
  
  # Get df with DNAm of new samples at all sites in the PB clock
  meth_betas <- readRDS(paste0('output/tbetas_PB_array', batch_no, '.rds')) %>%
    # Rename for joining
    rownames_to_column('chip.ID.loc') %>%
    # Select only sites in clock
    select(c(chip.ID.loc, clock[2:nrow(clock),]$cg)) %>%
    # Add batch names
    mutate(Batch = batch_no) %>%
    relocate(Batch, .after = chip.ID.loc) 
  
  # Predict ages of new samples
  age_preds <- meth_betas %>%
    rowwise() %>%
    # Add intercept and multiply meth betas by Cpg coefficients
    mutate(
      AgePredict = clock[1,]$beta + 
        sum(c_across(starts_with('cg')) * clock[2:nrow(clock),]$beta)
    ) %>%
    select(chip.ID.loc, Batch, AgePredict) %>%
    # Join with sample information
    right_join(sample_specs) %>%
    na.omit() %>%
    # Remove failed samples
    filter(! Sample_Name %in% failed_QC)
  
  # Calculate age to 0.25 yrs, assigning a birth date of January 1 of year, 
  # then get quarterly age by 3 mos.
  birth_dates <- age_preds %>%
    mutate(birth_date = as.Date(paste0(Born, '-01-01')),
           AgeDays = as.numeric(difftime(YMD, birth_date, units = 'days')),
           CorrectedAge = AgeDays/365,
           CorrectedAgeQu = case_when(CorrectedAge - age <= 0.25 ~ age + 0.25,
                                      CorrectedAge - age > 0.25 & 
                                        CorrectedAge - age <= 0.5 ~ age + 0.5,
                                      CorrectedAge - age > 0.5 & 
                                        CorrectedAge - age <= 0.75 ~ age + 0.75,
                                      CorrectedAge - age > 0.75 & 
                                        CorrectedAge - age <= 0.99 ~ age + 1)) %>%
    select(CorrectedAgeQu, Sample_Name) %>%
    rename(Age = CorrectedAgeQu)
  
  # Data frame with age acceleration (residuals of epi ~ chron age model)
  age_preds <- age_preds %>%
    filter(! is.na(age)) %>%
    left_join(birth_dates) %>%
    select(! age)
  
  age_accels <- age_preds %>%
    cbind(AgeAccel = lm(age_preds$AgePredict ~ age_preds$Age)$residuals) %>%
    rename(BearID = id, Sex = sex) %>%
    relocate(AgePredict, .before = AgeAccel) %>%
    select(! chip.ID.loc)
  
  return(age_accels)
  
} 
