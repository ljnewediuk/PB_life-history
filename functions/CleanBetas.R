
# Function to clean betas

# Load updated sample sheet, sample specs, and transposed normalized betas.
# Requires a list of failed samples ('failed_s'; any samples that failed 
# quality control), other samples to include ('excl_oth'; e.g. siblings),
# and probes to include ('keep_p'; usually the probes selected with EWAS)

cleanBetas <- function(batches, failed_s, excl_oth, keep_p, sep_train = F) {
  
  # Functions to load sample sheets, sample specs, and transposed betas
  load_sheets <- function(B) {
    S <- readRDS(paste0('output/updated_sample_sheet_PB_array', B, '.rds'))
    return(S)
  }
  load_specs <- function(B) {
    S <- readRDS(paste0('input/batch', B, '_samples.rds')) 
    return(S)
  }
  load_betas <- function(B) {
    S <- readRDS(paste0('output/tbetas_PB_array', B, '.rds'))
    return(S)
  }
  
  # Load sample sheets
  sample_sheet <- lapply(batches, function(x) load_sheets(x)) %>%
    bind_rows() %>%
    select(Sample_Name, chip.ID.loc) %>%
    dplyr::rename('SampleID' = Sample_Name) %>%
    distinct()
  
  # Load sample specs (i.e., sample type)
  sample_specs <- lapply(batches, function(x) load_specs(x)) %>%
    bind_rows() %>%
    select(sampleId, Spec, age, sex) %>%
    dplyr::rename('SampleID' = sampleId, 'Age' = age, 'Sex' = sex) %>%
    mutate(Spec = ifelse(Spec == '_Bloo', 'Blood', Spec)) %>%
    left_join(sample_sheet) %>% distinct()
  
  # Check which bears have multiple samples for excluding from training data
  multi_samps <- sample_sheet %>%
    mutate(BearID = substr(SampleID, 1, 6)) %>%
    group_by(BearID) %>%
    summarize(n_samps = n()) %>%
    filter(n_samps > 1) %>%
    pull(BearID)
  
  # Calculate age to 0.25 yrs, assigning a birth date of January 1 of year, 
  # then get quarterly age by 3 mos.
  birth_dates <- sample_specs %>%
    mutate(birth_year = as.numeric(substr(SampleID, 8, 11)) - Age) %>%
    mutate(sample_date = as.Date(substr(SampleID, 8, 17)),
           birth_date = as.Date(paste0(birth_year, '-01-01')),
           AgeDays = as.numeric(difftime(sample_date, birth_date, units = 'days')),
           CorrectedAge = AgeDays/365,
           CorrectedAgeQu = case_when(CorrectedAge - Age <= 0.25 ~ Age + 0.25,
                                      CorrectedAge - Age > 0.25 & 
                                        CorrectedAge - Age <= 0.5 ~ Age + 0.5,
                                      CorrectedAge - Age > 0.5 & 
                                        CorrectedAge - Age <= 0.75 ~ Age + 0.75,
                                      CorrectedAge - Age > 0.75 & 
                                        CorrectedAge - Age <= 0.99 ~ Age + 1)) %>%
    select(! c(Age, birth_year:CorrectedAge)) %>%
    dplyr::rename('Age' = CorrectedAgeQu)
  
  # Load sample matrices
  
  # Combine betas and exclude outliers
  meth_betas <- lapply(batches, function(x) load_betas(x)) %>%
    bind_rows() %>%
    rownames_to_column('chip.ID.loc') %>%
    left_join(birth_dates) %>%
    # Create column with ID to exclude related bears
    mutate(BearID = substr(SampleID, 1, 6)) %>%
    # Move ID columns to front
    relocate(c(SampleID, BearID, Age, Sex, Spec, chip.ID.loc)) %>%
    # Exclude failed samples (failed scans/wierd samples)
    filter(! SampleID %in% failed_s) %>%
    # Select only probes we want to include based on EWAS
    select(c(SampleID:chip.ID.loc, all_of(keep_p)))
  
  # If sep_train == TRUE, separate into testing and training data by removing
  # siblings and bears with multiple samples for training set, and return the 
  # testing and training data. Otherwise, just return all of the data
  
  if(isTRUE(sep_train)) {
    
    # Training betas
    meth_betas_train <- meth_betas %>%
      # Exclude other samples
      filter(! BearID %in% excl_oth) %>%
      # Exclude bears that have multiple samples
      filter(! BearID %in% multi_samps)
    
    # Testing betas
    meth_betas_test <- meth_betas %>%
      # Exclude bears in test set
      filter(! BearID %in% meth_betas_train$BearID)
    
    return(list(train = meth_betas_train, test = meth_betas_test))
    
  } else {
    
    return(meth_betas)
    
  }
  
}