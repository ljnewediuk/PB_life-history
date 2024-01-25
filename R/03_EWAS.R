
library(tidyverse)
library(limma)

# 03 - Epigenome-wide association survey

# Remove any probes correlated (p < 0.05) with tissue and sex, and keep top
# ~ 10% of probes correlated with age. We used this subset of probes for 
# designing the polar bear clock (started with the Cpgs most related to age 
# and least related to sex/tissue type in bears).

#   1) Two EWAS on age with blood and skin separately
#   2) One EWAS on tissue, age, sex, adjusting p-value for either sex or tissue

# The design matrix needs to be formatted with samples in rows and coefficients 
# in columns

# 1 Load data ===

# Load sample sheets
sample_sheets <- readRDS('output/updated_sample_sheet_PB_array1.rds') %>%
  rbind(readRDS('output/updated_sample_sheet_PB_array2.rds')) %>%
  rbind(readRDS('output/updated_sample_sheet_PB_array3.rds'))

# Load sample info
sample_info <- readRDS('input/batch1_samples.rds') %>%
  rbind(readRDS('input/batch2_samples.rds')) %>%
  rbind(readRDS('input/batch3_samples.rds')) %>%
  rename('Sample_Name' = sampleId)

# 2 Subset blood and skin for limma ====

# Get subset with skin only
sample_sheet_skin <- sample_sheets %>%
  filter(! grepl('Blood', Sample_Name))

# Blood only
sample_sheet_blood <- sample_sheets %>%
  filter(! grepl('Skin', Sample_Name))

# 3 Make function to run limma model ====
fit_limma <- function(betas = sample_sheets, samples = sample_info, 
                      mod_formula = c('Spec', 'age', 'sex')) {
  
  # Combine samples into design matrix
  design_mat <- betas %>%
    left_join(samples) %>%
    # Remove samples that failed or were poor quality
    filter(! Sample_Name %in% 
             c('X09304_2001-09-06_Blood', 'X09407_1988-09-07_Blood',
               'X12606_1997-08-28_Blood', 'X10776_1997-09-06_Blood', 
               'X10228_1998-08-31_Blood', 'X03292_2001-09-15_Blood',
               'X12697_2008-09-16_Skin', 'X09365_1988-09-29_Blood',
               'X12697_1997-09-22_Blood')) %>%
    # Fix mislabeled blood sample
    mutate(Spec = ifelse(Spec == '_Bloo', 'Blood', Spec))
  
  # Betas in rows, samples in columns
  beta_dat <- readRDS('output/nbetas_PB_array1.rds') %>%
    left_join(readRDS('output/nbetas_PB_array2.rds')) %>%
    left_join(readRDS('output/nbetas_PB_array3.rds')) %>%
    # Select only chip locations also in design matrix
    select(c(CGid, design_mat$chip.ID.loc)) %>%
    column_to_rownames('CGid') %>%
    as.matrix()
  
  # Create model matrix (for sex and tissue EWAS)
  mod <- model.matrix(reformulate(mod_formula, response = NULL), data = design_mat)
  
  # Run Limma for DNAm analysis.
  cb_fit <-  lmFit(beta_dat, mod)
  cb_ebayes <- eBayes(cb_fit)
  
  # pull out pval, coefficients from limma and combine with annotation
  cb_pvals_tech <- as.data.frame(cb_ebayes$p.value) %>%
    # Rename pvalue columns
    rename_with(.fn = ~paste(., 'pval', sep = '_')) %>%
    # Rename coefficient columns
    cbind(as.data.frame(cb_ebayes$coefficients)) %>%
    rename_with(.fn = ~paste(., 'coeff', sep = '_'), .cols = ! ends_with('pval')) %>%
    # Make Cg names into column
    rownames_to_column('CGid') 
  
  # Return data frame with p-values
  return(cb_pvals_tech)
  
}

# 4 Run limma for basic EWAS on tissue and sex ====

tissue_sex_limma <- fit_limma()

# Add adjusted pvalues for variable of interest
pvals_tissue_sex <- tissue_sex_limma %>%
  mutate(SpecSkin_pval_adj = p.adjust(SpecSkin_pval, method = 'BH'),
         sexM_pval_adj = p.adjust(sexM_pval, method = 'BH')) %>%
  mutate(sig_tissue = ifelse(SpecSkin_pval_adj < 0.05, 'yes', 'no'),
         sig_sex = ifelse(sexM_pval_adj < 0.05, 'yes', 'no'))

# Get pvlaues either > 0.05 or < 0.05 (depending on whether we want the Cgs
# correlated or uncorrelated with the variable)

# Tissue
cgs_cor_w_tissue <- pvals_tissue_sex %>%
  filter(SpecSkin_pval_adj < 0.05) %>%
  pull(CGid)

# Sex
cgs_cor_w_sex <- pvals_tissue_sex %>%
  filter(sexM_pval_adj < 0.05) %>%
  pull(CGid)

# 5 Run limma for EWAS on age separately for blood and skin ====

# Blood/age
blood_age_limma <- fit_limma(betas = sample_sheet_blood, mod_formula = c('age', 'sex')) %>%
  mutate(age_pval_adj = p.adjust(age_pval, method = 'BH')) %>%
  mutate(sig = ifelse(age_pval < 10e-6, 'yes', 'no'))

# Skin/age
skin_age_limma <- fit_limma(betas = sample_sheet_skin, mod_formula = c('age', 'sex')) %>%
  mutate(age_pval_adj = p.adjust(age_pval, method = 'BH')) %>%
  mutate(sig = ifelse(age_pval < 10e-6, 'yes', 'no'))

# Make table of p values correlated with age
data.frame(tissue = c(rep('skin', times = 3), rep('blood', times = 3)),
           pval = rep(c('10e-6', '10e-7', '10e-8'), times = 2),
           n = c(nrow(skin_age_limma[skin_age_limma$age_pval < 10e-6,]),
                 nrow(skin_age_limma[skin_age_limma$age_pval < 10e-7,]),
                 nrow(skin_age_limma[skin_age_limma$age_pval < 10e-8,]),
                 nrow(blood_age_limma[blood_age_limma$age_pval < 10e-6,]),
                 nrow(blood_age_limma[blood_age_limma$age_pval < 10e-7,]),
                 nrow(blood_age_limma[blood_age_limma$age_pval < 10e-8,])))

# Use all Cgs highly correlated with age (p < 10e-6) and exclude Cgs correlated 
# (p < 0.05) with sex
Cgs_sample <- skin_age_limma %>%
  # Cgs correlated with age
  rbind(blood_age_limma) %>%
  filter(age_pval_adj < 10e-6) %>%
  # Remove Cgs correlated with sex and tissue
  filter(! CGid %in% c(cgs_cor_w_sex)) %>%
  pull(CGid)  %>% 
  unique()

# 6 Save ====

# Cgs for clock
saveRDS(Cgs_sample, 'output/clock_Cpgs.rds')

# Limmas for plotting
saveRDS(blood_age_limma, 'output/limma_blood_age.rds')
saveRDS(skin_age_limma, 'output/limma_skin_age.rds')
saveRDS(pvals_tissue_sex, 'output/limma_tissue_sex.rds')


