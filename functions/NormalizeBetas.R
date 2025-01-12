
# Function to normalize betas and identify poor-quality samples

# Load sample sheet, idat files, then normalize the betas by batch. Remove 
# samples with high detection p-values. Save the updated sample sheet and the 
# normalized betas in two formats: long by wide and wide by long.

normBetas <- function(batch_no) {
  
  # 1 Load data ====
  
  # Load sample sheet
  sample_sheet_file_name <- paste0('input/PB_array', batch_no, '_sample_sheet', batch_no, '.rds')
  
  # 2 Update sample sheet ====
  
  # List idat file names to remove any chip positions without data
  chip.IDs <- substr(list.files(paste0('iscans/batch', batch_no, '/'), recursive = T, full.names = F), 1, 12)
  chip.positions <- substr(list.files(paste0('iscans/batch', batch_no, '/'), recursive = T, full.names = F), 27, 32)
  
  # Remove any files without both red and green (i.e., two files)
  run_samples <- data.frame(chip.ID.loc = paste(chip.IDs, chip.positions, sep = '_')) %>%
    group_by(chip.ID.loc) %>%
    summarize(n_samples = n()) %>%
    filter(! n_samples == 1)
  
  # Sample sheet needs a column 'Basename' pointing to the basename of a two-colour
  # .idat file (i.e., either _Red.idat or _Grn.idat). Need to add the col Basename
  # with the file path for each sample in the corresponding row. Then use 
  # 'read.metharray.exp'to find the corresponding files using the sample sheet.
  sample_sheet <- readRDS(sample_sheet_file_name) %>%
    mutate(chip.ID.loc = paste(chip.ID, stripe, sep = '_'),
           # Add basenames (i.e., file paths for iscan files)
           Basename = paste0('iscans/batch', 
                             batch_no, '/', chip.ID, '/', chip.ID, '_', stripe)) %>%
    # Filter only iscans with both red and green iscan files
    filter(chip.ID.loc %in% run_samples$chip.ID.loc)
  
  # 3 Create RG channel set ====
  
  # Create an RGChannelSet object containing raw red green channel data from .idat
  RGset <- minfi::read.metharray.exp(base = NULL, targets = sample_sheet, recursive = T) 
  
  # Annotate the RGset object with probe coordinates. 
  # This line is currently just a place holder as the annotation is empty, 
  # but needs to be run anyway.
  RGset@annotation <- c(array = 'HorvathMammalMethylChip40', annotation = "test.unknown")
  
  # 4 Run quality control check ====
  
  # Detection p-value compares the methylated and unmethylated channels to
  # background signals for every site. Large detection p-values indicate 
  # poor-quality samples.
  dpvals <- detectionP(RGset)
  
  # Get mean detection p-values
  mean_dpvals <- apply(dpvals, 2, mean) %>%
    as.data.frame() %>%
    rownames_to_column()
  
  # Make data.frame for saving and plotting
  colnames(mean_dpvals) <- c('chip.ID.loc', 'detection_p')
  
  # Make list of samples to remove
  # NOTE: This step  needs to be CHECKED manually. Look for samples with detection
  # p-values orders of magnitude greater than the rest.
  ord_dpvals <- mean_dpvals %>%
    arrange(desc(detection_p))
  
  # 5 Normalize betas ====
  
  # Calling preprocessNoob on RGset will return a MethylSet object, which contains 
  # normalized beta values along with a few other things
  Mset <- minfi::preprocessNoob(RGset)
  
  # Calling getBeta on Mset will return a data frame of normalized beta values for 
  # each CG Site
  n_betas <- as_tibble(minfi::getBeta(Mset), rownames = "CGid") 
  
  # 6 Transpose normalized betas for clock-fitting ====
  
  # Transpose matrix after removing CGid column
  n_betas_t <- n_betas %>%
    select(! CGid) %>%
    t() %>%
    as.data.frame()
  
  # Add back CGid by renaming rest of columns to CG sites
  colnames(n_betas_t) <- n_betas$CGid
  
  # 7 Save for clock ====
  
  # Save betas
  saveRDS(n_betas_t, paste0('output/tbetas_PB_array', batch_no, '.rds'))
  # Save normalized betas
  saveRDS(n_betas, paste0('output/nbetas_PB_array', batch_no, '.rds'))
  # Save sample sheet
  saveRDS(sample_sheet, paste0('output/updated_sample_sheet_PB_array', batch_no, '.rds'))
  # Save detection p-values
  saveRDS(ord_dpvals, paste0('output/detection_p_batch', batch_no, '.rds'))
}
