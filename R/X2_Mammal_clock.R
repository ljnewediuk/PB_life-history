
# X2 - Supplemental universal clock ====

# Author: Levi Newediuk

# Test universal clock on polar bear samples to confirm the predictions are
# similar to our clock (should be more variation in the UC clock because it
# includes multiple species, tissues, etc.)

# Original code comes from S1-S3 tables of Lu et al. 2021 preprint on bioRxiV
# (https://doi.org/10.1101/2021.01.18.426733)

# The commenting wasn't great, so not 100% clear on each function, but details 
# are in the pre-print. I made some changes for clarity, conciseness, or to fix 
# things that weren't working property.

# Load packages
library(tidyverse)
library(glmnet)
library(cowplot)

# 1 Load data ====

# Load samples that did not pass QC
failed_samps <- readRDS('output/failed_QC_samples.rds')

# Bind all batches into data frames
sample_specs <- data.frame()
sample_sheets <- data.frame()
sample_cpgs <- data.frame()
for(batch_no in c(1:3, 9)) {

  # Sample specs
  specs <- readRDS(paste0('input/batch', batch_no, '_samples.rds')) %>%
    select(sampleId, Spec, age, sex) %>%
    dplyr::rename('Sample_Name' = sampleId, 'Age' = age, 'Sex' = sex) %>%
    mutate(Year = substr(Sample_Name, 8, 11)) %>%
    # Filter bad samples
    filter(! Sample_Name %in% failed_samps)

  # Load updated sample sheet after normalizing and join with sample specs
  sheet <- readRDS(paste0('output/updated_sample_sheet_PB_array',
                                 batch_no, '.rds')) %>%
    select(Sample_Name, chip.ID.loc) %>%
    left_join(specs)

  # Load CpGs
  # These data include species characteristics needed for clocks (gestation
  # time, age at reproductive maturity, and maximum lifespan) and CGs
  cpgs <- readRDS(paste0('output/tbetas_PB_array', batch_no, '.rds')) %>%
    rownames_to_column('chip.ID.loc') %>%
    left_join(sheet) %>%
    # Calculate gestation time, maturity, and max lifespan in years using AnAge data
    mutate(GestationTimeInYears = 65/365,
           averagedMaturity.yrs = 1734/365,
           maxAgeCaesar = 43.8) %>%
    #  Move sample characteristics to the beginning of the dataframe for easy location
    relocate(c(Sample_Name, Age, Sex, Year, Spec, GestationTimeInYears:maxAgeCaesar),
             .after = chip.ID.loc)

  # Bind all specs
  sample_specs <- sample_specs %>% bind_rows(specs)
  # Bind sample sheets
  sample_sheets <- sample_sheets %>% bind_rows(sheet)
  # Bind all cpgs
  sample_cpgs <- sample_cpgs %>% bind_rows(cpgs)

}

# 2 Load clock sites and set clock names ====

# S3.1 to S3.3 tables in Supplmentary Data of Lu et al. 2023 from BioRxiv
glmnet.csv <- c('input/mammclock1.csv', 'input/mammclock2.csv', 'input/mammclock3.csv')

# Variable names used to refer to each clock in output
beta.name <- c('beta_clock1', 'beta_clock2', 'beta_clock3')
y.name <- c('Y.pred1', 'Y.pred2', 'Y.pred3') 
age.name <- c('DNAmAgeClock1', 'DNAmAgeClock2', 'DNAmAgeClock3')

# 2 Functions ====

# Clock 2 function (Max age and gestation)
F2_antitrans_clock2 <- function(y, y.maxAge, y.gestation, const = 1) {
  x0 <- const * exp(- exp(- 1 * y))
  x1 <- x0 * (y.maxAge + y.gestation)
  x <- x1 - y.gestation
  x
}

# Clock 3 function 
F1_logli <- function(age1, m1, m2 = m1, c1 = 1) {
  ifelse(age1 >= m1, (age1 - m1) / m2 , c1 * log((age1 - m1) / m2 / c1 + 1)) 
}

#Relative Adult Age
F2_revtrsf_clock3 <- function(y.pred, m1, m2 = m1, c1 = 1) {
  ifelse(y.pred < 0, (exp(y.pred / c1) - 1) * m2 * c1 + m1, y.pred * m2 + m1 ) 
}

# The `loglifn` function shows how to calculate m1 for the transformation 
# It is the `a_Logli` in the function
F3_loglifn <- function(dat1, b1 = 1, max_tage = 4, c1 = 5, c2 = 0.38, c0  = 0) {
  n <- nrow(dat1)
  # Calculate relative age?
  age1 <- (dat1$maxAgeCaesar + dat1$GestationTimeInYears) / 
    (dat1$averagedMaturity.yrs + dat1$GestationTimeInYears)
  a1 <- age1 / (1 + max_tage)
  dat1$a1_Logli <- a1 # x/m1 in manuscript
  a2 = (dat1$GestationTimeInYears + c0) / (dat1$averagedMaturity.yrs) 
  # m=5*(G/ASM)^0.38 from regression analysis/formula(7)
  dat1$a_Logli = a_Logli = c1*a2^c2
  x <- dat1$Age + dat1$GestationTimeInYears
  t2 = dat1$averagedMaturity.yrs * b1 + dat1$GestationTimeInYears 
  x2 = x/t2 # log(x/t2)
  F1_logli <- function(age1, m1, m2 = m1, c1 = 1) {
    ifelse(age1 >= m1, (age1 - m1) / m2 , c1 * log((age1 - m1) / m2 / c1 + 1))
  }
  y = F1_logli(x2, a_Logli, a_Logli)
  dat1$LogliAge <- y
  return(dat1)
}

# 3 Predict ages ====

# Add column for intercept
sample_cpgs$Intercept <- 1

# Factor by which age needs to be multiplied?
MYMAX = 1.3

# Max age x factor
sample_cpgs$HighmaxAgeCaesar = MYMAX * sample_cpgs$maxAgeCaesar 
sample_cpgs$HighmaxAgeCaesar[sample_cpgs$SpeciesLatinName == 'Ursus maritimus'] <- 
  sample_cpgs$maxAgeCaesar[sample_cpgs$SpeciesLatinName == 'Ursus maritimus']

# predict age for clocks 1-3
for(k in 1:3){
  glmnet <- read.csv(glmnet.csv[k])
  glmnet$beta <- glmnet[, beta.name[k]] 
  glmnet$var[1] <- ifelse(glmnet$var[1] == "(Intercept)", 'Intercept', glmnet$var[1]) 
  sample_cpgs[,y.name[k]] = as.numeric(as.matrix(subset(sample_cpgs,
                                                 select=as.character(glmnet$var))) 
                                %*% glmnet$beta)
}

#(1) Clock 1
sample_cpgs[,age.name[1]] <- exp(sample_cpgs[, y.name[k]]) - 2
#(2) Clock 2
sample_cpgs$DNAmRelativeAge <- exp(- exp(- 1 * sample_cpgs[, y.name[2]]))
sample_cpgs[, age.name[2]] <- F2_antitrans_clock2(sample_cpgs[, y.name[2]], 
                                                  sample_cpgs$HighmaxAgeCaesar, 
                                                  sample_cpgs$GestationTimeInYears, const = 1) 

sample_cpgs <- F3_loglifn(sample_cpgs)# to compute m estimate for tuning point in the ll transformation 

# Below m1 is in the original code
sample_cpgs$m1 = sample_cpgs$a_Logli
sample_cpgs$DNAmRelativeAdultAge <- F2_revtrsf_clock3(sample_cpgs[, y.name[3]], sample_cpgs$m1)
sample_cpgs[, age.name[3]]<- sample_cpgs$DNAmRelativeAdultAge * 
  (sample_cpgs$averagedMaturity.yrs + sample_cpgs$GestationTimeInYears) - sample_cpgs$GestationTimeInYears

# 4 Final output ====

# Output
output <- subset(sample_cpgs,select=c('Sample_Name','Age', 'Sex', 'Spec', 'Year', 
                               'DNAmRelativeAge','DNAmRelativeAdultAge',age.name)) %>%
  mutate(Year = as.numeric(Year)) %>%
  # Remove NAs and fix tissue type misspelling
  na.omit() %>%
  mutate(Spec = ifelse(Spec %in% c('Blood', 'Skin'), Spec, 'Blood'),
         # Add acceleration
         AccelClock2 = lm(sample_cpgs$DNAmAgeClock2 ~ sample_cpgs$Age)$residuals,
         AccelClock3 = lm(sample_cpgs$DNAmAgeClock3 ~ sample_cpgs$Age)$residuals,
         Born = Year - Age)

# 6 MAEs and correlations ====

clock_performance <- data.frame()
for(i in 1:3) {
  # Chronological age and predictions from each clock
  age_chron <- subset(output, select = Age) %>% pull()
  age_epi <- subset(output, select = paste0('DNAmAgeClock', i)) %>% pull()
  # Calculate MAE and correlation epi age ~ chronological age
  age_mae <- median(abs(age_epi - age_chron), na.rm = T)
  age_corr <- as.numeric(cor.test(age_epi, age_chron)$estimate)
  # Build df
  clock_results <- data.frame(Clock = i, MAE = age_mae, correlation = age_corr)
  clock_performance <- bind_rows(clock_performance, clock_results)
}

print(clock_performance)

# 5 Plot mammal clocks 2 & 3 ====

uc_clocks <- list()
for(i in 2:3) {
  # Get age predictions for clock
  plot_dat <- subset(output, select = c(Age, Spec, Born,
                                        get(paste0('DNAmAgeClock', i)),
                                        get(paste0('AccelClock', i))))
  colnames(plot_dat) <- c('Age', 'Spec', 'Born', 'AgePredict', 'AgeAccel')
    
  # Plot universal clock by tissue
  clock_plot <- ggplot(data = plot_dat, 
                       aes(x = Age, y = AgePredict, colour = Spec)) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    geom_point(size = 2) +
    geom_smooth(method = 'lm', linewidth = 1.5, colour = 'black', se = F) +
    scale_colour_manual(values = c('#d62d20', '#536878')) +
    ylab('Epigenetic age') + xlab('Chronological age') +
    theme(plot.margin = unit(c(0.5, 0.5, 1, 2), 'cm'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(),
          axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
          axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
          axis.text = element_text(colour = 'black', size = 18),
          legend.key = element_rect(colour = 'white'),
          legend.position = 'none')
  
  # Plot age accel ~ birth year relationship
  accel_plot <- ggplot() +
    geom_point(data = plot_dat, aes(x = Born, y = AgeAccel, colour = Spec), size = 2) +
    geom_smooth(data = plot_dat, aes(x = Born, y = AgeAccel, group = Spec, colour = Spec),
                method = 'lm', se = F, size = 1.5) +
    geom_smooth(data = plot_dat, aes(x = Born, y = AgeAccel),
                colour = 'black', method = 'lm', se = F, size = 1.5) +
    scale_colour_manual(values = c('#d62d20', '#536878')) +
    ylab('Epigenetic age acceleration') + xlab('Year of birth') +
    theme(plot.margin = unit(c(0.5, 0.5, 1, 2), 'cm'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(),
          axis.title.x = element_text(colour = 'black', size = 18, vjust = -5),
          axis.title.y = element_text(colour = 'black', size = 18, vjust = 5),
          axis.text = element_text(colour = 'black', size = 18),
          legend.position = 'none')
  
  # Add panel plot to list
  uc_clocks[[paste0('UC clock panel ', i)]] <- 
    plot_grid(clock_plot, accel_plot, ncol = 2, labels = c('A', 'B'), 
              label_size = 20)
  
  # Add basic clock plot to list
  uc_clocks[[paste0('UC basic clock ', i)]] <- clock_plot
  
}

# Make panel plot with universal clocks 2 & 3 together
uc_clock2 <- uc_clocks$`UC basic clock 2` + 
  labs(x = '', y = '') + 
  theme(plot.margin = unit(c(rep(0.25, 3), 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylim(0, 52.5)

uc_clock3 <- uc_clocks$`UC basic clock 3` + 
  labs(x = '', y = '') + 
  theme(plot.margin = unit(c(rep(0.25, 3), 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  ylim(0, 52.5)

# Plot panels
clock_panel <- plot_grid(uc_clock2, uc_clock3, 
                         ncol = 2, labels = c('A', 'B'), label_size = 22)

# X and Y labels for ewas panel plot
Ylab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Epigenetic age (years)', 
                             size = 7, angle = 90) + theme_void()
Xlab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Chronological age (years)', 
                             size = 7, hjust = 0.4) + theme_void()

# Add axis labels
clock_panel_y <- plot_grid(Ylab, clock_panel, rel_widths = c(0.1, 1))
plot_grid(clock_panel_y, Xlab, rel_heights = c(1, 0.08), ncol = 1)

# 6 Save plots ====
ggsave('uc_combined_plot.tiff', plot = last_plot(), path = 'figures/supplementary/', 
       device = 'tiff', dpi = 300, height = 12, width = 27, units = 'cm', bg = 'white')

# 7 Model age acceleration ~ birth year for both clocks ====

# Get data ready
uc_ages <- output %>%
  mutate(BearID = substr(Sample_Name, 1,6),
         across(c(Born, AccelClock2, AccelClock3), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Clock 2
accel_born_uc2_mod <-  brm(AccelClock2_sc ~ Born_sc + Sex + (Born_sc + Sex | BearID),
                           data = uc_ages, family = gaussian, 
                           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                           prior = prior(normal(0,1), class = b),
                           control = list(adapt_delta = 0.99, max_treedepth = 20),
                           backend = 'cmdstanr')

# Clock 3
accel_born_uc3_mod <-  brm(AccelClock3_sc ~ Born_sc + Sex + (Born_sc + Sex | BearID),
                           data = uc_ages, family = gaussian, 
                           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                           prior = prior(normal(0,1), class = b),
                           control = list(adapt_delta = 0.99, max_treedepth = 20),
                           backend = 'cmdstanr')
