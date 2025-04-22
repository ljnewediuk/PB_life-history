
# 08 - Animal models ====

# Author: Leah Kathan & Levi Newediuk

#===============================================================================
#POLAR BEAR EPIGENETICS
#HERITABILITY OF Lifetime Reproductive Success
#===============================================================================

#------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(MCMCglmm)


#------------------------------------------------------------------------------

#===============================================================================
#PEDIGREE DATA
#===============================================================================
#read in the pedigree file (id, dam, sire)
#bearPED.csv corrected for dams and sires issues (don't use WH_ped.txt)
#NOTE: bearPED.csv is an embargoed file due to restrictions on sharing the 
#population pedigree. We have provided a simulated version of the pedigree 
#(simBearPED.csv) with approximately the same number of individuals as the 
#western Hudson Bay pedigree. Note this is SIMULATED DATA and NOT 
#REPRESENTATIVE OF THE POPULATION, but does provide a reproducible example 
#for running through this script.
bearPED<- read.csv("input/bearPED.csv", header=T)
colnames(bearPED)[1] <- "animal"
bearPED$animal<-as.factor(bearPED$animal)
bearPED$dam<-as.factor(bearPED$dam)
bearPED$sire<-as.factor(bearPED$sire)
str(bearPED)
head(bearPED)
#*all individuals included in phenotypic data file must be included as 
# offspring in the pedigree file*
#------------------------------------------------------------------------------

#===============================================================================
#POLAR BEAR REPRODUCTIVE DATA
#===============================================================================
#------------------------------------------------------------------------------
#load year of first breeding data (BearID, Born, FirstRepro, LastRepro, NOffspringYr, LRS, Sex)
#NOTE: breeding.csv is an embargoed file due to restrictions on sharing the 
#population pedigree. We have provided a simulated version of these data 
#(simBreeding.csv) based on the simulated pedigree. Note this is SIMULATED 
#DATA and NOT REPRESENTATIVE OF THE POPULATION, but does provide a reproducible 
#example for running through this script.
breeding <-read.csv("input/breeding.csv", header=T)

# Add the dam information from bearPED to breeding
breeding$dam <- bearPED$dam[match(breeding$BearID, bearPED$animal)]

#subset data for time frame—birth years up to 1996 because we cannot fully
#capture lifetime reproductive success for bears born after 1996 (reproductive
#senescence is ~ age 20 and the pedigree was published in 2016)
breeding <- breeding[breeding$Born <= 1996, ]

#------------------------------------------------------------------------------
#match IDs in phenotypic data with pedigree data
#create new column to indicate if each BearID is in the pedigree
breeding$in_pedigree <- ifelse(breeding$BearID %in% bearPED$animal, 
                               "in_pedigree","na_pedigree")
table(breeding$in_pedigree) 
#add ID with phenotypic data but aren't in original pedigree and add to bearPED      
new_bearPED_rows <- breeding[breeding$in_pedigree == 
                               "na_pedigree", "BearID"]
unique_new_bearPED_rows <- unique(new_bearPED_rows)

new_rows_data <- data.frame(
  animal = unique_new_bearPED_rows,
  dam = rep(NA, length(unique_new_bearPED_rows)),
  sire = rep(NA, length(unique_new_bearPED_rows))
)
bearPED <- rbind(bearPED, new_rows_data)
#no IDs added to pedigree
#------------------------------------------------------------------------------
breeding_ped<-breeding
breeding_ped$BearID<-as.factor(breeding_ped$BearID)
breeding_ped$animal<-as.factor(breeding_ped$BearID) #duplicate BearID as "animal" for animal model breeding values
breeding_ped$Born<-as.factor(breeding_ped$Born) #birth year
breeding_ped$FirstRepro<-as.numeric(breeding_ped$FirstRepro) #age of first reproduction
breeding_ped$LastRepro<-as.numeric(breeding_ped$LastRepro) #age of last reproduction
breeding_ped$NOffspringYr<-as.numeric(breeding_ped$NOffspringYr) #ratio of number of offspring per year for years alive/not using
breeding_ped$LRS<-as.numeric(breeding_ped$LRS) #lifetime reproductive success
breeding_ped$Sex<-as.factor(breeding_ped$Sex) 
breeding_ped$dam<-as.factor(breeding_ped$dam) 
breeding_ped$in_pedigree<-as.factor(breeding_ped$in_pedigree) #BearID in pedigree
str(breeding_ped)
#BUT REMOVE the 4 ids with parentage issues
breeding_ped <- breeding_ped %>% subset(BearID != "X12780" & 
                                          BearID != "X17406" &
                                          BearID != "X17339" & 
                                          BearID != "X17059")
breeding_ped <- breeding_ped %>% subset(Sex != "NA" & Sex != "U")
breeding_ped <- breeding_ped %>% subset(dam != "NA" & dam != "U")                 
str(breeding_ped)
num_unique_ids <- breeding_ped %>% distinct(BearID) %>% nrow()
#456 individuals
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#OBSERVED DATA SCALE IS POISSON, USE LINK FUNCTION 
#MCMCglmm/GLMMs use 3 data scales: observed, expected, and latent data scales
# because phenotypic traits often cannot be modeled by a normally distributed 
# random error
#latent trait needs to be gaussian (normally distributed genetic 
# component)/varies -/+infinity
#observed data gets transformed with link-function and back with 
# inverse link-function
#error terms need to be normally distributed on latent scale but can follow 
# other distributions on observed data scale
#------------------------------------------------------------------------------
#use a poisson distribution for observed data
LRS.hist <- hist(breeding_ped$LRS,
                 main="Polar bear lifetime reproductive success",
                 xlab="Lifetime reproductive success",
                 ylab="Number of polar bears",
                 xlim=c(0, 20),
                 ylim=c(0,500))
qqnorm(breeding_ped$LRS)
qqline(breeding_ped$LRS)

# Filter the data for Females (F) and Males (M)
female_lrs <- breeding_ped[breeding_ped$Sex == "F", "LRS"]
male_lrs <- breeding_ped[breeding_ped$Sex == "M", "LRS"]

# Create histograms
hist(female_lrs,
     main="Polar bear lifetime reproductive success (Females)",
     xlab="Lifetime reproductive success",
     ylab="Number of polar bears",
     xlim=c(0, 20),
     ylim=c(0, 500),
     col="pink")

hist(male_lrs,
     main="Polar bear lifetime reproductive success (Males)",
     xlab="Lifetime reproductive success",
     ylab="Number of polar bears",
     xlim=c(0, 20),
     ylim=c(0, 500),
     col="blue")

#------------------------------------------------------------------------------

#===============================================================================
#MCMCglmm Animal models
#===============================================================================

#------------------------------------------------------------------------------
#PRIORS
#prior allows model to fit different variance structures
#G(additive genetic variance matrix), R(residual variance matrix)

#V = 1 nu = 1 is relatively uninformative, equivalent to an inverse-gamma prior
#with shape and scale equal to 0.001
#size of genetic and residual variance similar assumed to be similar in prior

#nu: degree of belief parameter (0.002 tells not believe prior expectation)
#V:prior expectation 
#For uninformative priors, nu = 0.002.

#Parameter expanded priors have better behaviour when when the posterior 
#distribution for the variances has support close to zero
#If using parameter expansion, specify prior covariance matrix for alpha

#poisson distribution
#expect values of latent trait and their variance to be small since exponential 

# Weakly informative inverse-Gamma priors
priorIG <- list(R = list(V = 1, nu = 1),
                G = list(G1 = list(V = 1, nu = 1),
                         G2 = list(V = 1, nu = 1),
                         G3 = list(V = 1, nu = 1)))

# Weakly informative inverse-Gamma priors with parameter expansion
priorPE <- list(R = list(V = 1, nu = 1),
                G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

# Very weakly informative near-flat priors with parameter expansion
priorPEW <- list(R = list(V = 1, nu = 0.002),
                G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
bearPED <- bearPED %>% dplyr::select(! Born)
#estimates the inverse relatedness matrix
Ainv <- inverseA(bearPED)$Ainv
#ginv specifies that animal is the random effect linked to the pedigree/inverse 
# relatedness matrix
#model_status has repeated measures for individuals, need to include ID as 
# random effect for repeatability and permanent environmental effects
#acceptance ratio for latent trait (liability set) should be around 0.44
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Model with inverse-Gamma priors
model_LRS.IG <- MCMCglmm(LRS ~ Sex,
                      random = ~animal + dam + Born,
                      ginv = list(animal = Ainv), 
                      family="poisson", data = breeding_ped, 
                      nitt = 65000, thin = 50, 
                      burnin = 15000, prior = priorIG)

# Model with parameter-expanded inverse-Gamma priors
model_LRS.PE <- MCMCglmm(LRS ~ Sex,
                      random = ~animal + dam + Born,
                      ginv = list(animal = Ainv), 
                      family="poisson", data = breeding_ped, 
                      nitt = 65000, thin = 50, 
                      burnin = 15000, prior = priorPE)

# Model with parameter-expanded very weak/flat priors
model_LRS.PEW <- MCMCglmm(LRS ~ Sex,
                         random = ~animal + dam + Born,
                         ginv = list(animal = Ainv), 
                         family="poisson", data = breeding_ped, 
                         nitt = 65000, thin = 50, 
                         burnin = 15000, prior = priorPEW)

# Save parameter-expanded model with flat priors
saveRDS(model_LRS.PEW, file = "mcmc_LRS.rds")
#------------------------------------------------------------------------------

#---------------------------CHECK MODEL PLOTS----------------------------------
model_LRS <- readRDS("output/mcmc_LRS.rds")

#model summary-----------------------------------------------------------------
summary(model_LRS)
#want consistent variation around largely unchanging mean value of intercept

#see if fixed effect (sex) has effect on trait 
#does posterior distribution for sex overlap zero? if yes, no effect
posterior.mode(model_LRS$Sol[, "SexM"])
HPDinterval(model_LRS$Sol[, "SexM"], 0.95)
posterior.mode(model_LRS$VCV)
#does not overlap zero here so sex has an effect on the trait

#posterior distribution of intercept/fixed effects, look for autocorrelation----
plot(model_LRS$Sol)

#posterior distributions of variance components, look for autocorrelation-------
plot(model_LRS$VCV)

#effective sample size from draws/need above 1000-------------------------------
effectiveSize(model_LRS[["Sol"]])
effectiveSize(model_LRS[["VCV"]])

#check for autocorrelation among samples for posterior distributions, ----------
#  should be close to zero except for lag 0
autocorr.diag(model_LRS$Sol)
autocorr.diag(model_LRS$VCV)

#diagnostic tests of convergence (Heidelberg stationary test)-------------------
heidel.diag(model_LRS$VCV)
#-------------------------------------------------------------------------------


# ---------------------------MODEL ESTIMATES ----------------------------------
#number of individuals in data------------------------------------------------
num_unique_individuals <- breeding_ped %>% 
  distinct(BearID) %>% 
  nrow()
cat("Number of unique individuals:", num_unique_individuals, "\n")

#number of observations in data-----------------------------------------------
num_observations <- nrow(breeding_ped)
cat("Number of observations:", num_observations, "\n")

# estimate phenotypic mean and sd---------------------------------------------
model_LRS <- readRDS(file="mcmc_LRS.rds")

mean_LRS <- mean(breeding_ped$LRS)
sd_LRS <- sd(breeding_ped$LRS)
n <- length(breeding_ped$LRS)
# Estimate the standard error
se_LRS <- sd_LRS / sqrt(n)
# Print the mean and standard error
cat("Mean LRS:", mean_LRS, "\n")
cat("Standard Error:", se_LRS, "\n")

# Calculate the sample variance of the Captured binary data
sample_variance_LRS <- var(breeding_ped$LRS)
# Print the sample variance
cat("Sample Variance of LRS Data (VObs):", sample_variance_LRS, "\n")

#No fixed effects
#heritability on liability scale-----------------------------------------------
model_LRS <- readRDS(file="mcmc_LRS.rds")
model_LRS <- model_LRS

# estimate means for random effects--------------------------------------------
#JUST CHANGE which variance looking at in numerator (1-4)
#1=animal/2=dam/3=cohort/4=residuals
posterior.1<-model_LRS.PEW$VCV[,1]
VX <- mean(posterior.1)   
HPDinterval(posterior.1) ##highest posterior density/ credible interval
standard_deviation <- sd(posterior.1) 
standard_error <- standard_deviation / sqrt(length(posterior.1))
#-------------------------------------------------------------------------------
