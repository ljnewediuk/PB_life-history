# PB_life-history

Data and code to accompany **Climate change, age acceleration, and the erosion of fitness in polar bears**

Levi Newediuk1, Evan S Richardson1,2, *, Brooke A. Biddlecombe1, Haziqa Kassim3, Leah Kathan1, Nicholas Lunn2, L Ruth Rivkin1,4,5, Ola E Salama3, Chloé Schmidt6, Meaghan J Jones3, *, Colin J Garroway1, *
1 Department of Biological Sciences, University of Manitoba, Canada
2 Environment and Climate Change Canada
3 Department of Biochemistry and Medical Genetics, University of Manitoba, and Children's Hospital Research Institute of Manitoba, Canada
4 Polar Bears International, Bozeman, MT, USA
5 San Diego Zoo Wildlife Alliance, Escondido, California, USA
6 German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig, Leipzig, Germany

**Abstract**
Climate change is increasingly disrupting evolved life history strategies and decreasing population viability in wild species1. The magnitude and pace at which environments will change mean the persistence of wild populations will depend substantially on their ability to adapt genetically. However, we know very little about the capacity for evolutionary change in response to climate warming. We mapped the effects of climate change, beginning with the decline of cellular function through to the erosion of fitness and adaptive potential in an intensively studied polar bear (Ursus maritimus) population in western Hudson Bay, Canada. Using estimates of epigenetic age acceleration, an indicator of declining cellular function associated with exposure to stress, we found that polar bears aged approximately one year faster, on average, for each degree Celsius temperature increase they experienced. Declining cellular function should reduce fitness3,4 and counter adaptive evolution in rapidly changing environments. Individuals who reproduced early had higher lifetime reproductive success; however, this was before the onset of rapid warming. Fitness benefits associated with early reproduction declined with warming, and today, bears have similar lifetime reproductive success regardless of when they first reproduce. Finally, using a large pedigree, we found no evidence for genetic variation associated with reproductive success in this population—the population is not evolving in response to the changing environment. The physiological costs of climate change accumulate across lifetimes to degrade cellular function and, ultimately, adaptive capacity. These findings warn that adaptive responses to warming could be the exception rather than the rule. 

In this paper we use epigenetic aging, a biomedical technique that tracks the deterioration of cellular function over lifetimes, to test the relationship between climate change and stress in polar bears (*Ursus maritimus*)

All code was run in R v4.3.1 and runs on macOS Monterey v12.7.3

Package versions: 

* tidyverse v2.0.0
* sesame v1.18.4
* minfi v1.46.0
* limma v3.56.2
* MCMCglmm v2.35
* ggdist v3.3.1
* cowplot v1.1.1
* bayesplot v1.10.0
* brms v2.20.4
* tidybayes v3.0.6
* glmnet v4.1-8

Two developmental versions of sesame and minfi packages are also required for script 1, normalizing betas. These packages are available in the input files.


Scripts include:

01 - Normalizing betas from methylation data

02 - Epigenome-wide association survey to select initial set of CpGs for building clock that are correlated with age but not with sex

03 - Building the polar bear clock using independent testing and training data with glmnet

04 - Estimating new ages from n = 94 bears in the later batch

05 - Clock performance metrics (median absolute error and correlation), creating a table to summarize testing and training samples for the supplement, and ploting the clock

06 - Summarizing lifetime reproductive success and ages at first reproduction for bears in pedigree

07 - Bayesian GLMs for testing relationships between age acceleration and birth year, age acceleration and age at first reproduction, and lifetime reproductive success and age at first reproduction

08 - Fitted effects from Bayesian GLM models

09 - Model plots in manuscript (Figures 3-5)

10 - Animal models for estimating heritability of lifetime reproductive success

X1 - (Supplement) Code used for volcano plots showing sites correlated with age, tissue, and sex

X2 - (Supplement) Test of polar bear samples aged using universal mammal clock

X3 - (Supplement) Test clock using early, less stressed bears to predict age in bears experiencing stress from climate change

X4 - (Supplement) Test clock using mature bears to predict age in mature bears

X5 - (Supplement) Quality control (detection p-values and PCAs)

X6 - (Supplement) Check for any difference in aging rates over time by tissue or sex

X7 - (Supplement) Resample individuals and sites and re-fit clock to test for overfitting

Folders:

"iscans" folder contains the raw R/G idat files for processing, separated by batch. Samples used for clock development came from batches 1-3, and validation was done with batch 9 and the remaining batch 1-3 samples not used for clock development

Input data includes:

* batch#_samples.rds: Information about samples included in each batch (3x96)
    * sampleId: unique ID including bear ID, date, and sample type
    * id: bear ID
    * Spec: tissue type (blood or skin)
    * YMD: year of collection
    * sex: M = male, F = female
    * age: age of bear at time of sample
    * Born: year of bear birth

* bear_capture_info.csv: Original bear capture Information
    * BearCode: bear ID
    * Date: date of sample
    * tbl_Bears_Sex: sex of bear; M = male, F = female
    * Born: year of birth
    * AQual: how chronological age was estimated; DERIVE = calculated based on known age for bears captured as COYs, TOOTH = age estimated from tooth cementum
    * VisitAge: age of bear at sample
    * NumCubs: number of cubs with bear
    * ABearCodeListX: Accompanying bears with bear ID, sex, and age

* bearPED.csv: Bear pedigree data from Malenfant et al. 2016 (temporarily embargoed)
    * id: bear ID
    * dam: female parent of bear
    * sire: male parent of bear

* breeding.csv: Life history traits for all years of captures
    * BearID: bear ID
    * Born: year of bear birth
    * FirstRepro: age of bear at first known reproduction
    * LastRepro: age of bear at last known reproduction
    * NOffspringYr: average number of offspring produced per year
    * LRS: lifetime reproductive success
    * Sex: F = female; M = male

* simBreeding.csv: Life history info from simulated pedigree to run reproducible example of animal models (script 10). NOTE: SIMULATED DATA. REAL DATA ARE EMBARGOED FOR DATA SHARING RESTRICTIONS.

* simBearPED.csv: Simulated pedigree to run reproducible example of animal models (script 10). NOTE: SIMULATED DATA. REAL DATA ARE EMBARGOED FOR DATA SHARING RESTRICTIONS.

* full_sibs.rds: Character vector of bear IDs that are full siblings or offspring of other bears in the aging set

* list_test_bears.rds: Character vector of bears included in clock testing set

* list_train_bears.rds: Character vector of bears included in clock training set

* PB_array#_sample_sheet#.rds: Technical microarray information about sample positions in 96-well plates and on microarray chips
    * Sample_Name: unique ID including bear ID, date, and sample type
    * Sample_Well: location on plate (row letter and column number)
    * Sample_Plate: plate number (in order run)
    * chip.No: unique chip number (8 per run)
    * chip.ID: unique chip id (8 per run)
    * Stripe: location of sample in rows and columns of chip, each with two columns (C) and 6 rows (R)
    * row: row of 96-well plate
    * column: column of 96-well plate

* HorvathMammal40.CanonicalManifest.3.2019.sesame.csv: Manifest file

* HorvathMammal40.Manifest.May2020.manifest.csv: Manifest file

* mamm_chip_probes_265275085cfa.bam: Mammal chip probes alignment

* mammclock#.csv: Universal mammal clocks (from Lu et al. 2023)

Output data include:

* PB_clock_ages.rds: Original bears (n = 134) aged using polar bear epigenetic clock
   * SampleID: unique sample ID including bear ID, date, and sample type
   * BearID: unique bear ID
   * Age: Age of bear at time of sampling
   * Sex: F = female; M = male
   * Spec: tissue type (blood or skin)
   * AgePredict: predicted epigenetic age based on polar bear clock
   * AgeAccel: Residual from lm(AgePredict ~ Age)
   * yr: year of collection
   
* WH_combined_ages.rds: Original n = 134 bears and n = 94 additional samples processed and aged in the later batch 9
  * Sample_Name: unique sample ID including bear ID, date, and sample type
  * BearID: unique bear ID
  * Spec: tissue type (blood or skin)
  * YMD: Year, month, and day of sample collection
  * Sex: F = female; M = male
  * Born: Birth year of the bear
  * Age: Age of bear at time of sampling
  * AgePredict: predicted epigenetic age based on polar bear clock
  * AgeAccel: Residual from lm(AgePredict ~ Age)
 
* clock_Cpgs.rds: CpG sites included in clock

* detection_p_batch#.rds: Sample detection p-values for quality control

* failed_QC_samples.rds: IDs of samples that failed all quality control tests (these samples are removed from the analysis)

* f_effects_####: Fitted effects from corresponding model

* updated_sample_sheet_PB_array#.rds: Original sample sheets with addition of locations for specific idat files in "iscans" folder for joining
    * chip.ID.loc: unique chip id and location of sample in rows and columns on chip
    * Basename: location of corresponding idat file in iscans folder


* nbetas_PB_array#.rds & tbetas_PB_array#.rds: Matrices of normalized betas and transposed normalized betas from raw idat files
  
* lh_info_epi.rds: Age acceleration data for bears with life history data
   * BearID: unique bear ID
   * Born: year of birth
   * Age: Age of bear at time of sampling
   * Sex: F = female; M = male
   * AgePredict: predicted epigenetic age based on polar bear clock
   * AgeAccel: Residual from lm(AgePredict ~ Age)
   * FirstRepro: age of bear at first known reproduction
   * LastRepro: age of bear at last known reproduction
   * NOffspringYr: average number of offspring produced per year
   * LostLitts: Estimated number of litters lost over lifetime
   * LRS: lifetime reproductive success
   
* lh_info_pop.rds: Life history data for whole population, not including age acceleration data
   * BearID: unique bear ID
   * Born: year of birth
   * Age: Age of bear at time of sampling
   * Sex: F = female; M = male
   * FirstRepro: age of bear at first known reproduction
   * LastRepro: age of bear at last known reproduction
   * NOffspringYr: average number of offspring produced per year
   * LostLitts: Estimated number of litters lost over lifetime
   * LRS: lifetime reproductive success

* mcmc_LRS.rds: MCMCglmm animal model produced using true pedigree data in script 10

* supplementary_bear_data.csv: Information about samples used for training and testing
  * ID: unique bear ID
  * DateSampled: Year, month, and day of sample collection
  * SampleType: tissue type (blood or skin)
  * Age: Age of bear at time of sampling
  * Sex: F = female; M = male
  * Testing: Was the sample used for predicting (yes/no)?
  * Training: Was the sample used to train the clock (yes/no)?

* limma_XXXX.rds: Limma models specifying relationship of each site with variables age, blood/skin, and tissue/sex

* bootstrap_clock_metrics: Median absolute errors and correlations from 500 resampled clocks (see script X7)

* bootstrap_posterior_pred.rds: Posterior predictive distribution combined from all resampled clocks in script X7
