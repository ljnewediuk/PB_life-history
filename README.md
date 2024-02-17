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

Versions: 

* R v4.3.1
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

Two developmental versions of sesame and minfi packages are also required for script 1. These packages are available in the input files.


Scripts include:

01 - Normalizing betas from methylation data

02 - Epigenome-wide association survey to select initial set of CpGs for building clock that are correlated with age but not with sex

03 - Building the polar bear clock using independent testing and training data with glmnet

04 - Summarizing lifetime reproductive success and ages at first reproduction for bears in pedigree

05 - Bayesian GLMs for testing relationships between age acceleration and birth year, age acceleration and age at first reproduction, and lifetime reproductive success and age at first reproduction

06 - Fitted effects from Bayesian GLM models

07 - Plots in manuscript

08 - Animal models for estimating heritability of lifetime reproductive success

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

* bearPED.csv: Bear pedigree data from Malenfant et al. 2016
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

* updated_sample_sheet_PB_array#.rds: Original sample sheets with addition of locations for specific idat files in "iscans" folder
    * chip.ID.loc: unique chip id and location of sample in rows and columns on chip
    * Basename: location of corresponding idat file in iscans folder


* nbetas_PB_array#.rds & tbetas_PB_array#.rds: Matrices of normalized betas and transposed normalized betas from raw idat files
