
library(ChIPseeker)
library(tidyverse)
library(Rsamtools)
library(data.table)

# 02 - Align probes

# Check alignment of polar bear probes with manifest after qAlign
# Modified from code by A. Shafer. The alignment helped us narrow down the
# initial set of probes by excluding those that bind multiple times or don't 
# bind at all in the polar bear genome.

# 1 Load files ====

aln <- BamFile('input/mamm_chip_probes_265275085cfa.bam')
aln <- scanBam(aln)
aln <- as.data.frame(aln[[1]])

# Load manifest
# Add a new column “targetCG” based on values in StrandTB. 
# If StrandTB is B targetCG = 49:50 and if T = 1:2
manifest <- read.csv('input/HorvathMammal40.Manifest.May2020.manifest.csv') %>%
  mutate(targetCG = ifelse(StrandTB == 'B', '49:50', '1:2'))

# 2 Organize ====

# Make sure targetCG is in df and is according to StrandTB
head(manifest, n=20)

# Select only columns needed in manifest, rename IlmnID column to qname, join 
# the manifest df to the aln df using qname, then change targetCG to character
aln_1 <- manifest %>% 
  dplyr::select(IlmnID, SourceSeq, targetCG) %>% 
  dplyr::rename(qname = IlmnID) %>% 
  right_join(aln) %>% 
  mutate(targetCG = as.character(targetCG))

# 3 Create Cgcount ====

CGcount <- rbindlist(lapply(1:nrow(aln_1), function(i){
  pattern <- DNAString(as.character(aln_1$SourceSeq[i]))
  subject <- DNAString(aln_1$seq[i])
  # Find hits between genome and manifest
  matches <- matchPattern(pattern, subject, max.mismatch = 0, algorithm = "naive-inexact")
  locations <- paste(start(matches), end(matches), sep=":")
  pattern2 <-reverseComplement(DNAString(as.character(aln_1$SourceSeq[i])))
  matches2 <- matchPattern(pattern2, subject, max.mismatch = 0, algorithm = "naive-inexact")
  locations2 <- paste(start(matches2), end(matches2), sep=":")
  hits <- data.frame(qname=aln_1$qname[i],
                     CGcount = length(start(matches))+length(start(matches2)), 
                     forward = paste(locations, collapse = " ; "),
                     reverse = paste(locations2, collapse = " ; "))
}))

aln_1$alignedStand <- ifelse(CGcount$forward!="", "forward", "complementReverse")
aln_1$targetCG <- ifelse(aln_1$alignedStand=="forward", aln_1$targetCG, 
                         ifelse(aln_1$alignedStand=="complementReverse"&aln_1$targetCG=="1:2", "49:50",
                                ifelse(aln_1$alignedStand=="complementReverse"&aln_1$targetCG=="49:50", "1:2",NA)))
aln_1$targetCG <- as.numeric(as.character(factor(aln_1$targetCG, levels = c("1:2", "49:50"), labels = c(0,48))))
aln_2 <- aln_1 %>% filter(!is.na(pos))

# Save mapped probes
saveRDS(aln_2, 'output/alignment.rds')


