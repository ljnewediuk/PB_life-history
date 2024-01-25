
# Load QuasR
library('QuasR')

# Sample and genome files
SampleFile1 <- 'SampleFile.txt'
genomeFile <- 'PB_genome.fna'

# Run alignment
alignment <- qAlign(SampleFile1, genomeFile, bisulfite = 'undir')
