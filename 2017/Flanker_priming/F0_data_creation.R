#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 04/30/17            #
#-----------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Index
# Lookup - 01:  Load in useful packages
# Lookup - 02:  Read in .csv file
# Lookup - 03:  Create/adjust variables
# Lookup - 04:  Variable key
# Lookup - 05:  Save new data

###
### Load in useful packages
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

###
### Read in .csv file
###
# Lookup - 02

# Change directory to where data is stored
setwd( 'Data' )

# Load in .csv file with data
rawDat = read.table( file = 'FlankerExp4.csv', sep = ',', header = T )
# Rename variables
colnames( rawDat ) = c(
  'Subjects', 'Condition', 'Correct', 
  'Accuracy', 'Resp', 'RT', 'Target', 'Preview1', 'Preview3' )

###
### Create/adjust variables
###
# Lookup - 03

# Convert responses to dummy coding
rawDat$Correct[ rawDat$Correct == 5 ] = 0
rawDat$Resp[ rawDat$Resp == 5 ] = 0

# Convert RT from ms to seconds
rawDat$RT = rawDat$RT/1000

# Create variable storing original subject ID numbers
rawDat$ID = rawDat$Subjects
# Set subject numbers to be incremental
rawDat$Subjects = createIncrement( rawDat$Subjects )

# Create variable for flanker type
# 1 = No flankers
# 2 = Identical (Same letter as target)
# 3 = Same response (Same class as target, i.e. 
#     matching vowel or consonant)
# 4 = Incompatible (consonant/vowel instead of target's vowel/consonant)

rawDat$FlankerLabel = 'None'
rawDat$Flanker = 1

sel = rawDat$Condition == 2 | rawDat$Condition == 8 | 
  rawDat$Condition == 14
rawDat$FlankerLabel[ sel ] = 'Identical'
rawDat$Flanker[ sel ] = 2

sel = rawDat$Condition == 3 | rawDat$Condition == 9 | 
  rawDat$Condition == 15
rawDat$FlankerLabel[ sel ] = 'Same'
rawDat$Flanker[ sel ] = 3

sel = rawDat$Condition == 4 | rawDat$Condition == 10 | 
  rawDat$Condition == 16
rawDat$FlankerLabel[ sel ] = 'Incompatible'
rawDat$Flanker[ sel ] = 4

# Create variable for prime type
# 1 = No prime
# 2 = Identical (Same letter as target)
# 3 = Same response (Same class as target, i.e. 
#     matching vowel or consonant)
# 4 = Incompatible (consonant/vowel instead of target's vowel/consonant)
rawDat$PrimeType = 1
rawDat$PrimeTypeLabel = 'None'
inc = 0
tmp = c( 'Identical', 'Same', 'Incompatible' )
for (i in 2:4) {
  sel = rawDat$Condition == (5 + inc) | 
    rawDat$Condition == (8 + inc) | 
    rawDat$Condition == (11 + inc) | 
    rawDat$Condition == (14 + inc)
  rawDat$PrimeType[ sel ] = i
  rawDat$PrimeTypeLabel[ sel ] = tmp[ i - 1 ]
  inc = inc + 1
}

# Create variable for prime duration (in s)
rawDat$PrimeDuration = 0
rawDat$PrimeDuration[ rawDat$Condition > 4 & 
                        rawDat$Condition < 11 ] = .1
rawDat$PrimeDuration[ rawDat$Condition > 10 ] = .8

# Create meaningful label for correct and choice options
rawDat$CorrectLabel = 'consonant'
rawDat$CorrectLabel[ rawDat$Correct == 0 ] = 'vowel'
rawDat$RespLabel = 'consonant'
rawDat$RespLabel[ rawDat$Choice == 0 ] = 'vowel'

# Extract number of subjects
N = length( unique( rawDat$Subjects ) )

###
### Variable key
###
# Lookup - 04

# Subject:         Indicates which subject completed a given trial
# Condition:       Indicate the type of condition, where...
# N = None, Id = Identical, In = Incompatible, SR = Same response
#                  #  | Prime (ms) | Prime (type) | Flanker
#                  1  | 0          | N            | N
#                  2  | 0          | N            | Id
#                  3  | 0          | N            | SR
#                  4  | 0          | N            | In
#                  5  | 100        | Id           | N
#                  6  | 100        | SR           | N
#                  7  | 100        | In           | N
#                  8  | 100        | Id   Same -> | Id
#                  9  | 100        | SR   Same -> | SR
#                  10 | 100        | In   Same -> | In
#                  11 | 800        | Id           | N
#                  12 | 800        | SR           | N
#                  13 | 800        | In           | N
#                  14 | 800        | Id   Same -> | Id
#                  15 | 800        | SR   Same -> | SR
#                  16 | 800        | In   Same -> | In
# Correct:         Indicates what the correct response is, where 
#                  1 = consonant, 0 = vowel
# Accuracy:        Indicates whether a response was correct (1) or 
#                  incorrect (0)
# Resp:            Indicates what response the subject made, where 
#                  1 = consonant, 0 = vowel
# RT:              The response time (seconds)
# Target:          The onscreen target presentation (following 
#                  a 1900 ms fixation and prime presentation)
# Preveiw1:        The fixation presentation
# Preview3:        The prime presentation
# FlankerLabel:    Interpretable labels for the flanker conditions
# Flanker:         The type of flanker, where...
#                  1 = No flankers
#                  2 = Identical (Same letter as target)
#                  3 = Same response (Same class as target, i.e. 
#                      matching vowel or consonant)
#                  4 = Incompatible (consonant/vowel instead of 
#                      target's vowel/consonant)
# PrimeType:       The type of prime, where...
#                  1 = No prime
#                  2 = Identical (Same letter as target)
#                  3 = Same response (Same class as target, i.e. 
#                      matching vowel or consonant)
#                  4 = Incompatible (consonant/vowel instead of 
#                      target's vowel/consonant)
# PrimeTypeLabel:  Interpretable labels for the prime conditions
# PrimeDuration:   The length of the prime in seconds
# CorrectLabel:    Interpretable labels for the correct response
# RespLabel:       Interpretable labels for the observed response

###
### Save new data
###
# Lookup - 05

save( rawDat, N, file = 'Flanker_priming.RData' )
setwd( orig_dir )
