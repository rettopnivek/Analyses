#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 01/03/17            #
#-----------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

###
### Load in useful packages
###

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

###
### Read in .csv file
###
# Lookup - 01

# Change directory to where data is stored
setwd( 'Data' )

# Load in .csv file with data
rawDat = read.table( file = 'FlankerExp4.csv', sep = ',', header = T )
# Rename variables
colnames( rawDat ) = c(
  'Subjects', 'Condition', 'Correct', 
  'Accuracy', 'Resp', 'RT', 'Target', 'Preview1', 'Preview3' )

# Convert responses to dummy coding
rawDat$Correct[ rawDat$Correct == 5 ] = 0
rawDat$Resp[ rawDat$Resp == 5 ] = 0

# Convert RT from ms to seconds
rawDat$RT = rawDat$RT/1000

# Subject:   Indicates which subject completed a given trial
# Condition: Indicate the type of condition, where...
#            #  | Prime (ms) | Prime (type)
#            1  | 0          | None
#            2  | 0          | None
#            3  | 0          | None
#            4  | 0          | None
#            5  | 100        | Identical
#            6  | 100        | Same response
#            7  | 100        | Incompatible
#            8  | 100        | Identical (flankers)
#            9  | 100        | Same response (flankers)
#            10 | 100        | Incompatible (flankers)
#            11 | 800        | Identical
#            12 | 800        | Same response
#            13 | 800        | Incompatible
#            14 | 800        | Identical (flankers)
#            15 | 800        | Same response (flankers)
#            16 | 800        | Incompatible (flankers)
# Correct:   Indicates what the correct response is, where 
#            1 = consonant, 0 = vowel
# Accuracy:  Indicates whether a response was correct (1) or 
#            incorrect (0)
# Resp:      Indicates what response the subject made, where 
#            1 = consonant, 0 = vowel
# RT:        The response time (seconds)
# Target:    The onscreen target presentation (following 
#            a 1900 ms fixation and prime presentation).
# Preveiw1:  The fixation presentation.
# Preview3:  The prime presentation.
