#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 12/15/16            #
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

# Index
# Lookup - 01:  Read in .csv files
# Lookup - 02:  Add in additional variables

###
### Read in .csv files
###
# Lookup - 01

# Change directory to where data is stored
setwd( 'Data' )

# Change directory to where original .csv files are stored
setwd('Raw_data')

# Create empty variable to store subsequent data
rawDat = c()

# Extract filenames
allFiles = dir()
# Isolate subject files
sel = grep( 'Subject', allFiles )
datFiles = allFiles[ sel ]
# Isolate .csv files
sel = grep( '.csv', datFiles )
datFiles = datFiles[ sel ]

# Extract sample size
N = length( datFiles )

# Loop over subjects
for (n in 1:N) {
  
  tmp = read.table( file = datFiles[n], 
                    sep = ',', header = T )
  sel = which( apply( tmp, 1, function(x) all( x[-1] == 0 ) ) )
  tmp = tmp[ -sel, ] # Trim out an empty row of zeros
  rawDat = rbind( rawDat, tmp )
  
}

# Clean up workspace
rm( allFiles, sel, datFiles, n, tmp )

###
### Add in additional variables
###
# Lookup - 02

# Define a dummy-coded variable indicating which angle 
# was primed (1 = left, 2 = right)
rawDat$PrimeSide = 0
sel = ( rawDat$PrimeType == 1 & rawDat$Target == 0 ) | 
  ( rawDat$PrimeType == 2 & rawDat$Target == 1 )
rawDat$PrimeSide[ sel ] = 2
sel = ( rawDat$PrimeType == 1 & rawDat$Target == 1 ) | 
  ( rawDat$PrimeType == 2 & rawDat$Target == 0 )
rawDat$PrimeSide[ sel ] = 1

# Define an index for the different conditions of interest
rawDat$Condition = 0
sel = rawDat$BlockType == 2
tmp = aggregate( rep(1,sum(sel)), list( 
   rawDat$Target[sel], rawDat$PrimeDuration[sel], 
   rawDat$PrimeType[sel] ),
  sum )
colnames( tmp ) = c('Co','D','PT','N')
for (i in 1:12 ) {
  sel = rawDat$Target == tmp$Co[i] & 
    rawDat$PrimeDuration == tmp$D[i] & 
    rawDat$PrimeType == tmp$PT[i] & 
    rawDat$BlockType == 2
  rawDat$Condition[sel] = i
}

# Define an index for conditions collapsed over choice
rawDat$DurxPri = 0
sel = rawDat$BlockType == 2
tmp = aggregate( rep(1,sum(sel)), list( 
  rawDat$PrimeDuration[sel], 
  rawDat$PrimeType[sel] ),
  sum )
colnames( tmp ) = c('D','PT','N')
for (i in 1:6 ) {
  sel = rawDat$PrimeDuration == tmp$D[i] & 
    rawDat$PrimeType == tmp$PT[i] & 
    rawDat$BlockType == 2
  rawDat$DurxPri[sel] = i
}

# Shift RTs by .5 seconds, so that they represent the 
# duration from the actual target onset, instead of 
# when the duration from when the mask was removed
rawDat$RT = rawDat$RT + .5

# Create an index for the current block number
rawDat$BlockNumber = 0
for ( n in 1:N ) {
  sel = rawDat$Subject == n & rawDat$BlockType == 2
  rawDat$BlockNumber[ sel ] = rep( 1:6, each = 60 )
}

# Clean up workspace
rm(sel,tmp,i,n)

# Save data set
setwd('..')
save(rawDat,N,file='Gabor_pilot_Dec.RData')

# Columns for rawDat
# 1:  Subject index
# 2:  Response times (in seconds)
# 3:  Choice (0 = left, 1 = right)
# 4:  Accuracy (0 = incorrect, 1 = correct)
# 5;  Angle of target stripes (0 = left, 1 = right)
# 6:  Type of prime (0 = no prime, 1 = foil primed, 
#     2 = target primed)
# 7:  Duration of prime (in ms)
# 8:  Contrast for the foil (Michelson contrast)
# 9:  Contrast for the target (Michelson contrast)
# 10: Illuminance for maximum peak of grey patches
# 11: Angle of stripes
# 12: Type of block (-1 = guided practice, 0 = practice,
#     1 = calibration, 2 = main study)
# 13: Indicates the angle that was primed (0 = no prime, 
#     1 = left primed, 2 = right primed)
# 14: An index of the conditions of interest, where...
#     L = left correct, R = right correct
#     SD = 50 ms prime, LD = 400 ms prime
#     N = neither primed, F = foil primed, T = target primed
#     1  -> L-SD-N
#     2  -> R-SD-N
#     3  -> L-LD-N
#     4  -> R-LD-N
#     5  -> L-SD-F
#     6  -> R-SD-F
#     7  -> L-LD-F
#     8  -> R-LD-F
#     9  -> L-SD-T
#     10 -> R-SD-T
#     11 -> L-LD-T
#     12 -> R-LD-T
# 15: An index of conditions collapsing over choice, where...
#     1  -> SD-N
#     2  -> SD-N
#     3  -> LD-N
#     4  -> LD-N
#     5  -> SD-F
#     6  -> SD-F
# 16: An index of the current block for the main study trials

