#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 09/28/2017          #
#-----------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicator for whether or not to save .RData files
saveFile = TRUE

# Index
# Lookup - 01:  Useful packages and functions
# Lookup - 02:  Pilot data
# Lookup - 03:  Current study
# Lookup - 04:  Data key

###
### Useful packages and functions
###
# Lookup - 01

### Load in useful packages

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

### Define useful functions

extract_data = function() {
  # Purpose: 
  # Extracts data from .csv files and converts to 
  # data frame and integer for sample size.
  # Returns: 
  # A data frame with the raw data and an integer 
  # with the number of subjects.
  
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
  
  return( list( N = N, rawDat = rawDat ) )
}

demographcs_extract = function() {
  # Purpose: 
  # Extracts the demographic information from the 
  # relevant text files and creates a data frame 
  # with the results.
  # Returns: 
  # A data frame with the demographics information.
  
  # Extract file names
  all_files = dir()
  
  # Determine which files start with "Demographics"
  sel = grep( "Demographics", all_files )
  
  # Isolate files for demographics info
  demo_files = all_files[ sel ]
  
  # Determine sample size
  N = length( demo_files )
  
  # Initialize character matrix
  out = cbind( 
    ID = rep( "blank", N ), 
    Sex = rep( "blank", N ),
    Ethnicity = rep( "blank", N ),
    Race = rep( "blank", N ) )
  
  for ( n in 1:N ) {
    
    # Subject ID
    tmp = strsplit( demo_files[n], '_' )[[1]][2]
    tmp = as.numeric( strsplit( tmp, '.txt' )[[1]] )
    out[n,1] = tmp
    info = scan( file = demo_files[n], what = "character",
                 sep = "\n")
    out[n,2] = strsplit( info[1], split = ", " )[[1]][2]
    out[n,3] = strsplit( info[2], split = ", " )[[1]][2]
    out[n,4] = strsplit( info[3], split = ", " )[[1]][2]
  }
  
  # Convert to data-frame
  out = as.data.frame( out )
  
  return( out )
}

add_variables = function( rawDat ) {
  # Purpose: 
  # Adds additional variables to the data frame of raw data.
  # Arguments: 
  # rawDat - The data frame of raw data generated from extract_data
  # Returns: 
  # A new data frame with additional variables.
  
  # Meaningful label for prime type
  rawDat$PrimeTypeLab = 'Neither primed'
  rawDat$PrimeTypeLab[ rawDat$PrimeType == 1 ] = 'Target primed'
  rawDat$PrimeTypeLab[ rawDat$PrimeType == 2 ] = 'Foil primed'
  
  # Meaningful label for target
  rawDat$TargetLab = 'Left'
  rawDat$TargetLab[ rawDat$Target == 1 ] = 'Right'
  
  # Meaningful label for choice
  rawDat$ChoiceLab = 'Left'
  rawDat$ChoiceLab[ rawDat$Choice == 1 ] = 'Right'
  
  # Define a dummy-coded variable indicating which angle 
  # was primed (0 = left, 1 = right, 2 = neither)
  rawDat$PrimeSide = 2; rawDat$PrimeSideLab = 'Neither'
  sel = ( rawDat$PrimeType == 1 & rawDat$TargetLab == 'Left' ) | 
    ( rawDat$PrimeType == 2 & rawDat$TargetLab == 'Right' )
  rawDat$PrimeSide[ sel ] = 0; rawDat$PrimeSideLab[sel] = 'Left'
  sel = ( rawDat$PrimeType == 1 & rawDat$TargetLab == 'Right' ) | 
    ( rawDat$PrimeType == 2 & rawDat$TargetLab == 'Left' )
  rawDat$PrimeSide[ sel ] = 1; rawDat$PrimeSideLab[ sel ] = 'Right'
  
  # Create an index for the current block number
  rawDat$BlockNumber = -2
  rawDat$BlockNumber[ rawDat$BlockType == 0 ] = -1
  rawDat$BlockNumber[ rawDat$BlockType == 1 ] = 0
  for ( n in 1:N ) {
    sel = rawDat$Subject == n & rawDat$BlockType == 2
    rawDat$BlockNumber[ sel ] = rep( 1:6, each = 60 )
  }
  
  # Define an index for the different conditions of interest
  rawDat$Condition = 0
  sel = rawDat$BlockType == 2
  tmp = aggregate( rep(1,sum(sel)), list( 
    rawDat$Target[sel], rawDat$PrimeDuration[sel], 
    rawDat$PrimeType[sel] ),
    sum )
  colnames( tmp ) = c('Co','D','PT','N')
  for (i in 1:nrow( tmp ) ) {
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
  for (i in 1:nrow( tmp ) ) {
    sel = rawDat$PrimeDuration == tmp$D[i] & 
      rawDat$PrimeType == tmp$PT[i] & 
      rawDat$BlockType == 2
    rawDat$DurxPri[sel] = i
  }
  
  return( rawDat )
}

###
### Pilot data
###
# Lookup - 02

# Relocate to location of raw files
setwd( 'Data/Raw_files/Pilot' )

# Read in data
tmp = extract_data()
N = tmp$N
rawDat = tmp$rawDat
rm( tmp )
demographics = demographcs_extract()
experiment_log = read.csv( 'Experiment_log.csv' )

# First subject completed wrong type of task
# (i.e., Psychophysics task), so remove
N = N - 1
rawDat = rawDat[ rawDat$Subject != 1, ]
demographics = demographics[ -1, ]
experiment_log = experiment_log[-1,]

# NOTE: Subject 4 is missing from the experiment log

# Add in additional variables

# Save original ID
rawDat$ID = rawDat$Subject
rawDat$Subject = createIncrement( rawDat$ID )
demographics$Subject = 1:N

rawDat = add_variables( rawDat )

setwd( orig_dir )

# Save data set
setwd('Data')
if (saveFile) save (rawDat, N, demographics, experiment_log, 
                    file = 'Gabor_pilot_9-17.RData' )
setwd( orig_dir )

###
### Current study
###
# Lookup - 03

# Relocate to location of raw files
setwd( 'Data/Raw_files/Current_study' )

# Read in data
tmp = extract_data()
N = tmp$N
rawDat = tmp$rawDat
rm( tmp )
demographics = demographcs_extract()
experiment_log = read.csv( 'Experiment_log.csv' )

# Add in additional variables

# Save original ID
rawDat$ID = rawDat$Subject
rawDat$Subject = createIncrement( rawDat$ID )
demographics$Subject = 1:N

rawDat = add_variables( rawDat )

setwd( orig_dir )

# Save data set
setwd('Data')
if (saveFile) save (rawDat, N, demographics, experiment_log, 
                    file = 'Gabor_current_9-17.RData' )
setwd( orig_dir )

###
### Data key
###
# Lookup - 04

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
# 13: Original identifying number for subjects
# 14: Meaningful label for prime type
# 15: Meaningful label for target angle
# 16: Meaningful label for choice
# 17: Indicates the angle that was primed (2 = no prime, 
#     0 = left primed, 1 = right primed)
# 18: Meaningful label for side that was primed
# 19: Block number
# 20: An index of the conditions of interest, where...
#     L = left correct, R = right correct
#     N = neither primed, F = foil primed, T = target primed
#     50, 100, 500, 1000 = Prime durations (in ms)
#     1   -> L - 50   - N
#     2   -> R - 50   - N
#     3   -> L - 100  - N
#     4   -> R - 100  - N
#     5   -> L - 500  - N
#     6   -> R - 500  - N
#     7   -> L - 1000 - N
#     8   -> R  -1000-  N
#     9   -> L - 50   - T
#     10  -> R - 50   - T
#     11  -> L - 100  - T
#     12  -> R - 100  - T
#     13  -> L - 500  - T
#     14  -> R - 500  - T
#     15  -> L - 1000 - T
#     16  -> R  -1000-  T
#     17  -> L - 50   - F
#     18  -> R - 50   - F
#     19  -> L - 100  - F
#     20  -> R - 100  - F
#     21  -> L - 500  - F
#     22  -> R - 500  - F
#     23  -> L - 1000 - F
#     24  -> R  -1000-  F
# 21: An index of conditions collapsing over choice, where...
#     1   -> 50   - N
#     2   -> 100  - N
#     3   -> 500  - N
#     4   -> 1000 - N
#     5   -> 50   - T
#     6   -> 100  - T
#     7   -> 500  - T
#     8   -> 1000 - T
#     9   -> 50   - F
#     10  -> 100  - F
#     11  -> 500  - F
#     12  -> 1000 - F
