#--------------------#
# Data creation      #
# Kevin Potter       #
# Updated 02/07/2017 #
#--------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Index
# Lookup - 01:  Load in useful packages and change directory
# Lookup - 02:  Read in demographic information
# Lookup - 03:  Combine participant logs
# Lookup - 04:  Extract data from .csv files
# Lookup - 05:  Exclude certain subjects
# Lookup - 06:  Unify subject numbers across data frames
# Lookup - 07:  Save results as .RData file

###
### Load in useful packages and change directory
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Change directory to location of data
setwd( 'Data' )
setwd( 'Wimber_dat' )

###
### Read in demographic information
###
# Lookup - 02

# Get all file names
allFiles = dir()

# Identify demographic output
sel = grep( 'demographics', allFiles )
demoFiles = allFiles[ sel ]

# Determine sample size
N = length( demoFiles )

allDemographics = matrix( NA, N, 5 )
colnames( allDemographics ) = c( 'ID', 'Sex', 'Ethnicity', 
                                 'Race', 'Age' )
allDemographics = as.data.frame( allDemographics )
# Determine duplicate subject numbers
dup = grep( '(2)', demoFiles, fixed = T )

for ( i in 1:length( demoFiles ) ) {
  
  # Extract ID number
  s = as.numeric( strsplit( demoFiles[i], split = '_' )[[1]][2] )
  if ( i %in% dup[1] ) s = s + .2
  allDemographics$ID[ i ] = s
  
  # Extract information
  txt = scan( file = demoFiles[i], what = 'character', sep = '\n' )
  allDemographics$Sex[ i ] = 
    strsplit( txt[1], split = ', ' )[[1]][2]
  allDemographics$Ethnicity[ i ] = 
    strsplit( txt[2], split = ', ' )[[1]][2]
  allDemographics$Race[ i ] = 
    strsplit( txt[3], split = ', ' )[[1]][2]
  allDemographics$Age[ i ] = 
    strsplit( txt[4], split = ', ' )[[1]][2]
  
}

###
### Combine participant logs
###
# Lookup - 03

# Identify logs
sel = grep( 'log', allFiles )
logFiles = allFiles[ sel ]

# Initialize variable to store participant logs
allLogs = c()

for ( i in 1:length( logFiles ) ) {
  
  logf = read.csv( file = logFiles[i] )
  logf[,1] = as.character( logf[,1] )
  
  # There are two subjects with duplicate file names.
  # The first duplicate is subject 12, room 201F
  if ( i == 3 ) {
    sel = which( logf$Output == "Subject_12.csv-Subject_12.mat" )
    logf$Output[sel] = "Subject_12 (2).csv-Subject_12 (2).mat"
  }
  
  allLogs = rbind( allLogs, logf )
  
}

# Extract subject ID number
allLogs$ID = NA
for ( i in 1:nrow( allLogs ) ) {
  
  tmp = strsplit( allLogs$Output[i], split = '.csv' )[[1]]
  if ( tmp[1] == "Subject_12 (2)" ) s = 12.2 else 
    s = strsplit( tmp[1], split = '_' )[[1]][2]
  allLogs$ID[ i ] = as.numeric( s )
}

# Save log as .csv file
setwd( '..' )
write.csv(allLogs,file='Participant_log.csv')
setwd( 'Wimber_dat' )

###
### Extract data from .csv files
###
# Lookup - 04

sel = grep( '.csv', allFiles )
tmp = allFiles[ sel ] # Isolate .csv files
sel = grep( '_Stim', tmp ) # Stimulus files
tmp = tmp[ -sel ]
sel = grep( 'log', tmp ) # Log files
dataFiles = tmp[ -sel ]

# Determine sample size
N = length( dataFiles )

# Initialize variable for data
rawDat = c()

# Loop through files
for ( i in 1:length( dataFiles ) ) {
  
  # Extract ID number
  tmp = strsplit( dataFiles[i], split = '_' )[[1]][2]
  s = as.numeric( strsplit( tmp, split = '.csv' )[[1]] )
  if ( is.na( s ) ) s = 12.2
  
  # Read in data
  dat = read.table( file = dataFiles[i], sep = ',', header = T )
  cat( c( s, nrow(dat), '\n' ) )
  dat$ID = s # Store subject ID number
  
  rawDat = rbind( rawDat, dat )
  
}

###
### Exclude certain subjects
###
# Lookup - 05

# Based on the participant logs...
# Subject  1: Responded too early for first 5 trials of SR phase.
# Subject  2: Computer error, keyboard stopped working during 
#             SR phase.
# Subject  7: Found SR phase to be too fast.
# Subject  9: Found SR phase to be too challenging.
# Subject 22: Confused about purpose of initial familiarization phase
#             (Either whether it was a recognition memory task
#              or indicating sense of familiarity).
# Subject 31: Provided almost no responses to final recognition 
#             memory phase.
# Subject 51: Provided incorrect input during SR phase.

# We choose to exclude subjects with computer errors or poor 
# understanding of the task instructions. Hence, subjects 2, 
# 22, 31, and 51 will be excluded.
# We will assume subjects 1, 7, and 9 are fine.
# Note this exclusion process is a priori, done before any 
# descriptive statistics or analyses were conducted.

allData = rawDat
excludeSubjects = c( 2, 22, 31, 51 )
for ( i in 1:length( excludeSubjects ) ) {
  sel = which( allData$ID == excludeSubjects[i] )
  allData = allData[ -sel, ]
}
allData$Subject = createIncrement( allData$ID )

# Define new sample size
N = max( allData$Subject )

###
### Unify subject numbers across data frames
###
# Lookup - 06

subjKey = aggregate( allData$ID, list( allData$Subject ), unique )
colnames( subjKey ) = c( 'S', 'ID' )

# Demographics information
allDemographics$Subject = NA
for ( i in 1:nrow( subjKey ) ) {
  id = subjKey$ID[ i ]
  sel = which( allDemographics$ID == id )
  allDemographics$Subject[ sel ] = subjKey$S[ i ]
}
allDemographics = allDemographics[ 
  order( allDemographics$Subject ), ]

# Participant logs
allLogs$Subject = NA
for ( i in 1:nrow( subjKey ) ) {
  id = subjKey$ID[ i ]
  sel = which( allLogs$ID == id )
  allLogs$Subject[ sel ] = subjKey$S[ i ]
}
allLogs = allLogs[ 
  order( allLogs$Subject ), ]

# Reorder data
allData = allData[ order( allData$Subject ), ]

###
### Save results as .RData file
###
# Lookup - 07

setwd( '..' )
save( N, rawDat, allData, allDemographics, allLogs, 
      file = 'Wimber_rep_3.RData' )

setwd( orig_dir )