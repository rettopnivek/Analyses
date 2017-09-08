#--------------------------#
# SDT parameter estimation #
# Kevin Potter             #
# Updated 02/28/2017       #
#--------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicate which code segments to run
runCode = c( T, T )

# Index
# Lookup - 01:  Load in useful packages, functions, and data
# Lookup - 02:  Estimates from data (13 conditions)
# Lookup - 03:  Estimates from data (14 conditions)

# Glossary
# URL = Unrelated lures
# RL = Related lures
# CL = Critical lures
# T = Target
# P2 = Phonetic (2 item list)
# P8 = Phonetic (2 item list)
# S2 = Semantic (8 item list)
# S8 = Semantic (8 item list)

###
### Load in useful packages, functions, and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Load in useful functions
source( 'F1_Useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'DRMdata.RData' )
setwd( orig_dir )

###
### Estimates from data (13 conditions)
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Initialize a matrix for the estimates
  parEst = matrix( NA, N, 14 )
  
  # Initialize matrix for data
  curDat = matrix( NA, 13, 3 )
  colnames( curDat ) = c('N','Y','P')
  curDat = as.data.frame( curDat )
  
  # Loop over subjects and calculate parameters
  for ( s in 1:N ) {
    
    # Convert 14 conditions into 13 conditions
    tmp = ratesDat[ ratesDat$S == s, ]
    curDat[2:13,] = tmp[3:14,c('N','Y','P')]
    curDat[1,] = c( sum( tmp$N[1:2] ),
                    sum( tmp$Y[1:2] ),
                    sum( tmp$Y[1:2] )/sum( tmp$N[1:2] ) )
    
    parEst[ s, ] = SDT_calc_4_dist( curDat,
                                     correct = T )
  }
  
  # Remove the column for unrelated lures, since the 
  # d' is fixed to 0
  parEst = parEst[,-1]
  
  # Label columns
  colnames( parEst ) = c( as.character( 
    ratesDat$label13[ ratesDat$S == 1 ][-(1:2)] ), 'k' )
  
  # Add in subject ID numbers and age-group
  tmp = aggregate( allDat[,c('ID','Age')], 
                   list( allDat$subject ), unique )
  colnames( tmp ) = c( 'S','ID','Age' )
  parEst = cbind( tmp, parEst )
  
  # Save results as a .csv
  setwd( 'Data' )
  write.table( parEst, file = 'SDT_estimates_13.csv', 
               sep = ',', row.names = F, col.names = T )
  setwd( orig_dir )
  
}

###
### Estimates from data (14 conditions)
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Initialize a matrix for the estimates
  parEst = matrix( NA, N, 16 )
  
  # Loop over subjects and calculate parameters
  for ( s in 1:N ) {
    
    curDat = ratesDat[ ratesDat$S == s, ]
    parEst[ s, ] = SDT_calc_4_dist( curDat,
                                    correct = T )
  }
  
  # Remove the columns for unrelated lures, since the 
  # d' values are fixed to 0
  parEst = parEst[,-(1:2)]
  
  # Label columns
  colnames( parEst ) = c( as.character( 
    ratesDat$label14[ ratesDat$S == 1 ][-(1:2)] ), 'kp', 'ks' )
  
  # Add in subject ID numbers and age-group
  tmp = aggregate( allDat[,c('ID','Age')], 
                   list( allDat$subject ), unique )
  colnames( tmp ) = c( 'S','ID','Age' )
  parEst = cbind( tmp, parEst )
  
  # Save results as a .csv
  setwd( 'Data' )
  write.table( parEst, file = 'SDT_estimates_14.csv', 
               sep = ',', row.names = F, col.names = T )
  setwd( orig_dir )
  
}

