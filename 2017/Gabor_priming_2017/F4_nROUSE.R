#-----------------------------#
# nROUSE parameter estimation #
# Kevin Potter                #
# Updated 11/08/2017          #
#-----------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicator if data from Experiment 1 or 2 should be loaded in
Exp = 1

# Indicate whether to create a PDF of figures
savePlot = F

# Indicate which code segments to run
runCode = c(
  T,  # 
  F,  # ?
  F,  # ?
  F
)

# Table of contents for figures
table_of_contents = c(
  '???'
)

# Index
# Lookup - 01:  Load in data, useful packages, and functions
# Lookup - 02:  
# Lookup - 03:  
# Lookup - 04:  

###
###
###
#


### Load in data

if ( Exp == 1 ) {
  dName = 'Gabor_Exp_1_11-17.RData'
} else {
  dName = 'Gabor_Exp_2_11-17.RData'
}

setwd( 'Data' )
load( dName )
rm( dName )
setwd( orig_dir )

# For easy variable access
rd = rawDat
colnames(rd) = c( 'S', 'RT', 'Ch', 'Ac', 'Co', 'PT', 
                  'PD', 'FC', 'TC', 'PP', 'A', 'BT', 
                  'PTL', 'CoL', 'ChL', 'PS', 'PSL', 
                  'BN', 'Cnd', 'DxP', 'ID', 'E' )

### Load in useful packages

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Functions for the nROUSE model
# install_github("rettopnivek/nROUSE")
library(nROUSE)

# Functions for optimization
# install.packages( 'optimx' )
library( optimx )

### Define useful functions

source( 'F1_useful_functions.R' )

# Define a function for MLE
mle_f = function( prm, dat, sum = T, predict = F, priors = NULL ) { 
  
  # Calculate the log-likelihoods 
  mp = c( 2, 6, 9:11 )
  p_index = NULL
  if ( any( dat$p_index != 0 ) ) {
    p_index = dat$p_index[ dat$p_index != 0 ]
  }
  
  if ( is.null( p_index ) ) { prm = 1; mp = NULL } else 
    mp = mp[ p_index ]
  
  if ( predict ) {
    theta = nROUSE_logLik( prm, dat, mapping = mp, predict = T ) 
    return( theta )
  }
  
  ll = nROUSE_logLik( prm, dat, mapping = mp, estimate = F ) 
  if ( !sum ) return( ll )
  
  # Sum the log-likelihoods 
  sll = sum( ll ) 
  
  # Check for NA values 
  if ( is.na( sll ) ) sll = -Inf 
  
  return( sll ) 
}

# Define a function to generate starting values
st_f = function() {
  
  sv = c( N = .0302, I = .9844, 
          SV = .0294, SO = .0609, SS = .015 )
  sv = log( sv )
  
  sv = sv + runif( 5, -3, 3 )
  if ( is.null( nDat$p_index ) ) out = NULL else 
    out = sv[nDat$p_index]
  return( out )
}

### Initialize plot

if ( savePlot ) {
  setwd( 'Plots' )
  pdf( 'nROUSE_results.pdf', width = 12 )
  setwd( orig_dir )
}

# Create table of contents
if ( !savePlot ) x11( width = 12 )
tableContents( table_of_contents, txtSz = 2 )

###
###
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Consider main trials only
  sel = rd$BT == 2
  cd = rd[ sel, ]
  
  # Determine aggregate performance
  nDat = aggregate( cd$Ac, list( cd$PD, cd$PTL, cd$PT ), mean )
  
  # Incorporate prime duration and type
  colnames( nDat ) = c( 'PrimeDur', 'PTL', 'PT', 'P' )
  nDat$Type = 0;
  nDat$Type[ nDat$PT == 1 ] = 1;
  nDat$Type[ nDat$PT == 2 ] = -1;
  
  # Incorporate target and mask durations
  nDat$TarDur = 84
  nDat$MaskDur = 1
  
  # Compute the frequency correct
  tmp = aggregate( rep(1,nrow(cd)), 
                   list( cd$PD, cd$PT[sel], cd$S[sel] ), 
                   sum )
  nDat$N = round( aggregate( tmp$x, list( tmp[,1], tmp[,2] ), mean )$x )
  nDat$Y = round( nDat$P * nDat$N )
  
  # Define a scaling adjustment parameter
  pscl = c( 3.571, 0.231, 0.166 )
  # Define index for parameters
  nDat$p_index = 0; nDat$p_index[1:5] = 1:5
  
  # Clean up workspace
  rm( tmp )
  
  stp = 10
  inc = 1
  chk = -Inf
  while( chk == -Inf & inc < stp ) {
    sv = st_f()
    chk = mle_f(sv,nDat)
  }
  
  sv = log( c( N = .025, I = .5, 
          SV = .01, SO = .01, SS = .01 ) )
  
  test = optim( sv, mle_f, dat = nDat, 
                 control = list( fnscale = -1, maxit = 10000 ) )
  
}


setwd( orig_dir )