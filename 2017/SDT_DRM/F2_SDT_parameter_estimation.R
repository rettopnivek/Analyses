#--------------------------#
# SDT parameter estimation #
# Kevin Potter             #
# Updated 02/06/2017       #
#--------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicate which code segments to run
runCode = c( T, T )

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
### Simulation example
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Define generating parameters
  gp = c( 0.00, # URL
          0.59, 0.41, 2.17, # P2
          0.53, 0.24, 2.19, # P8
          1.03, 0.64, 2.27, # S2
          0.99, 0.48, 2.26, # S8
          1.29 # Criterion
  )
  
  # Define probability of responding 'old'
  theta = 1 - pnorm( gp[14], gp[1:13] )
  
  # Define number of trials per each observation
  Ntrials = ratesDat$N[ ratesDat$S == 1 ]
  
  # Define number of iterations
  Nrep = 1000
  
  # Initialize matrix for difference between 
  # estimates and generating parameters
  diff_est_vs_gen = matrix( NA, Nrep, 14 )
  
  # Initialize data frame for simulated data
  dat = matrix( 1:13, 13, 6 )
  colnames( dat ) = c( 'Cond', 'S', 'N', 'Y', 'P', 'A' )
  dat = as.data.frame( dat )
  dat$S = 1; dat$A = 0; dat$N = Ntrials
  
  # Loop over iterations
  for ( nr in 1:Nrep ) {
    
    # Simulate observations
    dat$Y = rbinom( 13, Ntrials, theta )
    dat$P = dat$Y/dat$N
    
    # Estimate parameters
    est = SDT_calc_4_dist( dat, correct = T )
    
    diff_est_vs_gen[ nr, ] = est - gp
  }
  
  # Create plot of distribution of differences
  plotYes = T
  if ( plotYes ) {
    
    x11( width = 12 )
    plot( c(0,15), lowerUpper( .5, as.vector( diff_est_vs_gen ) ),
          type = 'n', xlab = 'Parameter', xaxt = 'n', 
          ylab = 'Estimate - generating', bty = 'n' )
    abline( h = seq(-1.5,1.5,.25), col = 'grey' )
    abline( h = 0 )
    axis( 1, 1:14, c( 'URL', paste(
      rep( c('RL','CL','T'), 4 ),
      rep( c('_P','_S'), each = 6 ),
      rep( c('2','8','2','8'), each = 3 ), sep = '' ), 'k' ),
      tick = F )
    
    xa = seq( -.4, .4, length = Nrep )
    for ( i in 1:14 ) {
      points( xa + i, diff_est_vs_gen[,i],
                              pch = 19, col = rgb( .5, .5, .5, .2 ) )
      # Uncertainty intervals
      if ( i > 1 ) {
        
        ui = quantile( diff_est_vs_gen[,i], 
                       prob = c( .025, .16, .25, .5, .75, .84, .975 ) )
        
        segments( rep( -.4, 4 ) + i, ui[ 1:4 ],
                  rep( .4, 4 ) + i, ui[ 1:4 ],
                  col = c( 'blue', 'orange', 'purple', 'green' ),
                  lwd = 2 )
        segments( rep( -.4, 3 ) + i, ui[ 7:5 ],
                  rep( .4, 3 ) + i, ui[ 7:5 ],
                  col = c( 'blue', 'orange', 'purple' ),
                  lwd = 2 )
      }
    }
    
    legend( 'topleft', c( '95%', '68%','50%', 'Median'),
            fill = c( 'blue', 'orange', 'purple', 'green' ),
            bty = 'n', horiz = T )
    
    mtext( 'Parameter recovery', side = 3, outer = T, line = -2 )
    
  }
  
}

###
### Estimates from data
###
# Lookup - 02

if ( runCode[2] ) {
  
  # Initialize a matrix for the estimates
  parEst = matrix( NA, N, 14 )
  
  # Loop over subjects and calculate parameters
  for ( s in 1:N ) {
    parEst[ s, ] = SDT_calc_4_dist( ratesDat[ ratesDat$S == s, ],
                                     correct = T )
  }
  
}

