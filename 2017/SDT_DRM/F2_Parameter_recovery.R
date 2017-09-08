#--------------------#
# Parameter recovery #
# Kevin Potter       #
# Updated 02/09/2017 #
#--------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicate which code segments to run
runCode = c( T, T )

# Indicate whether a pdf of the figures should be saved
savePlot = T
if ( savePlot ) {
  setwd( 'Plots' )
  pdf( 'Parameter_recovery.pdf', width = 12 )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages, functions, and data
# Lookup - 02:  Parameter recovery (13 conditions)
# Lookup - 03:  Parameter recovery (14 conditions)

###
### Load in useful packages and functions
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
### Parameter recovery (13 conditions)
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Define generating parameters
  gp = c( 0.59, 0.41, 2.17, # P2
          0.53, 0.24, 2.19, # P8
          1.03, 0.64, 2.27, # S2
          0.99, 0.48, 2.26, # S8
          1.29 # Criterion
  )
  
  # Define probability of responding 'old'
  theta = SDT_4_dist_sim( gp );
  
  # Define number of trials per each observation
  Ntrials = ratesDat$N[ ratesDat$S == 1 ]
  Ntrials = c( sum( Ntrials[1:2] ), Ntrials[3:14] )
  
  # Define number of iterations
  Nrep = 1000
  
  # Initialize matrix for difference between 
  # estimates and generating parameters
  diff_est_vs_gen = matrix( NA, Nrep, 14 )
  
  # Initialize data frame for simulated data
  dat = matrix( Ntrials, 13, 3 )
  colnames( dat ) = c( 'N', 'Y', 'P' )
  dat = as.data.frame( dat )
  dat$N = Ntrials
  
  # Loop over iterations
  for ( nr in 1:Nrep ) {
    
    # Simulate observations
    dat$Y = rbinom( 13, Ntrials, theta )
    dat$P = dat$Y/dat$N
    
    # Estimate parameters
    est = SDT_calc_4_dist( dat, correct = T )
    
    diff_est_vs_gen[ nr, ] = est - c(0,gp)
  }
  
  # Create plot of distribution of differences
  if (!savePlot) x11( width = 12 )
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
            pch = 19, 
            col = rgb( .5, .5, .5, .2 ) )
    # Uncertainty intervals
    if ( i > 1 ) {
      
      ui = quantile( diff_est_vs_gen[,i], 
                     prob = c( .025, .16, .25, .5, 
                               .75, .84, .975 ) )
      
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
  
  mtext( 'Parameter recovery (13 conditions)',
         side = 3, outer = T, line = -2 )
  
}

###
### Parameter recovery (14 conditions)
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Define generating parameters
  gp = c( 0.59, 0.41, 2.17, # P2
          0.53, 0.24, 2.19, # P8
          1.03, 0.64, 2.27, # S2
          0.99, 0.48, 2.26, # S8
          1.29, 1.29 # Criterion values
  )
  
  # Define probability of responding 'old'
  theta = SDT_4_dist_sim( gp );
  
  # Define number of trials per each observation
  Ntrials = ratesDat$N[ ratesDat$S == 1 ]
  
  # Define number of iterations
  Nrep = 1000
  
  # Initialize matrix for difference between 
  # estimates and generating parameters
  diff_est_vs_gen = matrix( NA, Nrep, 16 )
  
  # Initialize data frame for simulated data
  dat = matrix( Ntrials, 14, 3 )
  colnames( dat ) = c( 'N', 'Y', 'P' )
  dat = as.data.frame( dat )
  dat$N = Ntrials
  
  # Loop over iterations
  for ( nr in 1:Nrep ) {
    
    # Simulate observations
    dat$Y = rbinom( 14, Ntrials, theta )
    dat$P = dat$Y/dat$N
    
    # Estimate parameters
    est = SDT_calc_4_dist( dat, correct = T )
    
    diff_est_vs_gen[ nr, ] = est - c(0,0,gp)
  }
  
  # Create plot of distribution of differences
  if (!savePlot) x11( width = 12 )
  plot( c(0,17), lowerUpper( .5, as.vector( diff_est_vs_gen ) ),
        type = 'n', xlab = 'Parameter', xaxt = 'n', 
        ylab = 'Estimate - generating', bty = 'n' )
  abline( h = seq(-1.5,1.5,.25), col = 'grey' )
  abline( h = 0 )
  axis( 1, 1:16, c( 'URL_P','URL_S', paste(
    rep( c('RL','CL','T'), 4 ),
    rep( c('_P','_S'), each = 6 ),
    rep( c('2','8','2','8'), each = 3 ), sep = '' ), 'kp', 'ks' ),
    tick = F, cex.axis = .7 )
  
  xa = seq( -.4, .4, length = Nrep )
  for ( i in 1:16 ) {
    points( xa + i, diff_est_vs_gen[,i],
            pch = 19, 
            col = rgb( .5, .5, .5, .2 ) )
    # Uncertainty intervals
    if ( i > 2 ) {
      
      ui = quantile( diff_est_vs_gen[,i], 
                     prob = c( .025, .16, .25, .5, 
                               .75, .84, .975 ) )
      
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
  
  mtext( 'Parameter recovery (14 conditions)',
         side = 3, outer = T, line = -2 )
  
}

if ( savePlot ) dev.off()