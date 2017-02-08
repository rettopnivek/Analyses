#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 02/06/2017 #
#--------------------#

# Index
# Lookup - 01:  SDT_calc_4_dist

# References:
# Hautus, M. J. (1995). Corrections for extreme proportions and their 
#   biasing effects on estimated values of d'. Behavior Research Methods
#   Instruments, & Computers, 27(1), 46 - 51. DOI: 10.3758/BF03203619

# Lookup - 01
SDT_calc_4_dist = function( dat, correct = F ) {
  # Purpose:
  # Calculates parameter estimates for an equal-variance 
  # SDT model with 3 lure distributions and one target 
  # distribution with a fixed criterion over 4 conditions.
  # Arguments:
  # dat     - A data frame with
  #
  # correct - If true, applies a correction, the log-linear 
  #           approach, in which we add .5 to the hits and
  #           false alarm frequencies, then add 1 to the 
  #           total number of trials ( Hautus, 1995 )
  # Returns:
  # A vector of 13 estimated d' values (first value fixed to 0) and 
  # a single criterion estimate.
  
  # Check for inappropriate input
  if ( sum( dim( dat ) != c(13,6) ) > 0 ) 
    stop( 'Input must be 13 by 6 data frame' )
  
  # Apply desired corrections
  if ( correct ) {
    dat$P = ( dat$Y + .5 )/( dat$N + 1 )
  }
  
  # Initialize d' values
  dp = numeric( 13 )
  
  # Calculate criterion from unrelated lures
  k = -qnorm( dat$P[1] )
  
  # Calculate d' values using criterion
  dp[ 2:13 ] = qnorm( dat$P[ 2:13 ] ) + k
  
  # Return estimates
  out = c( dp, k )
  return( out )
}