#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 02/09/2017 #
#--------------------#

# Index
# Lookup - 01:  SDT_calc_4_dist
# Lookup - 02:  SDT_dist_sim
# Lookup - 03:  fancy_violin_plot

# References:
# Hautus, M. J. (1995). Corrections for extreme proportions and their 
#   biasing effects on estimated values of d'. Behavior Research Methods
#   Instruments, & Computers, 27(1), 46 - 51. DOI: 10.3758/BF03203619

# Lookup - 01
SDT_calc_4_dist = function( dat, correct = F ) {
  # Purpose:
  # Calculates parameter estimates for an equal-variance 
  # SDT model with 3 lure distributions and one target 
  # distribution with a fixed criterion (or pair of 
  # criterions) over 4 conditions.
  # Arguments:
  # dat     - A data frame with 13 or 14 rows giving the 
  #           frequencies, hits/false alarm rates, and the 
  #           total number of trials per rate (columns must 
  #           be labeled as 'Y', 'P', and 'N', respectively)
  # correct - If true, applies a correction, using the log-linear 
  #           approach, in which we add .5 to the hits and
  #           false alarm frequencies, then add 1 to the 
  #           total number of trials ( Hautus, 1995 )
  # Returns:
  # A vector of 13 estimated d' values (first value fixed to 0) 
  # and either a single criterion estimate or a pair of estimates
  # depending on the number of observations.
  
  # Check for inappropriate input
  if ( !is.data.frame( dat ) | 
       (nrow(dat) != 13 & nrow(dat) != 14 ) | 
       ( sum( c( 'Y', 'P', 'N' ) %in% colnames( dat ) ) != 3 ) )
    stop( 'Need data frame with 13 or 14 rows and correct variables' )
  
  # Apply desired corrections
  if ( correct ) {
    dat$P = ( dat$Y + .5 )/( dat$N + 1 )
  }
  
  # Check which version of model to estimate
  if ( nrow( dat ) == 13 ) {
    
    # Initialize d' values
    dp = numeric( 13 )
    
    # Calculate criterion from unrelated lures
    k = -qnorm( dat$P[1] )
    
    # Calculate d' values using criterion
    dp[ 2:13 ] = qnorm( dat$P[ 2:13 ] ) + k
    
  } else {
    
    # Initialize d' values
    dp = numeric( 14 )
    
    # Calculate criterion from unrelated lures
    k = c( -qnorm( dat$P[1] ),
           -qnorm( dat$P[2] ) )
    
    ord = rep( 1:2, each = 6 )
    
    # Calculate d' values using criterion
    dp[ 3:14 ] = qnorm( dat$P[ 3:14 ] ) + k[ ord ]
    
  }
  
  # Return estimates
  out = c( dp, k )
  return( out )
}

# Lookup - 02
SDT_4_dist_sim = function( prm, Ntrials = NULL ) {
  # Purpose:
  # A function for simulating responses for a equal-variance 
  # SDT model with 3 lure distributions and one target 
  # distribution with a fixed criterion or a pair of 
  # criterion values over 4 conditions.
  # Arguments:
  # prm     - A vector of 12 d' values and either a single 
  #           criterion or a pair of criterion values
  # Ntrials - An optional vector giving the total number of trials
  #           for each of the 13/14 observations
  # Returns:
  # Either a vector of probabilities or frequencies for picking 
  # "old" over the 4 distributions in each of the 4 conditions.
  
  # Extract parameters
  if ( length( prm ) == 13 ) {
    dp = c( 0, prm[ 1:12 ] )
    k = prm[13]
  } else {
    dp = c( 0, 0, prm[ 1:12 ] )
    k = c( prm[13:14], rep( prm[13:14], each = 6 ) )
  }
  
  # Calculate probability 
  theta = 1 - pnorm( k, dp, 1 )
  
  if ( length( Ntrials ) == 0 ) {
    out = theta
  } else {
    out = rbinom( length( dp ), Ntrials, theta )
  }
  
  return( out )
}

# Lookup - 03
fancy_violin_plot = function( x, pos, scaleH = .4, side = 1,
                              clr = 'black', pts = 19,
                              bclr = 'black', 
                              lnw = 1, lnt = 1 ) {
  # Purpose:
  # Adds a density plot to an already existing figure, 
  # drawing the density heights for each individual 
  # observation on the x-axis.
  # Arguments:
  # x      - A vector of observations
  # pos    - The base of the density on the x-axis 
  # scaleH - The maximum height of the density on the x-axis
  # side   - The side (left = 1, right = -1) to which the 
  #          density should extend
  # clr    - The color of the line and points
  # pts    - The type of point to use (i.e. pch)
  # bclr   - The background color of the point
  # lnw    - The width of the line (i.e. lwd)
  # lnt    - The type of line (i.e.lty)
  
  x = sort( x )
  d = density( x )
  df = approxfun( d )
  xa = df( x ) * scaleH
  
  lines( side * xa + pos, x, col = clr, lwd = lnw, lty = lnt )
  points( side * xa + pos, x, pch = pts, bg = bclr, col = clr )
  
}