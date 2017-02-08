#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 02/07/2017 #
#--------------------#

# Lookup - 01
S = function(x) {
  # Purpose:
  # Calculates descriptive statistics for a sample of values
  # Arguments:
  # x - A vector of values
  # Returns:
  # The mean, standard error of the mean, the inter-quartile range,
  # and the min/max of the sample.
  
  out = c( mean(x), sem(x), diff( quantile(x,prob=c(.25,.75)) ), 
           min(x), max(x) )
  
  return( out )
}

# Lookup - 02
violinPlot = function( x, pos, crit = NULL, scaleH = .5, 
                       type = 'greater', ... ) {
  # Purpose:
  # Adds a violin distribution to an already existing 
  # plot.
  # Arguments:
  # x      - A vector of numerical values
  # pos    - The position on the x-axis to draw the violin
  #          distribution
  # crit   - The value past which to calculate the posterior 
  #          probability.
  # scaleH - The width of the distribution
  # type   - Whether to calculate the probability 'greater' than 
  #          or 'less' than the critical value
  
  den = density( x )
  
  if ( length(crit) > 0 ) {
    if (type=='greater') {
      sel = den$x > crit
      PostProb = sum( x > crit )/length(x)
    }
    if (type=='less') {
      sel = den$x < crit
      PostProb = sum( x < crit )/length(x)
    }
  } else {
    sel = rep( T, length( den$x ) )
    PostProb = 1
  }
  
  den$y = den$y/max(den$y); den$y = den$y*scaleH;
  xa = c( -den$y[sel], rev(den$y[sel]) ) + pos
  ya = c( den$x[sel], rev(den$x[sel]) )
  polygon( xa, ya, ... )
  
  return( PostProb )
}