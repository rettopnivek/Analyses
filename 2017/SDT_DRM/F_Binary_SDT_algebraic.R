#----------------------------------#
# Algebraic methods for binary SDT #
# Kevin Potter                     #
# Updated 02/06/2017               #
#----------------------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define which code segments to run
runCode = c( F, T )

# Index
# Lookup - 01:  Attach useful packages and define useful functions
# Lookup - 02:  Derivation of algebraic method for d'/criterion 
#               estimation
# Lookup - 03:  Parameter recovery using algebraic method

# References:
# Hautus, M. J. (1995). Corrections for extreme proportions and their 
#   biasing effects on estimated values of d'. Behavior Research Methods
#   Instruments, & Computers, 27(1), 46 - 51. DOI: 10.3758/BF03203619
# Macmillan, N. A. & Kaplan, H. L. (1985). Detection theory analysis 
#   of group data: estimating sensitivity from average hit and 
#   false-alarm rates. Psychological Bulletin, 98(1), 185 - 199.

###
### Attach useful packages and define useful functions
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

SDT_calc_binary = function( dat, centered = F, 
                     correct = 0 ) {
  # Purpose:
  # Calculates estimates for d' and criterion values of  
  # equal-variance SDT for binary data using the algebraic method.
  # Arguments:
  # dat      - Either...
  #            1) A vector of four values, the frequencies for hits 
  #               and false alarms followed by the associated total 
  #               number of trials for each, or...
  #            2) A vector of two values, the proportion of hits and 
  #               false alarms
  # centered - Logical, indicates whether to use the parameterization
  #            in which the distributions are equidistant from 0
  # correct  - The type of correction to use, where...
  #            0 = none
  #            1 = The log-linear approach, add .5 to the hits and
  #                false alarm frequencies, then add 1 to the 
  #                total number of trials ( Hautus, 1995 )
  #            2 = The conditional approach, where only proportions 
  #                equal to 0 or 1 are adjusted by .5/N or (N-.5)/N
  #                respectively, where N is the associated number of 
  #                total trials for the given proportion 
  #                ( Macmillan & Kaplan, 1985 )
  # Returns:
  # A vector giving the estimates for d' and the criterion.
  
  if ( length( dat ) == 4 ) {
    
    # Extract frequencies
    S = dat[1]; N = dat[2];
    # Extract total number of trials
    Ns = dat[3]; Nn = dat[4]
    
    # Determine hits and false-alarm rates
    H = S/Ns; FA = N/Nn;
    
    # Apply corrections if indicated
    if ( correct == 1 ) {
      H = (S+.5)/(Ns+1)
      FA = (N+.5)/(Nn+1)
    }
    if ( correct == 2 ) {
      if ( H == 0 ) H = .5/Ns
      if ( FA == 0 ) FA = .5/Nn
      if ( H == 1 ) H = (Ns-.5)/Ns
      if ( FA == 1 ) FA = (Nn-.5)/Nn
    }
  }
  
  if ( length( dat ) == 2 ) {
    
    # Extract hits and false-alarm rates
    H = dat[1]; FA = dat[2]
    
  }
  
  if ( !centered ) {
    
    # Obtain estimate of d'
    dp_est = qnorm( H ) - qnorm( FA )
    
    # Obtain estimate of criterion
    k_est = qnorm( 1 - FA )
    
  } else {
    
    # Obtain estimate of criterion
    k_est = -.5*( qnorm( H ) + qnorm( FA ) )
    
    # Obtain estimate of d'
    dp_est = 2*( qnorm( H ) + k_est )
    
  }
  
  return( c( dp = dp_est, k = k_est ) )
}

###
### Derivation of algebraic method for d'/criterion estimation
###
# Lookup - 02

if ( runCode[1] ) {
  
  ### For when noise distribution is fixed to have a mean of 0
  
  # Define example parameters to test solutions
  dp = 0.8; k = -.1;
  
  # Determine the associated hits/false alarm rates
  H = 1 - pnorm( k, dp, 1 )
  FA = 1 - pnorm( k, 0.0, 1 )
  
  # Solve a pair of equations simultaneously
  
  # For d'
  
  # 1)
  # 1 - pnorm( k - dp ) == H
  # 1 - pnorm( k ) == FA
  # 2)
  # -pnorm( k - dp ) == H - 1
  # -pnorm( k ) == FA - 1
  # 3)
  # pnorm( k - dp ) == 1 - H
  # pnorm( k ) == 1 - FA
  # 4)
  # k - dp == qnorm( 1 - H )
  # k == qnorm( 1 - FA )
  # 5)
  # k - dp - k == qnorm( 1 - H ) - qnorm( 1 - FA )
  # 6)
  # -dp == qnorm( 1 - H ) - qnorm( 1 - FA )
  # 7)
  # dp == qnorm( 1 - FA ) - qnorm( 1 - H )
  # Equivalent to :
  # dp == qnorm( H ) - qnorm( FA )
  
  # For criterion
  
  # 1)
  # k - dp == qnorm( 1 - H )
  # k == qnorm( 1 - FA )
  # 5)
  # k == qnorm( 1 - H ) + dp
  # k == qnorm( 1 - FA )
  
  ### For when the distribution are set equidistant from 0
  
  # Define example parameters to test solutions
  dp = 0.8; k = -.1;
  
  # Determine the associated hits/false alarm rates
  H = 1 - pnorm( k, dp/2, 1 )
  FA = 1 - pnorm( k, -dp/2, 1 )
  
  # Solve a pair of equations simultaneously
  
  # For criterion
  
  # 1)
  # 1 - pnorm( k - dp/2 ) == H
  # 1 - pnorm( k + dp/2 ) == FA
  # 2)
  # -pnorm( k - dp/2 ) == H - 1
  # -pnorm( k + dp/2 ) == FA - 1
  # 3)
  # pnorm( k - dp/2 ) == 1 - H
  # pnorm( k + dp/2 ) == 1 - FA
  # 4)
  # k - dp/2 == qnorm( 1 - H )
  # k + dp/2 == qnorm( 1 - FA )
  # 5)
  # k - dp/2 == -qnorm( H )
  # k + dp/2 == -qnorm( FA )
  # 6)
  # k - dp/2 + k + dp/2 == -qnorm( H ) - qnorm( FA )
  # 7)
  # 2*k == -qnorm( H ) - qnorm( FA )
  # 8)
  # k == -.5*( qnorm( H ) + qnorm( FA ) )
  
  # For d'
  
  # 1)
  # k - dp/2 == qnorm( 1 - H )
  # 2)
  # -dp/2 == -qnorm( H ) - k
  # 3)
  # dp/2 == qnorm( H ) + k
  # 4)
  # dp == 2*( qnorm( H ) + k )
  
}

###
### Parameter recovery using algebraic method
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Define a set of generating parameters
  gp = c(
    dp = 1.5, # Memory strength
    k = -.3 # Bias
  )
  
  # Define version to use
  centered = T
  
  if ( centered ) {
    # Determine the associated hits/false alarm rates
    H = 1 - pnorm( gp['k'], gp['dp']/2, 1 )
    FA = 1 - pnorm( gp['k'], -gp['dp']/2, 1 )
  } else {
  # Determine the associated hits/false alarm rates
  H = 1 - pnorm( gp['k'], gp['dp'], 1 )
  FA = 1 - pnorm( gp['k'], 0, 1 )
  }
  
  # Simulate data and estimate parameters over 
  # multiple iterations
  Nrep = 1000
  est = matrix( NA, Nrep, 6 )
  colnames( est ) = c('dp','k','dp1','k1','dp2','k2')
  
  for ( nr in 1:Nrep ) {
    
    # Simulate a set of data
    dat = c( O = NA, N = NA, No = 20, Nn = 20 )
    dat['O'] = rbinom( 1, dat['No'], H )
    dat['N'] = rbinom( 1, dat['Nn'], FA )
    
    for ( v in 0:2 ) {
      tmp = SDT_calc_binary( dat, centered = centered,
                             correct = v )
      est[ nr, 1:2 + 2*v ] = tmp
    }
    
  }
  est = as.data.frame( est )
  
  # Plot the difference of the estimates and the generating 
  # parameters
  
  plotYes = T
  if ( plotYes ) {
    
    # Plot results for d'
    x11( width = 12 )
    layout( cbind( 1, 2, 3 ) )
    
    # No correction
    sel = est$dp != Inf & est$dp != -Inf & !is.na( est$dp )
    yl = lowerUpper( .5, est$dp[ sel ] - gp[1] )
    plot( c(0, Nrep + 1), yl, type = 'n', bty = 'l',
          xlab = 'Iteration', ylab = 'Estimate - generating' )
    points( 1:sum(sel), est$dp[sel] - gp[1], pch = 19 )
    abline( h = 0, lwd = 2 )
    abline( h = mean( est$dp[sel] - gp[1] ), 
            col = 'blue', lty = 2, lwd = 2 )
    
    # Log-linear approach
    yl = lowerUpper( .5, est$dp1 - gp[1] )
    plot( c(0, Nrep + 1), yl, type = 'n', bty = 'l',
          xlab = 'Iteration', ylab = 'Estimate - generating' )
    points( 1:Nrep, est$dp1 - gp[1], pch = 19 )
    abline( h = 0, lwd = 2 )
    abline( h = mean( est$dp1 - gp[1] ), 
            col = 'blue', lty = 2, lwd = 2 )

    # Conditional approach
    yl = lowerUpper( .5, est$dp2 - gp[1] )
    plot( c(0, Nrep + 1), yl, type = 'n', bty = 'l',
          xlab = 'Iteration', ylab = 'Estimate - generating' )
    points( 1:Nrep, est$dp2 - gp[1], pch = 19 )
    abline( h = 0, lwd = 2 )
    abline( h = mean( est$dp2 - gp[1] ), 
            col = 'blue', lty = 2, lwd = 2 )
    
  }
  
}

