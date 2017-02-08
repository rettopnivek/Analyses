#--------------------#
# WR model check     #
# Kevin Potter       #
# Updated 02/04/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, F )

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Load in useful packages and data
# Lookup - 03:  Determine generating values/priors
# Lookup - 04:  Simulation example (Hierarchical)

###
### Load in useful packages and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/rtplots")
library(rtplots)

# Package for sequential sampling models
# install_github("rettopnivek/seqmodels")
library(seqmodels)

# Load in package for Bayesian estimation
library(rstan)
# For parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load in useful functions
source( 'F1_useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'Recog_mem.RData' )
setwd( orig_dir )

###
### Pre-process data 
###
# Lookup - 02

# Remove missing responses
sel = which( is.na( d$RT ) )
d = d[ -sel, ]

# Determine sample size
N = length( unique( d$S ) )

# Remove overly fast responses
sel = which( d$RT < .2 )
d = d[ -sel, ]

###
### Determine generating values/priors
###
# Lookup - 03

if ( runCode[1] ) {
  
  # Determine average range
  xRange = c(
    mean( aggregate( d$RT, list( d$S ), min )$x ),
    mean( aggregate( d$RT, list( d$S ), max )$x ) )
  
  # Extract average CDF
  cdf1 = cdf_curve( d$RT, d$Ac, grp = d$S, 
             opt = list( draw = F, out = T ),
             xRange = xRange )
  cdf0 = cdf_curve( d$RT, d$Ac, sel = 0, grp = d$S, 
                    opt = list( draw = F, out = T ),
                    xRange = xRange )
  # Generate plot of CDF curves
  # x11();
  # blankRTplot( tDim = c( .4, 2.4) )
  # lines( cdf1$pv$x, cdf1$pv$y )
  # lines( cdf0$pv$x, cdf0$pv$y, lty = 2 )
  
  # Generate data from a meta-subject
  No = 1000
  Ysim = matrix( NA, No, 2 )
  Ysim[,2] = rbinom( No, 1, mean( d$Ac ) )
  u = runif( No )
  
  cdf1 = list( x = c(0,cdf1$pv$x),
               y = c(0,cdf1$pv$y/mean(d$Ac == 1 ) ) )
  cdf1$y[ length( cdf1$y ) ] = 1
  cdf0 = list( x = c(0,cdf0$pv$x),
               y = c(0,cdf0$pv$y/mean(d$Ac == 0 ) ) )
  cdf0$y[ length( cdf0$y ) ] = 1
  
  # Define a function for linear interpolation
  lin_interp = function( yN, y0, y1, x0, x1 ) {
    
    b1 = ( y1 - y0 ) / ( x1 - x0 ); # Slope
    b0 = y1 - b1*x1; # Intercept
    num = yN - b0;
    xN = ( num )/b1;
    
    return( xN );
  }
  
  # Generate RTs via inverse probability approximation
  for ( no in 1:No ) {
    
    if ( Ysim[no,2] == 1 ) {
      sel = max( which( cdf1$y < u[no] ) )
      Ysim[no,1] = lin_interp( u[no], cdf1$y[sel], cdf1$y[sel+1], 
                  cdf1$x[sel], cdf1$x[sel+1] )
    }
    
    if ( Ysim[no,2] == 0 ) {
      sel = max( which( cdf0$y < u[no] ) )
      Ysim[no,1] = lin_interp( u[no], cdf0$y[sel], cdf0$y[sel+1], 
                               cdf0$x[sel], cdf0$x[sel+1] )
    }
    
  }
  
  # Trim RTs less than minimum RT
  Ysim = Ysim[ Ysim[,1] > xRange[1], ]
  
  # Fit the Wald race model to the simulated data of the 
  # meta-subject
  mle_f = function( prm, dat, priors = NULL ) {
    
    prm[1:4] = exp( prm[1:4] )
    prm[5] = logistic( prm[5] ) * min( dat[,1] )
    
    sll = sum( 
      dwaldrace( dat[,1], dat[,2], prm[1], prm[2], prm[5],
                 prm[3], prm[4], prm[5], ln = 1 ) )
    if ( is.na( sll ) ) sll = -Inf
    
    return( sll )
  }
  st_f = function() {
    
    out = runif( 5, c(.5,.5,.5,.5,-1),
                 c(2,4,2,4,2) )
    out[1:4] = log( out[1:4] )
    
    return( out )
  }
  
  res = MLE( Ysim, mle_f, st_f )
  
  param = c( exp( res$param[ 1:2 ] ),
             logistic( res$param[5] ) * min( Ysim[,1] ),
             exp( res$param[ 3:4 ] ),
             logistic( res$param[5] ) * min( Ysim[,1] ) )
  
  plotYes = T
  if ( plotYes ) {
    
    # Determine RT values to plot
    rn = range( d$RT )
    xa = seq( rn[1], rn[2], length = 20 )
    
    # Extract empirical densities for each subject
    ya = array( NA, dim = c( N, length( xa ), 2 ) )
    
    # Loop over subjects
    for ( s in 1:N ) {
      
      # Estimate joint densities
      inc = 1
      for ( ac in 1:0 ) {
        
        t = d$RT[ d$S == s & d$Ac == ac ]
        p = mean( d$Ac[ d$S == s ] == ac )
        if ( length( t ) > 0 ) {
          ed = density( t )
          af = approxfun( ed )
          ya[ s, , inc ] = af( xa )*p
        }
        inc = inc + 1
      }
      
    }
    
    yl = lowerUpper( .5, na.omit( as.vector( ya ) ) )
    yl = c( -yl[2], yl[2] )
    
    # Generate plot of densities with overlaid average
    x11()
    clr = rainbow( N )
    blankRTplot( tDim = c(0,2.5), pDim = yl, ver = 'PDF' )
    
    # Add in individual densities
    for ( s in 1:N ) {
      lines( xa, ya[s,,1], col = clr[s] )
    }
    
    for ( s in 1:N ) {
      lines( xa, -ya[s,,2], col = clr[s], lty = 2 )
    }
    
    # Overlay averages
    ya_avg = numeric( length( xa ) )
    for ( i in 1:length( xa ) ) {
      ya_avg[i] = mean( na.omit( ya[,i,1] ) )
    }
    lines( xa, ya_avg, col = 'black', lwd = 2 )
    
    for ( i in 1:length( xa ) ) {
      ya_avg[i] = mean( na.omit( ya[,i,2] ) )
    }
    lines( xa, -ya_avg, col = 'black', lwd = 2 )
    
    dist_plot( param, 
               ver = 'PDF', col = 'blue', lty = 2, lwd = 2 )
    dist_plot( param, ch = 0, 
               ver = 'PDF', col = 'blue', lty = 2, lwd = 2,
               opt = list( flip = T ) )
    
  }
  
  # Maximum likelihood estimates from fits to 
  # meta-subject:
  
  # k1 = 3.36
  # x1 = 2.79
  # k0 = 3.59
  # x0 = 1.85
  # tau = .07
  
}

###
### Simulation example (Hierarchical)
###
# Lookup - 04

if ( runCode[2] ) {
  
  # Define generating parameters
  # Intercept 
  beta_x = c( .842, .5, 0, 0, 0 )
  beta_c = c( 1.245, 0 )
  sigma_zeta = .25
  Omega = rbind( c(1, -.3, .3), c(-.3, 1, -.3),
                 c(.3, -.3, 1) )
  tau = c(.71,.71,.71)
  
}


