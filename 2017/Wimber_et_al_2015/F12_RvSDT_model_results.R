#---------------------#
# RvSDT model results #
# Kevin Potter        #
# Updated 01/27/2017  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = 'RvSDT_model.pdf'
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages and data
# Lookup - 02:  Convergence check and initial parameter evaluation
# Lookup - 03:  Posterior retrodictive checks
# Lookup - 04:  Plots of marginal posterior distributions

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

# Load in package for Bayesian estimation
library(rstan)

# Load in useful functions
source( 'F1_useful_functions.R' )

# Load in data
setwd( 'Data' )
load( 'Recog_mem.RData' )
setwd( orig_dir )

# Load in posterior
fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
setwd(fName)
setwd( "Wimber_et_al" )
fname = 'RvSDT_model_post.RData'
load( fname )
setwd( orig_dir )

###
### Convergence check and initial parameter evaluation
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Check convergence of model
  plot_conv( fit, savePlot, 'sigma_lambda' )
  
  if (!savePlot) x11( width=12 )
  pairs( fit, pars = c('beta_dp','beta_c','mu_lambda','Omega') )
  
}

###
### Posterior retrodictive checks
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Extract data
  x = aggregate( d$Ch, list( d$B, d$IT, d$Co, d$S ), mean )
  colnames( x ) = c('B','IT','Co','S','P')
  
  ### Group-level performance
  
  # Generate retrodictive checks for group-level performance
  f = function(s) {
    yhat = rbinom( stanDat$No, 1, post$theta[s,] )
    out = aggregate( yhat, list( d$Cnd, d$Co ), mean )$x
    
    return( out )
  }
  rc = sapply( 1:nrow(post$beta_dp), f )
  
  # Create a plot of the proportions
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  plot( c(.8,2.2), c(0,1), type = 'n', xlab = 'Condition',
        ylab = 'P(Right)', xaxt = 'n', bty = 'n' )
  axis( 1, 1:2, c('Selective retrieval','Baseline'), tick = F )
  abline(h=.5)
  legend('topright','Right correct',bty='n')
  legend('bottomright','Left correct',bty='n')
  legend( 'topleft', c('Target','Competitor'),
          pch = c(19,22), pt.bg = 'white', bty = 'n' )
  
  ya = seq( .1, .9, .01 )
  points( rep(1,length(ya)), ya, pch = '.', cex = .5 )
  points( rep(2,length(ya)), ya, pch = '.', cex = .5 )
  
  pred = apply( rc, 1, findMode )
  
  ui = apply( rc, 1, quantile, prob = c( .025,.16,.25,.75,.84,.975 ) )
  ui = ui[ c(1:3,6:4), ]
  
  xa = c(
    c(.9,2.1),
    c(1.1,1.9),
    c(.9,2.1),
    c(1.1,1.9) )
  clr = c( rgb( .5, .5, .5, .4 ),
           rgb( .5, .5, .5, .4 ),
           rgb( .5, .5, .5, .4 ) )
  sz = c( .05, .04, .03 )
  
  for (i in 1:8) {
    
    for (j in 1:3) {
      
      drawEllipse( sz[j], diff( ui[c(j,j+3),i] ),
                   Xc = xa[i], Yc = ui[j,i] + diff( ui[c(j,j+3),i] )/2,
                   border = NA, col = clr[j] )
      
    }
    
  }
  points( xa, pred, pch = 19, col = 'grey40' )
  
  obs = aggregate( x$P, list( x$B, x$IT, x$Co ), 
                   function(x) c( mean(x), sd(x)/sqrt(length(x)) ) )
  colnames( obs ) = c( 'B', 'IT', 'Co', 'DS' )
  draw_figure( obs$DS[,1] )
  
  mtext( 'Group-level performance', side = 3, outer = T, line = -2 )
  
  ### Subject-level performance
  
  # Generate retrodictive checks for subject-level performance
  indexS = NULL
  f = function(s) {
    yhat = rbinom( stanDat$No, 1, post$theta[s,] )
    out = aggregate( yhat, list( d$Cnd, d$Co, 
                                 d$S ), mean )
    if ( length( indexS ) == 0 ) indexS <<- out[,3]
    
    return( out$x )
  }
  rc = sapply( 1:nrow(post$beta_dp), f )
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:12, 3, 4, byrow = T ) )
  
  for ( s in 1:12 ) {
    blankPlot( c(.8,2.2), c(0,1) )
    title( paste('Subject', s ) )
    if ( s == 1 | s == 5 | s == 9 ) 
      axis( 2, seq(0,1,.25) )
    if ( s == 9 | s == 10 | s == 11 | s == 12 ) 
      axis( 1, 1:2, c('SR','B' ), tick = F )
    
    
    draw_figure( x$P[ x$S == s ] )
    draw_figure( apply( rc[ indexS == s, ], 1, findMode ),
                 clr = c('blue','white'), lnTp = 2 )
  }
  
  if (!savePlot) x11(width=12)
  layout( matrix( 1:12, 3, 4, byrow = T ) )
  
  for ( s in 13:24 ) {
    blankPlot( c(.8,2.2), c(0,1) )
    title( paste('Subject', s ) )
    if ( s == 1 | s == 5 | s == 9 ) 
      axis( 2, seq(0,1,.25) )
    if ( s == 9 | s == 10 | s == 11 | s == 12 ) 
      axis( 1, 1:2, c('SR','B' ), tick = F )
    
    
    draw_figure( x$P[ x$S == s ] )
    draw_figure( apply( rc[ indexS == s, ], 1, findMode ),
                 clr = c('blue','white'), lnTp = 2 )
  }
  
}

###
### Plots of marginal posterior distributions
###
# Lookup - 05

if ( runCode[3] ) {

  if (!savePlot) x11(width=12)
  layout( matrix( 1:4, 1, 4, byrow = T ) )
  
  ttl = c( "Intercept", "SR-C" )
  for ( i in 1:2 ) {
    tmp = hist( post$beta_dp[,i], col = 'grey', border = 'white', 
                breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
    ui = quantile( post$beta_dp[,i], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
  }
  
  # Mean criterion value
  tmp = hist( post$beta_c, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'Bias', xlab = ' ' )
  ui = quantile( post$beta_c, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
  # Mean recall probability
  tmp = hist( post$mu_lambda, col = 'grey', 
              border = 'white', breaks = 40, freq = F, 
              main = 'Recall probability', xlab = ' ',
              xaxt = 'n' )
  axis( 1, logit( seq(.125,.275,.025) ), seq(.125,.275,.025) )
  ui = quantile( post$mu_lambda, c(.025,.975) )
  dn = c( max( which( tmp$mids < ui[1] ) ),
          min( which( tmp$mids > ui[2] ) ) )
  dn = tmp$density[ dn ]
  segments( ui, c(0,0), ui, dn, col = 'blue' )
  
  
}

if (savePlot) dev.off()
setwd( orig_dir )