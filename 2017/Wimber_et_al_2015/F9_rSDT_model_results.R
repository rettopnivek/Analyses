#--------------------#
# rSDT model results #
# Kevin Potter       #
# Updated 01/29/2017 #
#--------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Define code segments to run
runCode = c( T, T, T )

# Indicate model fit to plot
model_num = 3

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = paste( 'rSDT_model_M', model_num, '.pdf', sep = '' )
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
fname = paste( 'rSDT_model_M', model_num, '_post.RData', sep = '' )
load( fname )
setwd( orig_dir )

###
### Convergence check and initial parameter evaluation
###
# Lookup - 02

if ( runCode[1] ) {
  
  # Check convergence of model
  plot_conv( fit, savePlot, 'sigma_zeta' )
  
  if (!savePlot) x11( width=12 )
  pairs( fit, pars = c('Mu','tau','sigma_zeta') )
  
  if (!savePlot) x11( width=12 )
  pairs( fit, pars = c('Omega') )
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
  rc = sapply( 1:nrow(post$Mu), f )
  
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
  rc = sapply( 1:nrow(post$Mu), f )
  
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
  
  if ( model_num == 1 ) {
    
    if (!savePlot) x11(width=12)
    layout( matrix( 1:8, 2, 4, byrow = T ) )
    
    ttl = c( "d' for T-SR", "d' for T-B", "d' for C-SR", "d' for C-B" )
    for ( i in 1:4 ) {
      tmp = hist( post$Mu[,i], col = 'grey', border = 'white', 
                  breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
      ui = quantile( post$Mu[,i], c(.025,.975) )
      dn = c( max( which( tmp$mids < ui[1] ) ),
              min( which( tmp$mids > ui[2] ) ) )
      dn = tmp$density[ dn ]
      segments( ui, c(0,0), ui, dn, col = 'blue' )
    }
    
    # Mean criterion value
    tmp = hist( post$Mu[,5], col = 'grey', 
                border = 'white', breaks = 40, freq = F, 
                main = 'Bias', xlab = ' ' )
    ui = quantile( post$Mu[,5], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
    
    # Difference scores for first associates (targets)
    dT = post$Mu[,2] - post$Mu[,1]
    tmp = hist( dT, col = 'grey', 
                border = 'white', breaks = 40, freq = F, 
                main = 'B - SR for targets', xlab = ' ' )
    ui = quantile( dT, c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
    
    # Difference scores for second associates (competitors)
    dC = post$Mu[,4] - post$Mu[,3]
    tmp = hist( dC, col = 'grey', 
                border = 'white', breaks = 40, freq = F, 
                main = 'B - SR for competitors', xlab = ' ' )
    ui = quantile( dC, c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
    
    # Difference of differences
    tmp = hist( dC - dT, col = 'grey', 
                border = 'white', breaks = 40, freq = F, 
                main = 'Difference of differences', xlab = ' ' )
    ui = quantile( dC - dT, c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
    
    ### Subject-level effects (d')
    
    if (!savePlot) x11( width = 12 )
    layout( cbind( 1 ) )
    
    yl = lowerUpper( .5, as.vector( post$eta_dp ) )
    
    N = dim( post$eta_dp )[2]
    blankPlot( c(0, N + 2 ), yl )
    mtext( 'Subjects', side = 1, line = 2.5 )
    axis( 2, seq( yl[1], yl[2], .5 ) )
    mtext( "Subject effects (d')", side = 2, line = 2.5 )
    abline( h = 0, lty = 2 )
    
    ps = c( -.3, -.1, .1, .3 )
    clr = c('green','orange','blue','red')
    
    md = apply( post$eta_dp[,,3], 2, findMode )
    ord = order( md )
    
    for ( i in 1:4 ) {
      
      # Determine modes and rank-order 
      md = apply( post$eta_dp[,,i], 2, findMode )
      
      points( 1:length( md ) + ps[i], md[ ord ], 
              col = clr[i], pch = 19 )
      
    }
    
    legend( 'topleft', c('T-SR','T-B','C-SR','C-B'),
            fill = clr, bty = 'n' )
    
    ### Item-level effects (d')
    
    # Determine modes and rank-order 
    md = apply( post$zeta, 2, findMode )
    # Uncertainty intervals
    ui = apply( post$zeta, 2, quantile, prob = c(.025,.975) )
    
    # Determine number of times items appeared in different conditions
    res = aggregate( cbind( d$IT == 1 & d$B == 0, 
                            d$IT == 1 & d$B == 1,
                            d$IT == 2 & d$B == 0,
                            d$IT == 2 & d$B == 1 ), list( d$I, d$Cat ), sum )
    colnames( res ) = c('I','Cat','T_SR','T_B','C_SR','C_B')
    
    if (!savePlot) x11( width = 12 )
    layout( cbind( 1 ) )
    
    yl = lowerUpper( .5, as.vector( post$zeta ) )
    
    blankPlot( c(0, stanDat$Ni/3 + 1 ), yl )
    mtext( 'Images', side = 1, line = 2.5 )
    axis( 2, seq( yl[1], yl[2], .5 ) )
    mtext( "Item effects (criterion)", side = 2, line = 2.5 )
    abline( h = 0, lty = 2 )
    
    ps = c( -.3, 0, .3 )
    clr = c( 'black', 'red', 'blue' )
    for ( i in 1:3 ) {
      
      ya = list(
        md[ res$I[ res$Cat == i ] ],
        ui[ ,res$I[ res$Cat == i ] ] )
      ord = order( ya[[1]] )
      points( 1:length( ya[[1]] ) + ps[i], ya[[1]][ ord ], 
              col = clr[i], pch = 19 )
      
    }
    
    legend( 'topright', c('Faces','Objects','Scenes'),
            fill = clr, bty = 'n' )
  }
  
  if ( model_num == 2 ) {
    
    if (!savePlot) x11(width=12)
    layout( matrix( 1:4, 1, 4, byrow = T ) )
    
    ttl = c( "d' for Intercept", "d' for order", "d' for C-SR" )
    for ( i in 1:3 ) {
      tmp = hist( post$Mu[,i], col = 'grey', border = 'white', 
                  breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
      ui = quantile( post$Mu[,i], c(.025,.975) )
      dn = c( max( which( tmp$mids < ui[1] ) ),
              min( which( tmp$mids > ui[2] ) ) )
      dn = tmp$density[ dn ]
      segments( ui, c(0,0), ui, dn, col = 'blue' )
    }
    
    # Mean criterion value
    tmp = hist( post$Mu[,4], col = 'grey', 
                border = 'white', breaks = 40, freq = F, 
                main = 'Bias', xlab = ' ' )
    ui = quantile( post$Mu[,4], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
    
    ### Subject-level effects (d')
    
    if (!savePlot) x11( width = 12 )
    layout( cbind( 1 ) )
    
    yl = lowerUpper( .5, as.vector( post$eta_dp ) )
    
    N = dim( post$eta_dp )[2]
    blankPlot( c(0, N + 2 ), yl )
    mtext( 'Subjects', side = 1, line = 2.5 )
    axis( 2, seq( yl[1], yl[2], .5 ) )
    mtext( "Subject effects (d')", side = 2, line = 2.5 )
    abline( h = 0, lty = 2 )
    
    ps = c( -.3, -.1, .1, .3 )
    clr = c('green','orange','blue')
    
    md = apply( post$eta_dp[,,3], 2, findMode )
    ord = order( md )
    
    for ( i in 1:3 ) {
      
      # Determine modes and rank-order 
      md = apply( post$eta_dp[,,i], 2, findMode )
      
      points( 1:length( md ) + ps[i], md[ ord ], 
              col = clr[i], pch = 19 )
      
    }
    
    legend( 'topleft', c('I','O','C-SR'),
            fill = clr, bty = 'n' )
    
    ### Item-level effects (d')
    
    # Determine modes and rank-order 
    md = apply( post$zeta, 2, findMode )
    # Uncertainty intervals
    ui = apply( post$zeta, 2, quantile, prob = c(.025,.975) )
    
    # Determine number of times items appeared in different conditions
    res = aggregate( cbind( d$IT == 1 & d$B == 0, 
                            d$IT == 1 & d$B == 1,
                            d$IT == 2 & d$B == 0,
                            d$IT == 2 & d$B == 1 ), list( d$I, d$Cat ), sum )
    colnames( res ) = c('I','Cat','T_SR','T_B','C_SR','C_B')
    
    if (!savePlot) x11( width = 12 )
    layout( cbind( 1 ) )
    
    yl = lowerUpper( .5, as.vector( post$zeta ) )
    
    blankPlot( c(0, stanDat$Ni/3 + 1 ), yl )
    mtext( 'Images', side = 1, line = 2.5 )
    axis( 2, seq( yl[1], yl[2], .5 ) )
    mtext( "Item effects (criterion)", side = 2, line = 2.5 )
    abline( h = 0, lty = 2 )
    
    ps = c( -.3, 0, .3 )
    clr = c( 'black', 'red', 'blue' )
    for ( i in 1:3 ) {
      
      ya = list(
        md[ res$I[ res$Cat == i ] ],
        ui[ ,res$I[ res$Cat == i ] ] )
      ord = order( ya[[1]] )
      points( 1:length( ya[[1]] ) + ps[i], ya[[1]][ ord ], 
              col = clr[i], pch = 19 )
      
    }
    
    legend( 'topright', c('Faces','Objects','Scenes'),
            fill = clr, bty = 'n' )
    
  }
  
  if ( model_num == 3 ) {
    
    if (!savePlot) x11(width=12)
    layout( matrix( 1:4, 1, 4, byrow = T ) )
    
    ttl = c( "d' for Intercept", "d' for order" )
    for ( i in 1:2 ) {
      tmp = hist( post$Mu[,i], col = 'grey', border = 'white', 
                  breaks = 40, freq = F, main = ttl[i], xlab = ' ' )
      ui = quantile( post$Mu[,i], c(.025,.975) )
      dn = c( max( which( tmp$mids < ui[1] ) ),
              min( which( tmp$mids > ui[2] ) ) )
      dn = tmp$density[ dn ]
      segments( ui, c(0,0), ui, dn, col = 'blue' )
    }
    
    # Mean criterion value
    tmp = hist( post$Mu[,3], col = 'grey', 
                border = 'white', breaks = 40, freq = F, 
                main = 'Bias', xlab = ' ' )
    ui = quantile( post$Mu[,3], c(.025,.975) )
    dn = c( max( which( tmp$mids < ui[1] ) ),
            min( which( tmp$mids > ui[2] ) ) )
    dn = tmp$density[ dn ]
    segments( ui, c(0,0), ui, dn, col = 'blue' )
    
    blankPlot()
    
    ### Subject-level effects (d')
    
    if (!savePlot) x11( width = 12 )
    layout( cbind( 1 ) )
    
    yl = lowerUpper( .5, as.vector( post$eta_dp ) )
    
    N = dim( post$eta_dp )[2]
    blankPlot( c(0, N + 2 ), yl )
    mtext( 'Subjects', side = 1, line = 2.5 )
    axis( 2, seq( yl[1], yl[2], .5 ) )
    mtext( "Subject effects (d')", side = 2, line = 2.5 )
    abline( h = 0, lty = 2 )
    
    ps = c( -.3, -.1, .1, .3 )
    clr = c('green','orange')
    
    md = apply( post$eta_dp[,,2], 2, findMode )
    ord = order( md )
    
    for ( i in 1:2 ) {
      
      # Determine modes and rank-order 
      md = apply( post$eta_dp[,,i], 2, findMode )
      
      points( 1:length( md ) + ps[i], md[ ord ], 
              col = clr[i], pch = 19 )
      
    }
    
    legend( 'topleft', c('I','O','C-SR'),
            fill = clr, bty = 'n' )
    
    ### Item-level effects (d')
    
    # Determine modes and rank-order 
    md = apply( post$zeta, 2, findMode )
    # Uncertainty intervals
    ui = apply( post$zeta, 2, quantile, prob = c(.025,.975) )
    
    # Determine number of times items appeared in different conditions
    res = aggregate( cbind( d$IT == 1 & d$B == 0, 
                            d$IT == 1 & d$B == 1,
                            d$IT == 2 & d$B == 0,
                            d$IT == 2 & d$B == 1 ), list( d$I, d$Cat ), sum )
    colnames( res ) = c('I','Cat','T_SR','T_B','C_SR','C_B')
    
    if (!savePlot) x11( width = 12 )
    layout( cbind( 1 ) )
    
    yl = lowerUpper( .5, as.vector( post$zeta ) )
    
    blankPlot( c(0, stanDat$Ni/3 + 1 ), yl )
    mtext( 'Images', side = 1, line = 2.5 )
    axis( 2, seq( yl[1], yl[2], .5 ) )
    mtext( "Item effects (criterion)", side = 2, line = 2.5 )
    abline( h = 0, lty = 2 )
    
    ps = c( -.3, 0, .3 )
    clr = c( 'black', 'red', 'blue' )
    for ( i in 1:3 ) {
      
      ya = list(
        md[ res$I[ res$Cat == i ] ],
        ui[ ,res$I[ res$Cat == i ] ] )
      ord = order( ya[[1]] )
      points( 1:length( ya[[1]] ) + ps[i], ya[[1]][ ord ], 
              col = clr[i], pch = 19 )
      
    }
    
    legend( 'topright', c('Faces','Objects','Scenes'),
            fill = clr, bty = 'n' )
    
  }
  
}

if (savePlot) dev.off()
setwd( orig_dir )